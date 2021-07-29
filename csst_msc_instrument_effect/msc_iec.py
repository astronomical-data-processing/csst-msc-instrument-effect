import os
from typing import List
from astropy.io.fits import header
from datetime import datetime

import numpy as np
from astropy.io import fits

from csst_msc_instrument_effect.msc_crmask import CRMask
from csst_dfs_api.facility.level0 import Level0DataApi
from csst_dfs_api.facility.level0prc import Level0PrcApi


class InstrumentEffectCorrection:
    def __init__(self, data_path, bias_path, dark_path, flat_path, output_path, cray_path) -> None:
        self.data_path = data_path
        self.bias_path = bias_path
        self.dark_path = dark_path
        self.flat_path = flat_path
        self.cray_path = "/home/csstpipeline/data/bkg/MSC_CRD_210525121000_100000000_08_raw.fits"
        self.output = output_path
        # RDNOISE
        self.RDNOISE = "RDNOISE1"
        # GAIN
        self.GAIN = "GAIN1"
        # EXPTIME
        self.EXPTIME = "EXPTIME"

    def set_data(self):
        # data
        with fits.open(self.data_path) as hdul:
            self.primary_header = hdul[0].header
            self.data_raw = hdul[1].data.astype(int)
            self.image_header = hdul[1].header

    def set_bias(self):
        # bias
        bias = []
        paths = self.bias_path if isinstance(self.bias_path, list) else [self.bias_path]
        for path in paths:
            du = fits.getdata(path)
            du = du.astype(int)
            bias.append(du)
        self.bias = np.median(bias, axis=0)

    def set_dark(self):
        # dark
        dark = []
        paths = self.dark_path if isinstance(self.dark_path, list) else [self.dark_path]
        for path in paths:
            with fits.open(path) as hdul:
                du = hdul[1].data
                hu = hdul[0].header
                du = du.astype(int)
                du = du - self.bias
                du = du - self.cray / self.image_header[self.GAIN]
                du = du / hu[self.EXPTIME] * \
                    self.primary_header[self.EXPTIME]
                dark.append(du)
        self.dark = np.median(dark, axis=0)

    def set_flat(self):
        # flat
        flat = []
        paths = self.flat_path if isinstance(self.flat_path, list) else [self.flat_path]
        for path in paths:
            with fits.open(path) as hdul:
                du = hdul[1].data
                hu = hdul[0].header
                du = du.astype(int)
                du = du - self.bias
                du = du / hu[self.EXPTIME] * \
                    self.primary_header[self.EXPTIME]
                du = du - self.dark
                du = du / np.median(du)
                flat.append(du)
        self.flat = np.median(flat, axis=0)

    def set_cray(self):
        # cray
        self.cray = fits.getdata(self.cray_path).astype(int)

    def fix_data(self,):
        self.data_fix0 = (self.data_raw - self.bias - self.dark) / self.flat

    def set_flag(self,):
        flag = np.zeros_like(self.data_raw, dtype=np.uint16)
        # 00000001:   坏像元    因为探测器本身原因造成不可用于有效科学研究的像元, 像元响应小于0.5中值 大于1.5中值
        med = np.median(self.flat)
        flg = (self.flat < 0.5 * med) | (1.5 * med < self.flat)
        flag = flag | (flg * 1)
        # 00000010:   热像元    因为探测器本身原因造成的影响科学研究结果的像元. 像元在150秒积分时间内, 暗流计数大于探测器平均读出噪声的平方.
        dark = self.dark.copy()
        dark[dark < 0] = 0
        flg = 1 * self.image_header[self.RDNOISE] ** 2 <= dark  # 不确定是否包含 暂定包含
        flag = flag | (flg * 2)
        # 00000100:   暖像元    因为探测器本身原因造成的影响科学研究结果的像元. 像元的150秒积分时间内, 暗流计数大于0.5被读出噪声的平方, 但小于1被读出噪声的平方.
        flg = (0.5 * self.image_header[self.RDNOISE] ** 2 < dark) & (
            dark < 1 * self.image_header[self.RDNOISE] ** 2
        )
        flag = flag | (flg * 4)
        # 00001000:   饱和溢出像元  饱和像元及流量溢出污染的像元.
        flg = self.data_raw == 65535
        flag = flag | (flg * 8)
        del dark
        del flg
        del med
        # 00010000:   宇宙线像元    宇宙线污染的像元
        crobj = CRMask(self.data_fix0)
        flag_fits, data_fits = crobj.cr_mask()
        flag = flag | (flag_fits[1].data * 16)
        # 00100000:   卫星或者人造移动天体轨迹污染的像元.
        pass
        # 01000000:   鬼像污染的像元.
        pass
        # 10000000:   散射光污染的像元(包括dragon breath)
        pass
        self.flag = flag
        # self.data_fix1 = data_fits[1].data

    def set_weight(self):
        data = self.data_fix0.copy()
        data[self.data_fix0 < 0] = 0
        weight = 1 / (
            self.image_header[self.GAIN] * data +
            self.image_header[self.RDNOISE] ** 2
        )
        weight[self.flag > 0] = 0
        self.weight = weight

    def save(self):
        data_filename = self.data_path
        data_basename = os.path.basename(data_filename).replace("raw", "img")
        data_output = os.path.join(self.output, data_basename)
        data_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.data_fix0.astype(np.float32)
                    / self.primary_header[self.EXPTIME],
                    header=self.image_header,
                ),
            ]
        )
        data_fits.writeto(data_output)
        self.data_output = data_output

        flag_output = data_output.replace("img", "flg")
        flag_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.flag.astype(np.uint16), header=self.image_header
                ),
            ]
        )
        flag_fits.writeto(flag_output)
        self.flag_output = flag_output

        weight_output = data_output.replace("img", "wht")
        weight_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.weight.astype(np.float32), header=self.image_header
                ),
            ]
        )
        weight_fits.writeto(weight_output)
        self.weight_output = weight_output

        header_name = data_basename.replace('.fits', '.head')
        header_output = os.path.join(self.output, header_name)
        hu = self.image_header
        hu['INSTRU_S'] = (0, 'instrument effect status')
        hu['INST_TOL'] = (1234, 'instrument effect operation time')
        hu['INSTRU_V'] = ('0.1', 'instrument effect version')
        hu['INSTRU_P'] = ('test.conf', 'instrument effect config file')
        header_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=hu)
            ]
        )
        header_fits.writeto(header_output)
        self.header_output = header_output

    def run(self):
        self.set_data()
        self.set_cray()
        self.set_bias()
        self.set_dark()
        self.set_flat()
        self.fix_data()
        self.set_flag()
        self.set_weight()
        self.save()

if __name__ == "__main__":
    iec = InstrumentEffectCorrection(
        data_path="/home/csstpipeline/data/L0/MSC_MS_210525121500_100000001_08_raw.fits",
        bias_path=[
            "/home/csstpipeline/data/bkg/MSC_CLB_210525120000_100000000_08_raw.fits"],
        dark_path=[
            "/home/csstpipeline/data/bkg/MSC_CLD_210525121000_100000000_08_raw.fits"],
        flat_path=[
            "/home/csstpipeline/data/bkg/MSC_CLF_210525120500_100000000_08_raw.fits"],
        output_path="/home/csstpipeline/data/L05_test/",
    )

