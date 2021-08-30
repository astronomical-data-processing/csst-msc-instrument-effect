import os
from typing import List
from astropy.io.fits import header
from datetime import datetime

import numpy as np
from astropy.io import fits

from csst_msc_instrument_effect.msc_crmask import CRMask


class InstrumentEffectCorrection:
    def __init__(
        self,
        data_path,
        bias_path,
        dark_path,
        flat_path,
        output_path,
        config_path,
        cray_path,
    ) -> None:
        self.data_path = data_path
        self.bias_path = bias_path
        self.dark_path = dark_path
        self.flat_path = flat_path
        self.cray_path = cray_path
        self.output = output_path
        self.config_path = config_path
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

    def set_bias(self, mode):
        # bias
        def func(path):
            du = fits.getdata(path)
            return du.astype(int)

        if isinstance(self.bias_path, list) and len(self.bias_path) != 1:
            bias = np.array([*map(func, self.bias_path)])
            self.bias, var = array_combine(bias, mode)
        elif isinstance(self.bias_path, list) and len(self.bias_path) == 1:
            self.bias = func(self.bias_path[0])
        else:
            self.bias = func(self.bias_path)

    def set_dark(self, mode):
        # dark
        def func(path):
            with fits.open(path) as hdul:
                du = hdul[1].data
                hu = hdul[0].header
            du = du.astype(int)
            du = du - self.bias
            du = du - self.cray / self.image_header[self.GAIN]  # tmp
            return du / hu[self.EXPTIME] * self.primary_header[self.EXPTIME]

        if isinstance(self.dark_path, list) and len(self.dark_path) != 1:
            dark = np.array([*map(func, self.dark_path)])
            self.dark, var = array_combine(dark, mode)
        elif isinstance(self.dark_path, list) and len(self.dark_path) == 1:
            self.dark = func(self.dark_path[0])
        else:
            self.dark = func(self.dark_path)

    def set_flat(self, mode):
        # flat
        def func(path):
            with fits.open(path) as hdul:
                du = hdul[1].data
                hu = hdul[0].header
            du = du.astype(int)
            du = du - self.bias
            du = du / hu[self.EXPTIME] * self.primary_header[self.EXPTIME]
            du = du - self.dark
            return du / np.median(du)

        if isinstance(self.flat_path, list) and len(self.flat_path) != 1:
            flat = np.array([*map(func, self.flat_path)])
            self.flat, var = array_combine(flat, mode)
        elif isinstance(self.flat_path, list) and len(self.flat_path) == 1:
            self.flat = func(self.flat_path[0])
        else:
            self.flat = func(self.flat_path)

    def set_cray(self):
        # cray
        self.cray = fits.getdata(self.cray_path).astype(int)

    def fix_data(self,):
        self.data_fix0 = np.divide(
            self.data_raw - self.bias - self.dark,
            self.flat,
            out=np.zeros_like(self.data_raw, float),
            where=(self.flat != 0),
        )

    def set_flag(self,):
        flag = np.zeros_like(self.data_raw, dtype=np.uint16)
        # 00000001:   坏像元
        # 因为探测器本身原因造成不可用于有效科学研究的像元, 像元响应小于0.5中值 大于1.5中值
        med = np.median(self.flat)
        flg = (self.flat < 0.5 * med) | (1.5 * med < self.flat)
        flag = flag | (flg * 1)
        # 00000010:   热像元
        # 因为探测器本身原因造成的影响科学研究结果的像元. 像元在150秒积分时间内, 暗流计数大于探测器平均读出噪声的平方.
        dark = self.dark.copy()
        dark[dark < 0] = 0
        flg = 1 * self.image_header[self.RDNOISE] ** 2 <= dark  # 不确定是否包含 暂定包含
        flag = flag | (flg * 2)
        # 00000100:   暖像元
        # 因为探测器本身原因造成的影响科学研究结果的像元. 像元的150秒积分时间内, 暗流计数大于0.5被读出噪声的平方, 但小于1被读出噪声的平方.
        flg = (0.5 * self.image_header[self.RDNOISE] ** 2 < dark) & (
            dark < 1 * self.image_header[self.RDNOISE] ** 2
        )
        flag = flag | (flg * 4)
        # 00001000:   饱和溢出像元
        # 饱和像元及流量溢出污染的像元.
        flg = self.data_raw == 65535
        flag = flag | (flg * 8)
        # 00010000:   宇宙线像元
        # 宇宙线污染的像元
        crobj = CRMask(self.data_fix0, config_path=self.config_path)
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
            self.image_header[self.GAIN] * data + self.image_header[self.RDNOISE] ** 2
        )
        weight[self.flag > 0] = 0
        self.weight = weight

    def save(self):
        data_filename = self.data_path
        data_basename = os.path.basename(data_filename).replace("raw", "img")
        data_output = os.path.join(self.output, data_basename)
        data_header = self.image_header.copy()
        data_header["EXTNAME"] = "img"
        data_header[self.GAIN] = (
            data_header[self.GAIN] * self.primary_header[self.EXPTIME]
        )
        data_header["BUNIT"] = "e/s"
        data_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.data_fix0.astype(np.float32)
                    / self.primary_header[self.EXPTIME],
                    header=data_header,
                ),
            ]
        )
        data_fits.writeto(data_output, overwrite=True)
        self.data_output = data_output

        flag_output = data_output.replace("img", "flg")
        flag_header = self.image_header.copy()
        flag_header["EXTNAME"] = "flg"
        flag_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(data=self.flag.astype(np.uint16), header=flag_header),
            ]
        )
        flag_fits.writeto(flag_output, overwrite=True)
        self.flag_output = flag_output

        weight_output = data_output.replace("img", "wht")
        weight_header = self.image_header.copy()
        weight_header["EXTNAME"] = "wht"
        weight_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.weight.astype(np.float32), header=weight_header
                ),
            ]
        )
        weight_fits.writeto(weight_output, overwrite=True)
        self.weight_output = weight_output

        header_name = data_basename.replace(".fits", ".head")
        header_output = os.path.join(self.output, header_name)
        hu = self.image_header
        hu["INSTRU_S"] = (0, "instrument effect status")
        hu["INST_TOL"] = (1234, "instrument effect operation time")
        hu["INSTRU_V"] = ("0.1", "instrument effect version")
        hu["INSTRU_P"] = ("test.conf", "instrument effect config file")
        header_fits = fits.HDUList([fits.PrimaryHDU(header=hu)])
        header_fits.writeto(header_output, overwrite=True)
        self.header_output = header_output

    def run(self):
        self.set_data()
        self.set_cray()
        self.set_bias(mode="median_clip")
        self.set_dark(mode="mean_clip")
        self.set_flat(mode="mean_clip")
        self.fix_data()
        self.set_flag()
        self.set_weight()
        self.save()


def array_combine(ndarray, mode="mean"):
    """ Function to combine 3-D data array

    Parameters
    ----------
    ndarray: array, input data cube (3D)
    model: mean, median, sum, mean_clip, median_clip, default is mean
    """
    var = np.std(ndarray, axis=0)
    ndarray_sort = np.sort(ndarray, axis=0)  ### sort the flux for each pixel
    ndarray_cut = ndarray_sort[1:-1]  ### cut the min and max flux for each pixel

    if mode == "median":
        array = np.median(ndarray_sort, axis=0)
    elif mode == "median_clip":
        array = np.median(ndarray_cut, axis=0)
    elif mode == "sum":
        array = np.sum(ndarray_sort, axis=0)
    elif mode == "mean":
        array = np.mean(ndarray_sort, axis=0)
    elif mode == "mean_clip":
        array = np.mean(ndarray_cut, axis=0)
    return array, var

