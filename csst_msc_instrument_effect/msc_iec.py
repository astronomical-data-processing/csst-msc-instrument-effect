import os
import json
import argparse

import numpy as np
from astropy.io import fits

from csst_msc_instrument_effect.msc_crmask import CRMask


class InstrumentEffectCorrection:
    def __init__(self, json_file, output_path) -> None:
        with open(json_file) as f:
            self.json = json.load(f)
        if "file_data_fullname" not in self.json:
            raise "no file_data_fullname"
        self.output = output_path
        # RDNOISE
        self.RDNOISE = (
            self.json["fits_field_read_noise"]
            if "fits_field_read_noise" in self.json
            else "RDNOISE"
        )
        # GAIN
        self.GAIN = (
            self.json["fits_field_gain"] if "fits_field_gain" in self.json else "GAIN"
        )
        # EXPTIME
        self.EXPTIME = (
            self.json["fits_field_exposure_time"]
            if "fits_field_exposure_time" in self.json
            else "EXPTIME"
        )

    def set_data(self):
        # data
        with fits.open(self.json["file_data_fullname"]) as hdul:
            self.primary_header = hdul[0].header
            self.data_raw = hdul[1].data.astype(int)
            self.image_header = hdul[1].header

    def set_cray(self):
        # cray
        if "file_cray_fullname" not in self.json:
            self.cray = np.zeros_like(self.data_raw)
        else:
            self.cray = fits.getdata(self.json["file_cray_fullname"]).astype(int)

    def set_bias(self):
        # bias
        if "file_bias_fullname_list" not in self.json:
            self.bias = np.zeros_like(self.data_raw)
        else:
            bias = []
            for path in self.json["file_bias_fullname_list"]:
                du = fits.getdata(path)
                du = du.astype(int)
                bias.append(du)
            self.bias = np.median(bias, axis=0)

    def set_dark(self):
        # dark
        if "file_dark_fullname_list" not in self.json:
            self.dark = np.zeros_like(self.data_raw)
        else:
            dark = []
            for path in self.json["file_dark_fullname_list"]:
                with fits.open(path) as hdul:
                    du = hdul[1].data
                    hu = hdul[0].header
                    du = du.astype(int)
                    du = du - self.bias
                    du = du - self.cray / self.image_header[self.GAIN]
                    du = du / hu[self.EXPTIME] * self.primary_header[self.EXPTIME]
                    dark.append(du)
            self.dark = np.median(dark, axis=0)

    def set_flat(self):
        # flat
        if "file_flat_fullname_list" not in self.json:
            self.flat = np.ones_like(self.data_raw)
        else:
            flat = []
            for path in self.json["file_flat_fullname_list"]:
                with fits.open(path) as hdul:
                    du = hdul[1].data
                    hu = hdul[0].header
                    du = du.astype(int)
                    du = du - self.bias
                    du = du / hu[self.EXPTIME] * self.primary_header[self.EXPTIME]
                    du = du - self.dark
                    du = du / np.median(du)
                    flat.append(du)
            self.flat = np.median(flat, axis=0)

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
        self.data_fix1 = data_fits[1].data

    def set_weight(self):
        weight = 1 / (
            self.image_header[self.GAIN] * self.data_fix0
            + self.image_header[self.RDNOISE] ** 2
        )
        weight[self.flag > 0] = 0
        self.weight = weight

    def save(self):
        data_filename = self.json["file_data_fullname"]
        data_basename = os.path.basename(data_filename).replace("raw", "fix")
        data_output = os.path.join(self.output, data_basename)
        data_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.data_fix1.astype(np.float32), header=self.image_header
                ),
            ]
        )
        data_fits.writeto(data_output)

        flag_output = data_output.replace("fix", "flg")
        flag_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.flag.astype(np.uint16), header=self.image_header
                ),
            ]
        )
        flag_fits.writeto(flag_output)

        weight_output = data_output.replace("fix", "wht")
        weight_fits = fits.HDUList(
            [
                fits.PrimaryHDU(header=self.primary_header),
                fits.ImageHDU(
                    data=self.weight.astype(np.float32), header=self.image_header
                ),
            ]
        )
        weight_fits.writeto(weight_output)

        if "file_bias_fullname_list" in self.json:
            bias_filename = self.json["file_bias_fullname_list"][0]
            bias_basename = os.path.basename(bias_filename).replace("raw", "fix")
            bias_output = os.path.join(self.output, bias_basename)
            bias_fits = fits.HDUList([fits.PrimaryHDU(data=self.bias)])
            bias_fits.writeto(bias_output)

        if "file_dark_fullname_list" in self.json:
            dark_filename = self.json["file_dark_fullname_list"][0]
            dark_basename = os.path.basename(dark_filename).replace("raw", "fix")
            dark_output = os.path.join(self.output, dark_basename)
            dark_fits = fits.HDUList([fits.PrimaryHDU(data=self.dark)])
            dark_fits.writeto(dark_output)

        if "file_flat_fullname_list" in self.json:
            flat_filename = self.json["file_flat_fullname_list"][0]
            flat_basename = os.path.basename(flat_filename).replace("raw", "fix")
            flat_output = os.path.join(self.output, flat_basename)
            flat_fits = fits.HDUList([fits.PrimaryHDU(data=self.flat)])
            flat_fits.writeto(flat_output)


def main():
    args = argparse.ArgumentParser()
    args.add_argument("-j", "--json", required=True)
    args.add_argument("-o", "--output", required=True)
    args = args.parse_args()
    obj = InstrumentEffectCorrection(json_file=args.json, output_path=args.output)
    obj.set_data()
    obj.set_cray()
    obj.set_bias()
    obj.set_dark()
    obj.set_flat()
    obj.fix_data()
    obj.set_flag()
    obj.set_weight()
    obj.save()


if __name__ == "__main__":
    main()
