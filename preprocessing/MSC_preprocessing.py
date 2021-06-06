from preprocessing.MSC_crmask import CRMask
import numpy as np
from astropy.io import fits
import os


class csst_msc_instrument_effect:
    def __init__(self, **kwarg) -> None:
        # data
        if "data" in kwarg:
            with fits.open(kwarg["data"]) as hdul:
                self.phu = hdul[0].header
                self.du = hdul[1].data.astype(int)
                self.hu = hdul[1].header
        else:
            raise "no data"
        # RDNOISE
        self.RDNOISE = kwarg["RDNOISE"] if "RDNOISE" in kwarg else "RDNOISE"
        # GAIN
        self.GAIN = kwarg["GAIN"] if "GAIN" in kwarg else "GAIN"
        # EXPTIME
        self.EXPTIME = kwarg["EXPTIME"] if "EXPTIME" in kwarg else "EXPTIME"
        # cray
        self.cray = self.fix_cray(kwarg["cray"]) if "cray" in kwarg else 0
        # bias
        self.bias = self.fix_bias(kwarg["bias"]) if "bias" in kwarg else 0
        # dark
        self.dark = self.fix_dark(kwarg["dark"]) if "dark" in kwarg else 0
        # flat
        self.flat = self.fix_flat(kwarg["flat"]) if "flat" in kwarg else 1


    def fix_bias(
        self, paths,
    ):
        bias = []
        for path in paths:
            du = fits.getdata(path)
            du = du.astype(int)
            bias.append(du)
        return np.median(bias, axis=0)

    def fix_dark(
        self, paths,
    ):
        dark = []
        for path in paths:
            with fits.open(path) as hdul:
                du = hdul[1].data
                hu = hdul[0].header
                du = du.astype(int)
                du = du - self.bias
                du = du - self.cray / self.hu[self.GAIN]
                du = du / hu[self.EXPTIME] * self.phu[self.EXPTIME]
                dark.append(du)
        return np.median(dark, axis=0)

    def fix_flat(
        self, paths,
    ):
        flat = []
        for path in paths:
            with fits.open(path) as hdul:
                du = hdul[1].data
                hu = hdul[0].header
                du = du.astype(int)
                du = du - self.bias
                du = du / hu[self.EXPTIME] * self.phu[self.EXPTIME]
                du = du - self.dark
                du = du / np.median(du)
                flat.append(du)
        return np.median(flat, axis=0)

    def fix_cray(
        self, path,
    ):
        du = fits.getdata(path)
        du = du.astype(int)
        return du

    def fix_data(self,):
        return (self.du - self.bias - self.dark) / self.flat

    def fix_mask(self,):
        mask = np.zeros_like(self.du, dtype=np.uint16)
        # 00000001:   坏像元    因为探测器本身原因造成不可用于有效科学研究的像元, 像元响应小于0.5中值 大于1.5中值
        med = np.median(self.flat)
        msk = (self.flat < 0.5 * med) | (1.5 * med < self.flat)
        mask = mask | (msk * 1)
        # 00000010:   热像元    因为探测器本身原因造成的影响科学研究结果的像元. 像元在150秒积分时间内, 暗流计数大于探测器平均读出噪声的平方.
        dark = self.dark.copy()
        dark[dark < 0] = 0
        msk = 1 * self.hu[self.RDNOISE] ** 2 <= dark  # 不确定是否包含 暂定包含
        mask = mask | (msk * 2)
        # 00000100:   暖像元    因为探测器本身原因造成的影响科学研究结果的像元. 像元的150秒积分时间内, 暗流计数大于0.5被读出噪声的平方, 但小于1被读出噪声的平方.
        msk = (0.5 * self.hu[self.RDNOISE] ** 2 < dark) & (
            dark < 1 * self.hu[self.RDNOISE] ** 2
        )
        mask = mask | (msk * 4)
        # 00001000:   饱和溢出像元  饱和像元及流量溢出污染的像元.
        msk = self.du == 65535
        mask = mask | (msk * 8)
        # 00010000:   宇宙线像元    宇宙线污染的像元
        crobj = CRMask(self.data)
        msk, data = crobj.cr_mask()
        msk = msk[1].data
        data = data[1].data
        mask = mask | (msk * 16)
        # 00100000:   卫星或者人造移动天体轨迹污染的像元.
        pass
        # 01000000:   鬼像污染的像元.
        pass
        # 10000000:   散射光污染的像元(包括dragon breath)
        pass
        return mask, data

    def fix_weight(self, data):
        weight = 1 / (self.hu[self.GAIN] * data + self.hu[self.RDNOISE] ** 2)
        weight[self.mask > 0] = 0
        return weight

    def run(self):
        # fix data
        self.data = self.fix_data()
        # mask
        self.mask, self.data = self.fix_mask()
        # weight
        data = self.data.copy()
        data[data < 0] = 0
        self.weight = self.fix_weight(data)

if __name__ == "__main__":
    args = {
        "data": "/mnt/storage-data/xiezhou/TEST_20210602/MSC_MS_210525121500_100000001_08_raw.fits",
        "bias": [
            "/mnt/storage-data/xiezhou/TEST_20210602/MSC_CLB_210525120000_100000000_08_raw.fits",
        ],
        "dark": [
            "/mnt/storage-data/xiezhou/TEST_20210602/MSC_CLD_210525121000_100000000_08_raw.fits",
        ],
        "flat": [
            "/mnt/storage-data/xiezhou/TEST_20210602/MSC_CLF_210525120500_100000000_08_raw.fits",
        ],
        "cray": "/mnt/storage-data/xiezhou/TEST_20210602/MSC_CRD_210525121000_100000000_08_raw.fits",
        "RDNOISE": "RDNOISE1",
        "GAIN": "GAIN1",
    }
    pl = csst_msc_instrument_effect(**args)
    pl.run()
    if not os.path.exists(".local"):
        os.mkdir(".local")
    data = fits.HDUList(
        [
            fits.PrimaryHDU(header=pl.phu),
            fits.ImageHDU(pl.data.astype(np.float32), pl.hu),
        ]
    )
    filename = os.path.basename(args["data"]).replace("raw", "fix")
    filename = os.path.join('.local', filename)
    data.writeto(filename)
    bias = fits.HDUList(
        [
            fits.PrimaryHDU(header=pl.phu),
            fits.ImageHDU(data=pl.bias.astype(np.int16), header=pl.hu),
        ]
    )
    filename = os.path.basename(args["bias"]).replace("raw", "fix")
    filename = os.path.join('.local', filename)
    bias.writeto(filename)
    dark = fits.HDUList([fits.PrimaryHDU(pl.dark.astype(np.int16))])
    filename = os.path.basename(args["dark"]).replace("raw", "fix")
    filename = os.path.join('.local', filename)
    dark.writeto(filename)
    flat = fits.HDUList([fits.PrimaryHDU(pl.flat.astype(np.float32))])
    filename = os.path.basename(args["flat"]).replace("raw", "fix")
    filename = os.path.join('.local', filename)
    flat.writeto(filename)
    flg = fits.HDUList([fits.PrimaryHDU(data=pl.mask.astype(np.uint16), header=pl.hu)])
    filename = os.path.basename(args["data"]).replace("raw", "flg")
    filename = os.path.join('.local', filename)
    flg.writeto(filename)
    weight = fits.HDUList(
        [fits.PrimaryHDU(data=pl.weight.astype(np.float32), header=pl.hu)]
    )
    filename = os.path.basename(args["data"]).replace("raw", "wht")
    filename = os.path.join('.local', filename)
    weight.writeto(filename)
