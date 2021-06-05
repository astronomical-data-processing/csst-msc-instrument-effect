import unittest

from lib.MSC_preprocessing import Pipeline
import os


class MSCPipelineTestCase(unittest.TestCase):
    def setUp(self):
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
        self.api = Pipeline(**args)

    def test_run(self):
        self.api.run()
        print('run end')