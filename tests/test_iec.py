import json
import unittest

from csst_msc_instrument_effect.msc_iec import InstrumentEffectCorrection


class MSCPipelineTestCase(unittest.TestCase):
    def setUp(self):
        self.iec = InstrumentEffectCorrection(json_file="request.json", output_path=".local")

    def test_all(self):
        self.iec.set_data()
        self.iec.set_cray()
        self.iec.set_bias()
        self.iec.set_dark()
        self.iec.set_flat()
        self.iec.fix_data()
        self.iec.set_flag()
        self.iec.set_weight()
        from glob import glob
        test_file_list = glob(".local/*.fits")
        if len(test_file_list) > 0:
            import os
            for path in test_file_list:
                os.remove(path)
        self.iec.save()

