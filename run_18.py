import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
from glob import glob
from astropy.io.fits import getdata, getheader
from concurrent.futures import ThreadPoolExecutor
from time import time

from csst_msc_instrument_effect.msc_iec import InstrumentEffectCorrection


def run(data_path, output, ref, conf):
    header = getheader(data_path)
    number = header['DETECTOR'][-2:]
    number_list = ['06', '07', '08', '09', '11', '12', '13', '14', '15',
                   '16', '17', '18', '19', '20', '22', '23', '24', '25']
    if number not in number_list:
        print(number, 'not in')
        return
    bias = getdata(glob(os.path.join(ref, "*CLB*_" + number + '_*'))[0])
    dark = getdata(glob(os.path.join(ref, "*CLD*_" + number + '_*'))[0])
    flat = getdata(glob(os.path.join(ref, "*CLF*_" + number + '_*'))[0])
    iec = InstrumentEffectCorrection(
        data_input=data_path,
        bias_input=bias,
        dark_input=dark,
        flat_input=flat,
        output_path=output,
        config_path=conf
    )
    iec.run()
    print(number, 'done')


if __name__ == '__main__':
    input_path = 'D:/Desktop/data/'
    output_path = 'D:/Desktop/data/output/'
    ref_path = r'D:/Desktop/data/ref/'
    conf_path = 'D:/Desktop/test/csst-msc-instrument-effect/CSST_2021-12-30_CCD23_epoch20.pth'
    with ThreadPoolExecutor(18) as tpe:
        for path in glob(os.path.join(input_path, '*_MS_SCI_*.fits')):
            tpe.submit(run, path, output_path, ref_path, conf_path)
