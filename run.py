import os
from glob import glob
import argparse
from astropy.io.fits import getdata, getheader
from time import time

from csst_msc_instrument_effect.msc_iec import InstrumentEffectCorrection


def run(data_path, output, ref, conf):
    header = getheader(data_path)
    number = header['DETECTOR'][-2:]
    # number = os.path.basename(data_path).split('_')[4]  # 原来版本
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', help='fits file path or folder path', required=True)
    parser.add_argument(
        '-o', '--output', help='output path', default='/data/L05/')
    parser.add_argument(
        '-r', '--ref', help='ref path', default='/data/ref/')
    parser.add_argument(
        '-c', '--conf', help='config path', default="/home/csstpipeline/csst-msc-instrument-effect-main/MSC_crmask.ini")
    args = parser.parse_args()
    if os.path.isdir(args.input):
        for path in glob(os.path.join(args.input.strip(), '*_MS_SCI_*.fits')):
        #for path in glob(os.path.join(args.input.strip(), '*_MS_*.fits')):
            run(path, args.output, args.ref, args.conf)
    else:
        print('run test')
        begin = time()
        run(args.input.strip(), args.output.strip(), args.ref.strip(), args.conf.strip())
        print(time() - begin)


if __name__ == '__main__':
    # python run.py -i "/data/test20211012/30s/MSC_MS_210525220000_100000020_13_raw.fits" -o "/data/test20211012/output/150s/"
    # python run.py -i "/data/test20211012/150s/" -o '/data/test20211012/output/150s/'
    # python run.py -i "/data/cali_20211012/L0/150s/" -o '/data/cali_20211012/output/'
    # python run.py -i '/data/test20211012/30s/' -o '/data/test20211012/output/30s/'
    # python run.py -i "/data/sim_data/20220413/MSC_0000100/CSST_MSC_MS_CRS_20270810081950_20270810082220_100000100_13_L0_1.fits" -o "/home/csstpipeline/test/"
    main()
