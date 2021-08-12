from csst_msc_instrument_effect.msc_iec import InstrumentEffectCorrection



iec = InstrumentEffectCorrection(
    data_path="/data/L0/MSC_MS_210525121500_100000001_08_raw.fits",
    bias_path="/data/BKG/MSC_CLB_210525120000_100000000_08_raw.fits",
    dark_path="/data/BKG/MSC_CLD_210525121000_100000000_08_raw.fits",
    flat_path="/data/BKG/MSC_CLF_210525120500_100000000_08_raw.fits",
    cray_path="/data/BKG/MSC_CRD_210525121000_100000000_08_raw.fits",
    output_path="/data/test/",
    config_path="MSC_crmask.ini",
)
iec.run()