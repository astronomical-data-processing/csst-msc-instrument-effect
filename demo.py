from csst_msc_instrument_effect.msc_iec import InstrumentEffectCorrection



iec = InstrumentEffectCorrection(
    data_input="../data/test/MSC_MS_210525160000_100000022_24_raw.fits",#"../data/test/MSC_MS_210525123000_100000001_07_raw.fits",#"/data/L0/MSC_MS_210525121500_100000001_08_raw.fits",
    bias_input="/data/BKG/MSC_CLB_210525120000_100000000_24_raw.fits",
    dark_input="/data/BKG/MSC_CLD_210525121000_100000000_24_raw.fits",
    flat_input="/data/BKG/MSC_CLF_210525120500_100000000_24_raw.fits",
    cray_input="/data/BKG/MSC_CRD_210525121000_100000000_24_raw.fits",
    output_path="/data/test/",
    config_path="MSC_crmask.ini",
)
iec.run()
