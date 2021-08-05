import sys
from datetime import datetime

from csst_dfs_api.facility.level0 import Level0DataApi
from csst_dfs_api.facility.level0prc import Level0PrcApi
from csst_dfs_api.facility.calmerge import CalMergeApi
from csst_msc_instrument_effect.msc_iec import InstrumentEffectCorrection

l0Api = Level0DataApi()
prcApi = Level0PrcApi()
calApi = CalMergeApi()

# -----------------------------------------------------
# get data
l0_id = sys.argv[1]

rec = l0Api.get(level0_id=l0_id)
print(rec)
print("Find Level0 Data :" + rec.data.file_path)

bias_rec = calApi.get_latest_by_l0(level0_id=l0_id, ref_type="bias")
dark_rec = calApi.get_latest_by_l0(level0_id=l0_id, ref_type="dark")
flat_rec = calApi.get_latest_by_l0(level0_id=l0_id, ref_type="flat")
cray_rec = calApi.get_latest_by_l0(level0_id=l0_id, ref_type="cray")  # temp

print("Find Level0 Bias :" + bias_rec.data.file_path)
print("Find Level0 Dark :" + dark_rec.data.file_path)
print("Find Level0 Flat :" + flat_rec.data.file_path)
# ------------------------------------------------------
# to-do data processing here
iec = InstrumentEffectCorrection(
    data_path=rec.data.file_path,
    bias_path=bias_rec.data.file_path,
    dark_path=dark_rec.data.file_path,
    flat_path=flat_rec.data.file_path,
    output_path="/home/csstpipeline/data/L05_test/",
    config_path="/home/csstpipeline/csst-msc-instrument-effect-master/MSC_crmask.ini",
    cray_path=cray_rec.data.file_path,  # temp
)
iec.run()

params_file_path = iec.config_path
result_file_path = ', '.join([
    iec.data_output,
    iec.flag_output,
    iec.weight_output,
    iec.header_output,
])
dt = datetime.now()
# ------------------------------------------------------
# write result
rec = prcApi.write(
    level0_id=l0_id,
    pipeline_id="P1",
    prc_module="MSC-IE",  # csst-msc-instrument-effect
    params_file_path=params_file_path,
    prc_status=0,
    prc_time=dt.strftime('%Y-%m-%d %H:%M:%S'),
    result_file_path=result_file_path)
print('Write Result:', rec)
