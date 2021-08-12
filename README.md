# csst_msc_instrument_effect 仪器效应改正

用于改正CSST的0级数据中的仪器效应问题，如本底、暗流、平场，并标记出坏像元、热像元、暖像元、溢出像元、宇宙线。也会运算出权重weight。

## 在命令行运行

```bash
sh csst_msc_iec.sh [id]
```

## 使用python运行
```python
from csst_msc_instrument_effect.msc_iec import InstrumentEffectCorrection


iec = InstrumentEffectCorrection(
    data_path=data_filename,
    bias_path=bias_filename,
    dark_path=dark_filename,
    flat_path=flat_filename,
    cray_path=cray_filename,
    output_path=output_path,
    config_path=config_path,
)
iec.run()
```
