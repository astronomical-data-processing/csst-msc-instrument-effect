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
    data_path=data_filename,    # 科学图fits文件路径
    bias_path=bias_filename,    # 本底fits文件路径(可以是路径列表)
    dark_path=dark_filename,    # 暗场fits文件路径(可以是路径列表)
    flat_path=flat_filename,    # 平场fits文件路径(可以是路径列表)
    cray_path=cray_filename,    # 模拟宇宙线文件路径(因为模拟数据的问题, 暂时添加, 之后会去除)
    output_path=output_path,    # 输出路径
    config_path=config_path,    # crmask配置文件路径, 指的就是MSC_crmask.ini
)
iec.run()
```
生成的文件会有img, flag, weight, header四个

## 使用BKGscript.py

运行BKGscript.py 更改main函数内input_path与save_path与number_list三个变量
input_path为参考源文件路径
save_path为输出路径
number_list为需要进行合并的参考文件编号
