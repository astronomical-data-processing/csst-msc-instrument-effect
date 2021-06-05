# MSC_INS_EFFECT 仪器效应改正

用于改正CSST的0级数据中的仪器效应问题，如本底、暗流、平场，并标记出坏像元、热像元、暖像元、溢出像元、宇宙线。也会运算出权重weight。

## 使用方法：

首先需要准备传入参数
```python
args = {
    "data": "data.fits",  # 科学图像路径
    "bias": ["bias0.fits","bias1.fits","bias2.fits",],  # 本底图像路径列表
    "dark": ["dark0.fits","dark1.fits","dark2.fits",],  # 暗场图像路径列表
    "flat": ["flat0.fits","flat1.fits","flat2.fits",],  # 平场图像路径列表
    "cray": "cray.fits",  # 宇宙线图像路径列表（暂时）
    "RDNOISE": "RDNOISE",  # 头文件中读出噪声字段名
    "GAIN": "GAIN",  # 头文件中增益字段名
    "EXPTIME": "EXPTIME",  # 头文件中曝光时间字段名
}
from preprocessing.MSC_preprocessing import Pipeline  # 载入包

pl = Pipeline(**args)  # 创建管线对象
pl.run()  # 运行管线
pl.data  # 科学图像改正后的结果
pl.mask  # 像元标记结果
pl.weight  # 权重结果
```