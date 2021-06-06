# csst_msc_instrument_effect 仪器效应改正

用于改正CSST的0级数据中的仪器效应问题，如本底、暗流、平场，并标记出坏像元、热像元、暖像元、溢出像元、宇宙线。也会运算出权重weight。

## 使用方法：

首先需要准备传入参数
```json
{
    "file_data_fullname": "data.fits",  # 科学图像路径
    "file_bias_fullname_list": ["bias0.fits","bias1.fits","bias2.fits",],  # 本底图像路径列表
    "file_dark_fullname_list": ["dark0.fits","dark1.fits","dark2.fits",],  # 暗场图像路径列表
    "file_flat_fullname_list": ["flat0.fits","flat1.fits","flat2.fits",],  # 平场图像路径列表
    "file_cray_fullname": "cray.fits",  # 宇宙线图像路径列表（暂时）
    "fits_field_read_noise": "RDNOISE",  # 头文件中读出噪声字段名
    "fits_field_gain": "GAIN",  # 头文件中增益字段名
    "fits_field_exposure_time": "EXPTIME"  # 头文件中曝光时间字段名
}
```
在命令行运行
```bash
csst-msc-iec -j json_file --json=json_file -o output_path --output=output_path
```