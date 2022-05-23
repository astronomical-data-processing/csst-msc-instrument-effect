# 运行: sh /home/csstpipeline/csst-msc-instrument-effect-main/run.sh [输入文件路径] [输出文件路径] [背景文件路径] [配置文件路径]
# 结果默认保存在/data/cali_20211012/output/
# 背景文件路径默认在/data/test20211012/ref/
# 配置文件路径默认在/home/csstpipeline/csst-msc-instrument-effect-main/MSC_crmask.ini
input=${1}
output=${2:-'/data/cali_20211012/output/'}
ref=${3:-'/data/test20211012/ref/'}
conf=${4:-"/home/csstpipeline/csst-msc-instrument-effect-main/MSC_crmask.ini"}
echo ${input}
echo ${output}
echo ${ref}
echo ${conf}
python /home/csstpipeline/csst-msc-instrument-effect-main/run.py -i ${input} -o ${output} -r ${ref} -c ${conf}
