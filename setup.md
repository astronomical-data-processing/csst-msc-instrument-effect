## 探测器效应修正

### 启动容器
```sh
docker run -it --entrypoint bash --network=host --rm  -e GRPC_SERVER=192.168.122.185:50051  -v /data/gateway/log:/var/log/openccf -v /data/C3_pipeline/backup/CAL:/home/csstpipeline/data/bkg:ro -v /data/C3_pipeline/backup/CAL:/CAL:ro -v /data/C3_pipeline/backup:/hidden -v /data/single-image-reduction/tasks/detector-effect-correction-input:/tasks -v /data/single-image-reduction/entity/L00:/L0 -v /data/single-image-reduction/entity/L05:/home/csstpipeline/data/L05_test -v /data/single-image-reduction/tasks/WCS-calibration-input:/next-stage -w /input csst/detector-effect-correction
```

- 0级数据映射为 /L0
- 0.5级数据映射为 /home/csstpipeline/data/L05_test 

### 容器内运行命令
```sh

run.sh 10000000106

```

### 注意事项及存在问题
- 代码与数据分目录存放；
- 不同类型的数据也分目录存放；
- 临时文件与最终文件也分目录存放

### github推送
