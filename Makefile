IMAGE=hub.cstcloud.cn/csst/detector-effect-correction:dev

build:
	docker build --network=host -t $(IMAGE) .

push:
	docker push $(IMAGE)

clean:
	docker rmi $(IMAGE)

rmi:
	for H in h0 h1 h2 ; do echo $$H ; ssh $$H docker rmi $(IMAGE) ; done
