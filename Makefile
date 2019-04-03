DOCKER_IMAGE ?= quay.io/pypa/manylinux1_x86_64

darwin:
	python setup.py bdist_wheel

linux:
	docker run --rm -v `pwd`:/io ${DOCKER_IMAGE} /io/build-wheels.sh 
