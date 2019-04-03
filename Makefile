DOCKER_IMAGE ?= quay.io/pypa/manylinux1_x86_64

darwin:
	python setup.py bdist_wheel

linux:
	docker run --rm -v `pwd`:/io ${DOCKER_IMAGE} /io/build-wheels.sh 

local-install:
	pip install h3cy --ignore-installed --no-index --find-links ./dist
.PHONY: local-install

test: local-install
	pytest --tb=line test/

clean-all:
	-rm -r _skbuild CMakeFiles CMakeCache.txt
