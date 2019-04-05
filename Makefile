DOCKER_IMAGE ?= quay.io/pypa/manylinux1_x86_64

darwin:
	-rm -r test/__pycache__
	python setup.py bdist_wheel

linux:
	-rm -r test/__pycache__
	docker run --rm -v `pwd`:/io ${DOCKER_IMAGE} /io/build-wheels.sh 

local-install:
	pip install h3cy --ignore-installed --no-index --find-links ./dist
.PHONY: local-install

test: local-install
	pytest --verbose --tb=line -s test/

clean-cache:
	find . -iname '__pycache__'

clean-all:
	-rm -r _skbuild CMakeFiles CMakeCache.txt
