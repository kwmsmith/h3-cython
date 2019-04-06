DOCKER_IMAGE ?= quay.io/pypa/manylinux1_x86_64

darwin: clean-cache
	python setup.py bdist_wheel

linux: clean-cache
	docker run --rm -v `pwd`:/io ${DOCKER_IMAGE} /io/build-wheels.sh 

local-install:
	pip install h3cy --ignore-installed --no-index --find-links ./dist
.PHONY: local-install

test: local-install
	pytest test/test_h3.py
	# pytest --verbose --tb=line -s test/

clean-cache:
	find . -iname '__pycache__' -exec rm -r '{}' +
.PHONY: clean-cache

clean-all:
	-rm -r _skbuild CMakeFiles CMakeCache.txt
