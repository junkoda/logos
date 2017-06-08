default: build install

.PHONY: build install

build:
	python3 setup.py config build_ext --inplace

install:
	python3 setup.py install

