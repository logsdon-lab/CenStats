.PHONY: build test venv clean install

BIN=$(shell pwd)/venv/bin/
PROJECT_NAME=centromere_status_checker

# test:

build:
	$(MAKE) clean
	$(BIN)python3 -m pip install --upgrade build
	$(BIN)python3 -m build

install:
	$(BIN)python3 -m pip install $(shell find dist -name "*.whl" | sort -r | head -1)

venv:
	python3 -m venv venv

clean:
	rm -rf dist/

upload:
	$(BIN)python3 -m pip install --upgrade twine
	$(BIN)python3 -m twine upload --repository pypi dist/*
