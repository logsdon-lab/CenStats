.PHONY: build test venv dev clean install upload

BIN=$(shell pwd)/venv/bin/
PROJECT_NAME=censtats

test:
	$(BIN)python3 -m pip install pytest
	$(BIN)python3 -m pytest -vv

build:
	$(MAKE) clean
	$(BIN)python3 -m pip install --upgrade build
	$(BIN)python3 -m build

install:
	$(BIN)python3 -m pip uninstall -y $(PROJECT_NAME)
	$(BIN)python3 -m pip install $(shell find dist -name "*.whl" | sort -r | head -1)

venv:
	python3 -m venv venv

dev:
	$(BIN)python3 -m pip install -r requirements.txt
	$(BIN)python3 -m pip install pytest

clean:
	rm -rf dist/ *.egg-info .*_cache

upload:
	$(BIN)python3 -m pip install --upgrade twine
	$(BIN)python3 -m twine upload --repository pypi dist/*
