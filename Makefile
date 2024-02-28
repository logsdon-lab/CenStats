.PHONY: build test venv clean

BIN=venv/bin/
PROJECT_NAME=centromere_status_checker

# test:

build:
	$(BIN)python3 -m pip install --upgrade build
	$(BIN)python3 -m build

venv:
	python3 -m venv venv

clean:
	rm -rf dist/ $(PROJECT_NAME)*/ venv/

upload:
	$(BIN)python3 -m pip install --upgrade twine
	$(BIN)python3 -m twine upload --repository pypi dist/*
