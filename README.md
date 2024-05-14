# `CenStats`
[![CI](https://github.com/logsdon-lab/centromere-status-checker/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/centromere-status-checker/actions/workflows/main.yml)

Centromere statistics toolkit.

* `status`
    * Determine the status of centromeric contigs based on [`RepeatMasker`](https://www.repeatmasker.org/) annotations.

### Setup
```bash
pip install censtats
```

### Usage
```bash
usage: censtats [-h] {status} ...

Centromere statistics tool kit.

positional arguments:
  {status}

options:
  -h, --help  show this help message and exit
```

Read the docs [here](https://github.com/logsdon-lab/CenStats/wiki/Usage).

### Build
```bash
make venv && make build && make install
source venv/bin/activate && censtats -h
```

To run tests:
```bash
source venv/bin/activate && pip install pytest
pytest -s -vv
```
