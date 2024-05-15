# `CenStats`
[![CI](https://github.com/logsdon-lab/centromere-status-checker/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/centromere-status-checker/actions/workflows/main.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/CenStats)](https://pypi.org/project/CenStats/0.0.1/)

Centromere statistics toolkit.

* `status`
    * Determine the status of centromeric contigs based on [`RepeatMasker`](https://www.repeatmasker.org/) annotations.
* `length`
    * Estimate HOR array length from [`stv`](https://github.com/fedorrik/stv) bed file and [`HumAS-HMMER`](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) output.

### Setup
```bash
pip install censtats
```

### Usage
```bash
usage: censtats [-h] {status,length} ...

Centromere statistics toolkit.

positional arguments:
  {status,length}

options:
  -h, --help       show this help message and exit
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
