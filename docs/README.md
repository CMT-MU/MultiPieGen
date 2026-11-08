# MultiPieGen

[![test](https://github.com/hiro-k-phys/mp-newbasis/actions/workflows/test.yml/badge.svg)](https://github.com/hiro-k-phys/mp-newbasis/actions/workflows/test.yml)

Generation code and data for the new multipole basis of [MultiPie](https://github.com/CMT-MU/MultiPie).

## Install (development)

```bash
pip install -e ".[dev]"
```

## Layout

```
multipiegen/        package
  data/             generated data tables (see data/readme.md)
  util/             utilities (multipolar harmonics)
tests/              pytest tests (consistency of parent-group mappings)
others/             non-package material
  creation/         scripts that generate the data tables
  etc/              harmonics data generation (draft)
  notebooks/        notebooks for plotting and trials
memo/               notes (TeX/PDF)
claude/             verification and parent-SG pipeline
```

The generation scripts in `others/creation/` write their output next to themselves;
copy the result into `multipiegen/data/` when it is final.

## Test

```bash
pytest
```
