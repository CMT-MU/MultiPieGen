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
  util/             utilities (multipolar harmonics, cluster basis)
  binary_data/      binary data (cluster_basis.pkl: orthonormal Oh/D6h cluster basis)
docs/examples/      sample code for each feature
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

## Branches

- Branch from `develop`, and open pull requests against `develop`.
- When a pull request is merged into `develop`, the workflow `delete merged branch` deletes its branch on GitHub.
- To delete the merged branch locally as well, enable the hook once per clone:

```bash
git config core.hooksPath .githooks
git config fetch.prune true
```

  After that, `git pull` on `develop` deletes the local branches that are merged and gone on the remote.
- `develop` is merged into `main` manually on GitHub; `develop` is never deleted.
