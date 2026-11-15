# MultiPieGen

[![PyPI](https://img.shields.io/pypi/v/multipiegen?label=pip%20install%20multipiegen)](https://pypi.org/project/multipiegen/)
[![Python](https://img.shields.io/pypi/pyversions/multipiegen)](https://pypi.org/project/multipiegen/)
[![License](https://img.shields.io/github/license/CMT-MU/MultiPieGen)](https://github.com/CMT-MU/MultiPieGen/blob/main/LICENSE)
[![docs](https://github.com/CMT-MU/MultiPieGen/actions/workflows/docs.yml/badge.svg)](https://cmt-mu.github.io/MultiPieGen/)

Generation code and data for the new multipole basis of [MultiPie](https://github.com/CMT-MU/MultiPie) Ver.3 and later.
MultiPie in this package and its documentation means MultiPie Ver.2 (whose general utilities are used here) unless
another version is named.

The new basis is built on a common root cluster (Oh and D6h clusters). Each point group and space group is related to it
through its parent group. This repository holds:

- the data tables (`multipiegen/data`): groups, Wyckoff positions, parent groups, operation and harmonics mappings, ...
- the code that generates the binary data and the PDFs from these tables

The input data are independent of MultiPie Ver.2: the binaries and PDFs are made from these tables alone. MultiPie Ver.2 is
used only for its general functions and utilities (expression parsing, `BinaryManager`, LaTeX/PDF output, tags and types),
and in the tests to compare with it.

The generated binaries and PDFs are delivered to other repositories: part of the binaries to MultiPie Ver.3 and later, and the rest
of the binaries and the PDFs to a separate repository.

## Requirements

- Python 3.12 or later
- [MultiPie](https://github.com/CMT-MU/MultiPie) Ver.2.4 or later, NumPy 2, SymPy
- for the PDFs only: a TeX distribution with `ptex2pdf` (e.g., TeX Live)

## Install

```bash
pip install multipiegen
```

The package includes the binary data for MultiPie (`multipiegen/binary_data`), which can be loaded with
`multipiegen.util.util_binary_data.load_binary_data`. The other binary data (`multipiegen/binary`) and the PDFs
(`multipiegen/pdf`) are kept in the repository but are **not** included in the package.

For development (to generate and load all the binaries), clone the repository and use the editable install:

```bash
pip install -e ".[dev]"
```

## Layout

```
multipiegen/        package
  data/             data tables (see the "Data tables" page)
  util/             utilities (multipolar harmonics, cluster basis, vector cluster basis, bond cluster basis, cluster SAMB, harmonics, harmonics symbols, multipolar harmonics, atomic multipole, info, group, symmetry operation, representation matrix, response tensor)
  create_binary/    code that generates the binary data
  create_pdf/       code that generates the PDFs
  binary/           generated binary data (info.pkl, cluster_scalar_basis.pkl, cluster_vector_basis.pkl, cluster_bond_basis.pkl, cluster_samb.pkl, harmonics.pkl, multipolar_harmonics.pkl, atomic_multipole.pkl, atomic_multipole_group.pkl, group.pkl, response_tensor.pkl) [not packaged]
  binary_data/      binary data for MultiPie (info.pkl, PG/, SG/, MPG/, MSG/: one file for each group, made from binary/) [packaged]
  pdf/              generated LaTeX and PDF (PG/032-Oh/harmonics.pdf, info/PG.pdf, ...) [not packaged]
docs/
  examples/         sample code for each feature
  src/guide.md      user guide: definitions of the provided quantities and their generation rules (guide_ja.md: Japanese)
  src/changes_from_multipie_v2.md  changes from MultiPie Ver.2, in the Details chapter (_ja.md: Japanese)
tests/              pytest tests
others/             non-package material
  creation/         scripts that generate the data tables in multipiegen/data
```

## Workflow

1. **Data tables**: the scripts in `others/creation/` write their output next to themselves. Copy the result into
   `multipiegen/data/` when it is final.
2. **Binary data**: run the modules in `multipiegen/create_binary/`. Each writes to `multipiegen/binary/`. Then
   `create_binary_data` regroups them into `multipiegen/binary_data/` for MultiPie.
3. **PDFs**: run the modules in `multipiegen/create_pdf/`. Each writes to `multipiegen/pdf/`.
4. **Check**: run `pytest`. The tests also check that the saved binaries agree with the current data tables.

When a feature is added, its generator goes in `create_binary/` or `create_pdf/`, and a sample goes in `docs/examples/`.

## Documentation

The documentation (this page, data, binary contents, harmonics PDFs, examples and API) is built with
[Jupyter Book](https://jupyterbook.org) from `docs/`:

```bash
make -C docs build   # -> docs/_build/html (needs the dev install)
```

The workflow `docs` builds it for pull requests and pushes, and publishes it to the `gh-pages` branch on `main`.

## Test

```bash
pytest
```

The tests run on GitHub Actions (ubuntu/macos × Python 3.12/3.13) for pull requests and for pushes to `main` and `develop`.

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

## License

[MPL-2.0](../LICENSE)
