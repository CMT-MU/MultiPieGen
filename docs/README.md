# MultiPieGen

[![test](https://github.com/hiro-k-phys/mp-newbasis/actions/workflows/test.yml/badge.svg)](https://github.com/hiro-k-phys/mp-newbasis/actions/workflows/test.yml)

Generation code and data for the new multipole basis of [MultiPie](https://github.com/CMT-MU/MultiPie).

The new basis is built on a common root cluster (Oh and D6h clusters). Each point group and space group is related to it
through its parent group. This repository holds:

- the data tables (`multipiegen/data`): groups, Wyckoff positions, parent groups, operation and harmonics mappings, ...
- the code that generates the binary data and the PDFs from these tables

The generated binaries and PDFs are delivered to other repositories: part of the binaries to MultiPie, and the rest of
the binaries and the PDFs to a separate repository.

## Requirements

- Python 3.12 or later
- [MultiPie](https://github.com/CMT-MU/MultiPie) 2.4 or later, NumPy 2, SymPy
- for the PDFs only: a TeX distribution with `ptex2pdf` (e.g., TeX Live)

## Install (development)

```bash
pip install -e ".[dev]"
```

The binary data (`multipiegen/binary_data`) and the PDFs (`multipiegen/pdf`) are kept in the repository but are
**not** included in the package; the package contains only `.py` files. Use the editable install above to load the binaries.

## Generated data

| output | content | generate | load |
|---|---|---|---|
| `binary_data/info.pkl` | overall information: group id/tag/name, character alias and compatibility relation, root cluster, orbital basis and basis functions, response tensor | `python -m multipiegen.create_binary.create_info` | `multipiegen.util.util_info.load_info()` |
| `binary_data/cluster_basis.pkl` | orthonormal cluster basis of the Oh (48) and D6h (24) clusters | `python -m multipiegen.create_binary.create_cluster_basis` | `multipiegen.util.util_cluster_basis.load_cluster_basis()` |
| `binary_data/group.pkl` | data of all PG (47), SG, MPG and MSG: symmetry operations, product table, mapping to the Oh/D6h cluster, character, Wyckoff site and bond, active multipoles | `python -m multipiegen.create_binary.create_group` | `multipiegen.util.util_group.load_group()` |
| `pdf/PG/{no}-{tag}/harmonics.{tex,pdf}` | harmonics of all 47 point groups, with polar and axial symbols side by side | `python -m multipiegen.create_pdf.create_pdf_harmonics` | |

The binaries use the format of MultiPie's `BinaryManager` (gzip-compressed pickle with a header). `str(load_info())`
shows the description of the contents.

```python
from multipiegen.util.util_info import load_info
from multipiegen.util.util_cluster_basis import load_cluster_basis

info = load_info()
info["harmonics"]["basis_function"]["fc"]  # (expression, <m|gamma>.T)
basis = load_cluster_basis()["Oh"]["basis"]  # {(rank, irrep, mul, comp): values on the 48 cluster points}
```

See [examples](examples/README.md) for runnable sample code for each feature.

## Layout

```
multipiegen/        package
  data/             data tables (see data/readme.md)
  util/             utilities (multipolar harmonics, cluster basis, harmonics symbols, info, group, symmetry operation)
  create_binary/    code that generates the binary data
  create_pdf/       code that generates the PDFs
  binary_data/      generated binary data (info.pkl, cluster_basis.pkl, group.pkl) [not packaged]
  pdf/              generated LaTeX and PDF (PG/032-Oh/harmonics.pdf, ...) [not packaged]
docs/
  examples/         sample code for each feature
tests/              pytest tests
others/             non-package material
  creation/         scripts that generate the data tables in multipiegen/data
  etc/              harmonics data generation (draft)
  notebooks/        notebooks for plotting and trials
  memo/             notes (TeX/PDF)
  claude/           verification and parent-SG pipeline
```

## Workflow

1. **Data tables**: the scripts in `others/creation/` write their output next to themselves. Copy the result into
   `multipiegen/data/` when it is final.
2. **Binary data**: run the modules in `multipiegen/create_binary/`. Each writes to `multipiegen/binary_data/`.
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
