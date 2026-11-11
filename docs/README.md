# MultiPieGen

[![test](https://github.com/hiro-k-phys/MultiPieGen/actions/workflows/test.yml/badge.svg)](https://github.com/hiro-k-phys/MultiPieGen/actions/workflows/test.yml)

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

## Install (development)

```bash
pip install -e ".[dev]"
```

The binary data (`multipiegen/binary`) and the PDFs (`multipiegen/pdf`) are kept in the repository but are
**not** included in the package; the package contains only `.py` files. Use the editable install above to load the binaries.

## Generated data

| output | content | generate | load |
|---|---|---|---|
| `binary/info.pkl` | overall information: group id/tag/name, character alias and compatibility relation, root cluster, orbital basis and basis functions, response tensor | `python -m multipiegen.create_binary.create_info` | `multipiegen.util.util_info.load_info()` |
| `binary/cluster_scalar_basis.pkl` | orthonormal (scalar) cluster basis of the Oh (48) and D6h (24) clusters, keyed by PGMultipoleType as MultiPie Ver.2's site SAMB | `python -m multipiegen.create_binary.create_cluster_scalar_basis` | `multipiegen.util.util_cluster_basis.load_cluster_basis()` |
| `binary/cluster_vector_basis.pkl` | orthonormal vector cluster basis (internal rank s=1, polar q / axial g internal basis) of the Oh (144) and D6h (72) clusters, keyed by PGMultipoleType as MultiPie Ver.2's vector SAMB; intermediate basis for cluster_bond_basis.pkl | `python -m multipiegen.create_binary.create_cluster_vector_basis` | `multipiegen.util.util_cluster_vector_basis.load_cluster_vector_basis()` |
| `binary/cluster_bond_basis.pkl` | orthonormal bond cluster basis of the Oh (48) and D6h (24) clusters, from the vector cluster basis and the bond vectors (vector0), keyed by PGMultipoleType as MultiPie Ver.2's bond_a | `python -m multipiegen.create_binary.create_cluster_bond_basis` | `multipiegen.util.util_bond_cluster_basis.load_bond_cluster_basis()` |
| `binary/response_tensor.pkl` | response tensors (rank 0-4, Q/G/T/M) of all MPG in terms of active multipoles | `python -m multipiegen.create_binary.create_response_tensor` | `multipiegen.util.util_response_tensor.load_response_tensor()` |
| `binary/atomic_multipole.pkl` | matrices of atomic multipoles (rank l <= 7, s = 0, 1) in the LMs/JML/LM basis of s, p, d, f orbitals | `python -m multipiegen.create_binary.create_atomic_multipole` | `multipiegen.util.util_atomic_multipole.load_atomic_multipole()` |
| `binary/atomic_multipole_group.pkl` | atomic multipole SAMB of all 47 point groups in the JML/Lgs/Lg basis | `python -m multipiegen.create_binary.create_atomic_multipole_group` | `multipiegen.util.util_atomic_multipole.load_atomic_multipole_group()` |
| `binary/cluster_samb.pkl` | cluster SAMB of each Wyckoff site and bond (site; bond = symmetric Q and anti-symmetric T/M) of the real point groups and the space groups | `python -m multipiegen.create_binary.create_cluster_samb` | `multipiegen.util.util_cluster_samb.load_cluster_samb()` |
| `binary/harmonics.pkl` | expression and u-matrix of the harmonics of all 47 point groups (rank l <= 11); the other multipole-keyed binaries look them up by the key | `python -m multipiegen.create_binary.create_harmonics` | `multipiegen.util.util_harmonics.load_harmonics()` |
| `binary/multipolar_harmonics.pkl` | multipolar harmonics (s = 1, 2, 3, rank l <= 3) of all 47 point groups with polar (q) and axial (g) internal basis (data_basis orbitals) | `python -m multipiegen.create_binary.create_multipolar_harmonics` | `multipiegen.util.util_harmonics_multipole.load_multipolar_harmonics()` |
| `binary/group.pkl` | data of all PG (47), SG, MPG and MSG: symmetry operations, product table, mapping to the Oh/D6h cluster, character, Wyckoff site and bond, active multipoles | `python -m multipiegen.create_binary.create_group` | `multipiegen.util.util_group.load_group()` |
| `pdf/PG/{no}-{tag}/harmonics.{tex,pdf}` | harmonics of all 47 point groups, with polar and axial symbols side by side | `python -m multipiegen.create_pdf.create_pdf_harmonics` | |
| `pdf/{PG,SG,MPG,MSG}/*/{character_table,symmetry_operation,wyckoff_site,wyckoff_bond,representation_matrix,multipolar_harmonics_s{1,2,3}_{q,g},atomic_multipole_{lg,lgs,jml},response_tensor_X}.{tex,pdf}` | each group: character table, representation matrices, multipolar harmonics and atomic multipoles (PG), symmetry operations, Wyckoff site and bond (PG, SG), response tensors (MPG) | `python -m multipiegen.create_pdf.create_pdf_group` | |
| `pdf/info/*.{tex,pdf}` | overall information: lists of PG, SG, MPG and MSG, compatibility relation with orbitals, active multipoles | `python -m multipiegen.create_pdf.create_pdf_info` | |

The binaries use the format of MultiPie Ver.2's `BinaryManager` (gzip-compressed pickle with a header). `str(load_info())`
shows the description of the contents.

```python
from multipiegen.util.util_info import load_info
from multipiegen.util.util_cluster_basis import load_cluster_basis

info = load_info()
info["harmonics"]["basis_function"]["fc"]  # (expression, <m|gamma>.T)
basis = load_cluster_basis()["Oh"]["basis"]  # {PGMultipoleType: values (dim, 48) on the cluster points}
```

See [examples](examples/README.md) for runnable sample code for each feature.

## Layout

```
multipiegen/        package
  data/             data tables (see the "Data tables" page)
  util/             utilities (multipolar harmonics, cluster basis, vector cluster basis, bond cluster basis, cluster SAMB, harmonics, harmonics symbols, multipolar harmonics, atomic multipole, info, group, symmetry operation, representation matrix, response tensor)
  create_binary/    code that generates the binary data
  create_pdf/       code that generates the PDFs
  binary/           generated binary data (info.pkl, cluster_scalar_basis.pkl, cluster_vector_basis.pkl, cluster_bond_basis.pkl, cluster_samb.pkl, harmonics.pkl, multipolar_harmonics.pkl, atomic_multipole.pkl, atomic_multipole_group.pkl, group.pkl, response_tensor.pkl) [not packaged]
  pdf/              generated LaTeX and PDF (PG/032-Oh/harmonics.pdf, info/PG.pdf, ...) [not packaged]
docs/
  examples/         sample code for each feature
  src/notes/        notes on the multipole data (generation rules, changes from MultiPie Ver.2; English and Japanese)
tests/              pytest tests
others/             non-package material
  creation/         scripts that generate the data tables in multipiegen/data
```

## Workflow

1. **Data tables**: the scripts in `others/creation/` write their output next to themselves. Copy the result into
   `multipiegen/data/` when it is final.
2. **Binary data**: run the modules in `multipiegen/create_binary/`. Each writes to `multipiegen/binary/`.
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
