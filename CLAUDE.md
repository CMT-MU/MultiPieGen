# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MultiPieGen (`multipiegen`, repo `hiro-k-phys/MultiPieGen`) is the generator for the new multipole basis of
[MultiPie](https://github.com/CMT-MU/MultiPie) Ver.3 and later. It holds data tables (`multipiegen/data`) and code that turns them into
binaries (`multipiegen/binary/*.pkl`) and PDFs (`multipiegen/pdf/`). Those outputs are delivered elsewhere (part of
the binaries to MultiPie Ver.3 and later, the rest plus the PDFs to a separate repo), so keep binaries in MultiPie Ver.2's `BinaryManager`
format and readable without code from this repo. The installed `multipie` package is a runtime dependency for its general
functions and utilities only (`str_to_sympy`, `BinaryManager`, `PDFviaLaTeX`, tags, info/multipole types). **The input data
must be independent of MultiPie Ver.2:** never read its data (`Group`, its binaries or data modules) in `multipiegen/` or in
the table generators in `others/creation/`; take groups and symmetry operations from our own tables
(`util_group_data`, `util_symmetry_operation.group_symmetry_operation`). Tests may use `Group` to compare with
MultiPie Ver.2. A sibling checkout `../MultiPie` (tools in
`../MultiPie/others/tools/`) is the reference implementation that many generators here were ported from. MultiPie here means
Ver.2 (`multipie>=2.4`) unless another version is named; write "MultiPie Ver.2" when citing it in code, docs, or PRs.

## Commands

```bash
pip install -e ".[dev]"                              # editable install is required (see below)
pytest                                               # all tests (~20 s with -n auto; parent-group tests dominate)
pytest -n auto -p no:cacheprovider tests             # what CI runs
pytest tests/test_info.py::test_tesseral_tag         # single test
black -l 130 <files>                                 # formatter (line length 130, from pyproject)

python -m multipiegen.create_binary.create_info            # -> binary/info.pkl
python -m multipiegen.create_binary.create_cluster_scalar_basis   # -> binary/cluster_scalar_basis.pkl
python -m multipiegen.create_binary.create_cluster_vector_basis  # -> binary/cluster_vector_basis.pkl (s=1, intermediate for cluster_bond_basis, ~1 min)
python -m multipiegen.create_binary.create_cluster_bond_basis    # -> binary/cluster_bond_basis.pkl (needs cluster_vector_basis.pkl)
python -m multipiegen.create_binary.create_cluster_samb          # -> binary/cluster_samb.pkl (site/bond of real PG + SG, bond = Q + i×(T/M); needs group, cluster_scalar_basis, cluster_bond_basis)
python -m multipiegen.create_binary.create_harmonics            # -> binary/harmonics.pkl (expression and u-matrix of the PG harmonics, l<=11; ~10 s)
python -m multipiegen.create_binary.create_multipolar_harmonics  # -> binary/multipolar_harmonics.pkl (s=1..3, l<=3, 47 PGs, parallel ~2.5 min)
python -m multipiegen.create_binary.create_atomic_multipole        # -> binary/atomic_multipole.pkl (l<=7, s=0,1, LMs/JML/LM; ~10 s)
python -m multipiegen.create_binary.create_atomic_multipole_group  # -> binary/atomic_multipole_group.pkl (47 PGs, needs atomic_multipole.pkl; ~30 s)
python -m multipiegen.create_binary.create_group           # -> binary/group.pkl (parallel ~35 s; --serial ~2 min)
python -m multipiegen.create_binary.create_response_tensor # -> binary/response_tensor.pkl (needs group.pkl; ~20 s)
python -m multipiegen.create_pdf.create_pdf_harmonics      # -> pdf/PG/{no}-{tag}/harmonics.{tex,pdf} (needs ptex2pdf, ~1.5 min)
python -m multipiegen.create_pdf.create_pdf_group         # -> pdf/{PG,SG,MPG,MSG}/*/{character_table,symmetry_operation,wyckoff_*,representation_matrix,multipolar_harmonics_s*_*,atomic_multipole_*,response_tensor_X}.{tex,pdf} (~5200 PDFs, parallel; --group/--kind to limit)
python -m multipiegen.create_pdf.create_pdf_info          # -> pdf/info/{PG,SG,MPG,MSG,compatibility_relation,active_multipole}.{tex,pdf} (~30 s)
python docs/examples/<name>.py                             # runnable samples
make -C docs build                                         # Jupyter Book -> docs/_build/html (gh-pages on main)
```

`binary/` and `pdf/` are in git but deliberately **not packaged** (`include-package-data = false`; the wheel is
`.py` only). Anything that loads a binary needs the editable install, and CI uses `pip install -e .`. Do not re-add
package-data.

Regenerating a binary or PDF whose content is unchanged still produces a git diff (gzip/PDF timestamps). Compare the
content, and `git restore` the file if only the timestamp changed. `tests/test_*` include "binary is up to date"
checks that rebuild in memory and compare with the saved `.pkl`, so after changing a data table, rerun the matching
`create_binary` module.

## Architecture

**Top groups and mapping.** All harmonics are defined once for the two top groups, Oh and D6h
(`data_top_harmonics`; `util/util_top_harmonics.top_harmonics_symmetry` / `top_harmonics_u_matrix`). Real and complex combinations are
mixed there: complex irreps end in `a`/`b`, e.g. `Ega`, `T1ua`. Every one of the 47 point groups
(32 + 5 settings `-1` + 10 complex `-c`) gets its harmonics through `data_harmonics_mapping[pg] = {"top", "Q": {key:
(coeff, top_key)}, "G": {...}}`. Polar `Q` and axial `G` of a group come from the same top harmonic, so the axial
partner of a polar harmonic is found by matching `top_key`, with sign `cQ*cG`
(`util/util_harmonics_symbol.polar_to_axial`). `data_compatibility_relation` alone is not enough for subgroups: the
mapping depends on the parity of the top harmonic.

**Key conventions.**
- Harmonics key: `(rank, irrep, multiplicity, component)`; multiplicity is 1-based and component is 0-based, with `-1`
  meaning unique.
- Tesseral key: `(l, m)` with `m < 0` meaning `S`.
- Binaries keyed by `PGMultipoleType` do not hold the xyz expression or u-matrix of the harmonics. Look them up by
  `(X, l, Gamma, n)` of the key in `harmonics.pkl` (`util_harmonics.harmonics_expression` / `harmonics_u_matrix`),
  with T -> Q and M -> G.
- Symmetry-operation mappings in `data_parent_pg/sg` and `util/util_op_mapping` are 1-based.
- `util_op_mapping.op_mapping_d6h_oh_cluster` (PG only) maps to the Oh/D6h cluster. An SG reaches the cluster
  through its PG: its operations have the PG's point parts in the same order. `op_mapping_parent_sg` maps an SG to a
  top parent SG (221-230, 191-194) with changes of axes, so its point parts can differ. Do not use it to index the
  cluster.
- The irrep order in `data_character_table` is authoritative for sorting. g before u, and `'` block before `''` block
  (C3h, D3h, D3h-1, C3h-c were reordered this way). This order differs from MultiPie Ver.2, which sorts irreps by name string.

**Orbitals.** `data_basis.py` replaces MultiPie Ver.2's `data_definition.py`. Orbital tags were renamed relative to MultiPie Ver.2:
f2→fc, f1→fs, f3x→fx, f3y→fy, ga→gc, gb→gs, gc→gc4, i→i3, l→l3, i0→ixy, j→j3. Note that `gc` now means old `ga`.
Orbitals carry no expressions; `harmonics_key[name]` points into the top harmonics. Intentional sign changes versus MultiPie Ver.2:
`gv` and D6h rank-4 `E2g` multiplicity 2 are multiplied by −1. Both carry into D6h-family subgroups and the response
tensor.

**Root cluster.** Oh (48 points) and D6h (24 points) clusters are written in Cartesian coordinates with `x,y,z`
placeholders (`data_root_cluster`). `site0` and `vector0` are substituted into them. The D6h representative point is
the Oh one with the principal axis rotated to [111]. `util/util_cluster_basis` builds the exact (sympy)
Gram–Schmidt basis per (irrep, component), from lower rank, and must yield 48 / 24 vectors (regular representation).

**Where code goes.**
- Binary generators go in `multipiegen/create_binary/`, PDF generators in `multipiegen/create_pdf/`. Run them as
  `python -m`, with an `if __name__ == "__main__"` guard. Loaders and reusable logic go in `multipiegen/util/`.
- `multipiegen/data` holds only source tables. Data derivable from them are computed in `multipiegen/util`, not stored:
  `util_group_data` (tags, ids, id sets, group names), `util_top_harmonics`, `util_harmonics_mapping`
  (`data_harmonics_mapping`), `util_op_mapping`, and `create_group.character_plus`. Keep it that way: do not add a
  table that can be computed from the others.
- The remaining `data_*.py` table generators live in `others/creation/`. They write next to themselves; copy the result into
  `multipiegen/data/` when final.
- **Never overwrite `others/creation/data_*.py` or anything in `multipiegen/data/` without asking the maintainer
  first.** This includes running an `others/creation/create_*.py` script, copying a result into `multipiegen/data/`,
  and editing or reformatting a data table, even when only the order or formatting changes. These tables are
  curated by hand; regenerating them can silently undo manual fixes.
- `others/` holds non-package material: `creation/` (table generators). `others/archive/` keeps old material locally
  only (old memo, etc, notebooks, earlier verification pipelines, the MultiPie Ver.2 comparison script, old scripts);
  it is in `.gitignore` and will be deleted. Do not add files there to git, and do not depend on them.
- `docs/src/notes/` holds the notes on the multipole data (generation rules and changes from MultiPie Ver.2, in
  English and Japanese), the "Notes" chapter of the docs. Update them when the algorithms or keys change.
- Every new feature also gets a runnable sample `docs/examples/<feature>.py` and a row in `docs/examples/README.md`
  (the docs page `docs/src/examples.md` is generated from that table by `docs/create_pages.py`).

## Branch flow

The default branch is `main`, but work is cut from `origin/develop`, and every PR targets `develop`
(`gh pr create --base develop`). develop→main PRs and merges are done manually by the maintainer on GitHub; never push
to main. Repo-level auto-delete is off so that develop survives develop→main merges. Instead,
`.github/workflows/delete-merged-branch.yml` deletes head branches merged into develop, and
`.githooks/post-merge` (enable with `git config core.hooksPath .githooks`) deletes local merged branches on
`git pull` in develop.
