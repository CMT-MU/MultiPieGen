---
orphan: true
---

# Changes from MultiPie Ver.2 (multipole data)

This note lists what MultiPieGen does differently from MultiPie Ver.2 for the multipole-keyed data: conventions,
algorithms, and the resulting differences. The algorithms of MultiPieGen themselves are described in
[Generation of the multipole data](generation_rules.md). The numbers below come from a comparison with MultiPie Ver.2
made on 2026-11-11 (develop 369c8bf). A Japanese version is [here](changes_from_multipie_v2_ja.md).

Contents:

1. Conventions (keys, irreps, orbital names, signs)
2. Harmonics of the point groups
3. Multipolar harmonics
4. Cluster bases (new)
5. Cluster SAMB
6. Atomic multipoles
7. Symmetry-operation mapping used by the cluster SAMB
8. Other data (representation matrices, active multipoles, response tensors)

---

## 1. Conventions

| item | MultiPie Ver.2 | MultiPieGen |
|---|---|---|
| multiplicity `n` | renumbered in each data set (multipolar harmonics, atomic multipole group, cluster SAMB) | always the multiplicity of the harmonics `(X, l, Gamma)` of the point group, which fixes its xyz polynomial; duplicates are numbered by `p` |
| irrep order | sorted by the irrep name | the order of `data_character_table` (`g` before `u`, `'` block before `''` block) |
| orbital names | f2, f1, f3x, f3y, ga, gb, gc, i, l, i0, j | fc, fs, fx, fy, gc, gs, gc4, i3, l3, ixy, j3 (`gc` now means the old `ga`) |
| orbital definitions | `data_definition.py` | `data_basis.py`; each orbital points to a harmonics key, expressions come from the top harmonics |
| sign of D6h rank-4 `E2g` multiplicity 2 | as in MultiPie | multiplied by -1 (both components) |
| sign of the orbital `gv` | as in MultiPie | multiplied by -1 (`data_response_tensor`) |
| expression of the harmonics in the values | the multipole-keyed binaries hold the cartesian expression `ex` of the harmonics | not held; looked up by `(X, l, Gamma, n)` of the key in `harmonics.pkl` (values: multipolar harmonics `(internal_coefficient, multipole_expression)`, atomic multipole group `matrix`, cluster SAMB `samb`) |

- The `E2g` sign change is a change of the top harmonics. It carries into the D6h-family subgroups (rank-4 E
  multiplets of 16 hexagonal/trigonal groups) and into everything built on the harmonics.
- The `gv` sign change is a change of the orbital (and of the response-tensor multipole); the top harmonics of Oh are
  unchanged, and the harmonics of the Oh family agree with MultiPie Ver.2.
- In the data compared so far, the harmonics multiplicity `n` happened to coincide with MultiPie Ver.2's renumbered
  `n` (multipolar harmonics, atomic multipole group). The difference of rule matters for the bond_a cluster SAMB.

## 2. Harmonics of the point groups

- Defined once for Oh and D6h (`data_top_harmonics`, `util_top_harmonics`) and mapped to the 47 point groups by
  `data_harmonics_mapping` (`c = ±1` times a top harmonics). MultiPie Ver.2 stores the harmonics of each point group
  in `harmonics.pkl`.
- Result: the keys `(X, l, Gamma, n)` are identical for all 47 groups (rank <= 11). 32 multiplets have the opposite
  sign, all from the `E2g` sign change (section 1).

## 3. Multipolar harmonics

| item | MultiPie Ver.2 (`harmonics_multipole.pkl`) | MultiPieGen (`multipolar_harmonics.pkl`) |
|---|---|---|
| internal basis | tesseral basis for all groups (`Q_v`, `Q_xy`, …, old f names) | cubic orbitals (`cubic_orbital_basis[s]`) for the Oh family, hexagonal orbitals (`hexagonal_orbital_basis[s]`) for the D6h family |
| internal symbols | `Q_v`, `Q_2`, … | `Q + orbital_latex_name`, e.g. `Q_{u}`, `Q_{x}^{\alpha}` |
| `n` | renumbered per `(x, l, X, Gamma, k)` | harmonics multiplicity; `p` for duplicates |
| rank | l <= 11 | l <= 3 |

- Result: the keys are identical (14,302 entries, l <= 3). After converting the internal-basis coefficients to the
  tesseral basis, all entries are identical, so only the choice of the internal basis (orbitals, order, symbols)
  differs.

## 4. Cluster bases (new)

`cluster_scalar_basis.pkl`, `cluster_vector_basis.pkl`, and `cluster_bond_basis.pkl` have no counterpart in MultiPie
Ver.2. The closest one is `harmonics_root_cluster.pkl`.

| item | MultiPie Ver.2 `harmonics_root_cluster` | MultiPieGen cluster bases |
|---|---|---|
| cluster | one root cluster of 64 points for Oh and D6h | separate Oh (48 points) and D6h (24 points) clusters |
| representative point | MultiPie Ver.2's point | site0 = (1, 5, 6), vector0 = unit vector of (6, 1, 5) |
| content | values of the harmonics of each point group up to rank 11 (not orthonormalized, dependent ones included) | Gram-Schmidt per (irrep, component) of the top group from lower rank, dependent ones dropped: an orthonormal basis (regular representation) |
| internal type | x = q | x = q and x = g (vector and bond bases) |

- The key format follows MultiPie Ver.2's site SAMB (scalar), vector SAMB (vector), and bond_a (bond).
- The bond basis is built by contracting the vector basis with the bond vectors, then Gram-Schmidt in the order rank,
  k, multiplicity. MultiPie Ver.2 keys each anti-symmetric bond function by the first overlapping vector SAMB; that
  rule is not used.

## 5. Cluster SAMB

| item | MultiPie Ver.2 (`cluster_samb.pkl`) | MultiPieGen (`cluster_samb.pkl`) |
|---|---|---|
| source | harmonics of each point group evaluated on the 64-point root cluster | the cluster bases of Oh/D6h restricted to the operations of the group (coset 0 of `op_mapping_d6h_oh_cluster`) |
| special points | sum over the Wyckoff mapping, then orthogonalization | the same (average of the general-point SAMB by the Wyckoff mapping, sign for bond_a) |
| site, bond_s | Gram-Schmidt in [Gamma, l, Q/G, n, p] | Gram-Schmidt in [Gamma, l, n, p] per irrep; the same keys |
| bond_a | orthogonalized in the monopole order, then keyed by the first overlapping vector SAMB, `p` renumbered per Wyckoff | `cluster_bond_basis` (x = q) restricted, Gram-Schmidt in [Gamma, l, k, Q/G, n, p], keys of the top group kept (including `p`) |
| vector SAMB | made | not made |
| bond kinds | `"bond_s"` and `"bond_a"` separately | `"bond"` = `bond_s` (Q) + i `bond_a` with Q -> T, G -> M |
| groups | real PG and SG | real PG and SG (complex PG not included) |

- Result: for every Wyckoff site and bond, each irrep spans the same space as MultiPie Ver.2 (1889 sites, 6663 bonds).
  The site and bond_s keys are identical.
- bond_a keys differ in 5990 of 6663 bond Wyckoff positions (keys before the change of type to T/M, i.e. X = Q/G):
  - p only (4866): the top key's `p` is kept (e.g. Oh `3a@1o`: Q₁(T1u) p = 1 here, p = -1 in MultiPie Ver.2).
  - n (93): another harmonics multiplicity of the same (X, l, Gamma) is used (e.g. C4 `4c@1a`: Q₂(B) n = 1 here,
    n = 2 in MultiPie Ver.2).
  - (X, l, Gamma) (1031, 168 groups including 15 point groups): another key for the same space, often of lower rank
    here (e.g. C4v `4f@1a`: G₁(E) here, Q₁(E) in MultiPie Ver.2).
- The values differ because the representative general point (site0) differs.

## 6. Atomic multipoles

- `atomic_multipole.pkl`: identical to MultiPie Ver.2 (keys, order, values). Only the computation differs: only the
  matrix elements allowed by Jz conservation are evaluated (about 10 s).
- `atomic_multipole_group.pkl`:

| item | MultiPie Ver.2 | MultiPieGen |
|---|---|---|
| construction | for each point group, project the atomic multipoles with its harmonics, Gram-Schmidt per (X, irrep), times sqrt(dim) and 1/sqrt(2) for off-diagonal blocks | the normalized SAMB of each top harmonics component, computed once per top group; each point group takes `c` times it (equivalent normalization) |
| `n` | renumbered | harmonics multiplicity |
| groups | 37 real PG | 47 PG (including complex `-c`) |
| Lg basis | MultiPie Ver.2's lg orbitals (old names) | `data_basis.atomic_lgs_basis` |

- Result (37 real PG): the keys are identical; 896 matrices have the opposite sign (the `E2g` sign change); the other
  60,757 matrices are identical.

## 7. Symmetry-operation mapping used by the cluster SAMB

| item | MultiPie Ver.2 | MultiPieGen |
|---|---|---|
| `symmetry_operation["mapping"]` of PG | right cosets `H g` from the product table | left cosets `g H` along the parent chain (`util_op_mapping.op_mapping_d6h_oh_cluster`) |
| of SG | none | the mapping of its point group (the operations of an SG have the point parts of its PG in the same order) |
| SG to parent SG | none | `op_mapping_parent_sg` (with changes of axes; not used for the cluster) |

## 8. Other data

- Representation matrices: identical for all 241 irreps of the 47 point groups.
- Products of irreps in the character data of `group.pkl` (polar-axial conversion, symmetric and anti-symmetric
  products): computed from the characters (`create_group.character_plus`) instead of being stored. Two differences
  from MultiPie Ver.2:
  - Corrected: for E of C4 and S4 and Eg, Eu of C4h (real 2D irreps of complex type), MultiPie Ver.2 has the symmetric
    product 2A + B and the anti-symmetric product B. The correct values are A + 2B and A: for the basis (x, y),
    x² + y² is A, x² − y² and xy are B, and x₁y₂ − y₁x₂ is invariant under C4. (E = Ea + Eb with Ea⊗Ea = Eb⊗Eb = B and
    Ea⊗Eb = A.)
  - The order of the polar-axial conversion follows the character table (10 groups had another order).
- Active multipoles of MPG (`group.pkl`): identical after renaming the orbitals (51 of 157 groups have renamed
  orbitals, e.g. `Qf1` → `Qfs`).
- Response tensors (`response_tensor.pkl`): identical to MultiPie Ver.2's `Group.response_tensor_all` except for the
  multipole names (checked when the binary was made, not in the comparison of 2026-11-11).
