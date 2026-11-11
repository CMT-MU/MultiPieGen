# Generation of the multipole data: algorithms, keys, and orbital basis

This note collects, in one place, how MultiPieGen builds its multipole-keyed data, and the rules for keys and orbital
bases. It describes the algorithms of this package only. Each section names the code that does the work. A Japanese version is `generation_rules_ja.md`. The changes
from MultiPie Ver.2 are collected separately in `changes_from_multipie_v2.md`.

Contents:

1. Common conventions (keys, irreps, top groups, orbital names)
2. Harmonics of the point groups
3. Multipolar harmonics (`multipolar_harmonics.pkl`)
4. Cluster bases of the top groups (`cluster_scalar_basis.pkl`, `cluster_vector_basis.pkl`, `cluster_bond_basis.pkl`)
5. Cluster SAMB of each Wyckoff site and bond (`cluster_samb.pkl`)
6. Atomic multipoles (`atomic_multipole.pkl`) and their point-group SAMB (`atomic_multipole_group.pkl`)
7. Exact arithmetic and simplification

---

## 1. Common conventions

### 1.1 Multipole key

All multipoles are keyed by `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` of the `multipie` package. `SphericalMultipoleType(X, l, s, k,
x)` is used for the atomic multipoles before symmetry adaptation.

| field | meaning |
|---|---|
| `X` | type: `Q` (electric), `G` (electric toroidal), `T` (magnetic toroidal), `M` (magnetic). The harmonics and the cluster data use `Q`/`G` only. |
| `l` | rank |
| `Gamma` | irrep of the point group (name of `data_character_table`) |
| `n` | multiplicity of the harmonics `(X, l, Gamma)` of the point group, **1-based**, `-1` if unique. It always identifies the xyz polynomial of the harmonics (never renumbered). |
| `p` | additional multiplicity, used only to number entries that share `(X, l, Gamma, n)` and the other fields; `-1` if unique |
| `s` | internal rank (0 for scalar functions) |
| `k` | internal component, `-s <= k <= s` |
| `x` | internal type, `q` (polar internal basis) or `g` (axial internal basis) |

All components of a multiplet are stored together in the value, as rows in the order of the components.

Other keys:

- Harmonics key `(l, irrep, n, comp)`: `n` 1-based (`-1` unique), `comp` 0-based (`-1` for a single component).
- Tesseral key `(l, m)`: `m >= 0` is the cosine-type `C_{l,m}`, `m < 0` the sine-type `S_{l,|m|}`.
- Symmetry-operation numbers in mappings are 1-based.

### 1.2 Irreps

- The order of irreps in `data_character_table` is authoritative for sorting: `g` before `u`, and the `'` block
  before the `''` block (C3h, D3h, D3h-1, C3h-c are in this order).
- Complex irreps (the `-c` point groups) end in `a`/`b`, e.g. `Ega`, `T1ua`. Real 2D irreps of complex type (e.g. `E`
  of C3, C4, C6, S4, C3i, C4h, C6h, T, Th) are used in the real point groups.

### 1.3 Top groups

- All harmonics are defined once for the two top groups, Oh and D6h, and mapped to every one of the 47 point groups
  (32 crystallographic, 5 other settings `-1`, 10 complex `-c`).
- Cubic point groups use Oh. Trigonal and hexagonal point groups use D6h. Lower-symmetry groups use the top group
  that their parent chain reaches (`data_harmonics_mapping[pg]["top"]`).

### 1.4 Orbital names

`data_basis` is the source of the orbital names (`data_harmonics` holds tables with older names; do not use them for
names).

Each orbital name points to a harmonics key (`data_basis.harmonics_key`): a tesseral key, or an Oh/D6h top key.
Its expression and u-matrix are `info["harmonics"]["basis_function"][name] = (expression, u-matrix)`, where the
u-matrix is the column `<m|orbital>` for `m = l, l-1, ..., -l`.

| list in `data_basis` | orbitals (rank 0..3) | used for |
|---|---|---|
| `atomic_orbital_basis` (tesseral order) | s; px, py, pz; dv, dxy, dxz, dyz, du; fc, fs, fbz, f3, fx, fy, faz | atomic (Lg/Lgs) basis |
| `cubic_orbital_basis` | s; px, py, pz; du, dv, dyz, dxz, dxy; f3, fax, fay, faz, fbx, fby, fbz | internal basis of cubic groups (Oh family) |
| `hexagonal_orbital_basis` | s; px, py, pz; du, dxz, dyz, dxy, dv; faz, fs, fc, fx, fy, f3, fbz | internal basis of trigonal/hexagonal groups (D6h family) |

LaTeX names: `atomic_latex_name` (e.g. `p_{x}`, `f_{z}^{\alpha}`) for the atomic basis, and `orbital_latex_name`
(e.g. `_{x}`, `_{3c}`) with the prefix Q/G for the internal basis symbols of the multipolar harmonics
(`info["harmonics"]["internal_basis"]`).

Spinful bases (`data_basis`): spin `["u", "d"]`; LM `M = L..-L`; JM `J` ascending, `M` descending; LMs and Lgs are
`(m or orbital, spin)` with spin fastest.

---

## 2. Harmonics of the point groups

Data: `data_top_harmonics_info` (top harmonics), `data_harmonics_mapping` (mapping to each point group).
Code: `util_harmonics_multipole.pg_harmonics`, and the table generators `others/creation/create_top_harm.py`,
`create_harmonics_mapping.py`.

1. **Top harmonics.** For Oh and D6h, every harmonics up to rank 11 is given as a cartesian polynomial
   (`top_harmonics_symmetry[top][top_key]`) and as its u-matrix to the spherical harmonics
   (`top_harmonics_u_matrix[top][top_key]`, `m = l..-l`). Only polar harmonics are stored; the axial one of the same
   function has the irrep with `g` and `u` interchanged.
   - Real and complex combinations are mixed: complex irreps end in `a`/`b` and have a single component.
   - The basis of each 2D irrep is chosen so that the simplest `(x, y)` basis gives representation matrices of `C_n`
     equal to rotation matrices.
2. **Mapping.** `data_harmonics_mapping[pg] = {"top": top, "Q": {pg_key: (c, top_key)}, "G": {...}}` with `c = ±1`.
   A harmonics of the point group is `c` times a top harmonics, `X(pg_key) = c T(top_key)`, and its u-matrix is
   `c U_T`.
   - It is created from the compatibility relation (`data_compatibility_relation`): each top harmonics component
     `(irrep, comp)` goes to a polar irrep/component and an axial irrep/component of the subgroup with a coefficient.
     For a complex subgroup irrep, the complex top harmonics (`irrep + "a"/"b"`) is referred to instead.
   - The multiplicity `n` of the subgroup is assigned in the order of `top_harmonics_tag` (`-1` if unique).
   - Polar `Q` and axial `G` of a group come from the same top harmonics; the axial partner of a polar harmonics is
     found by matching `top_key`, with the sign `cQ cG` (`util_harmonics_symbol.polar_to_axial`). The compatibility
     relation alone is not enough for subgroups, because the mapping depends on the parity of the top harmonics.
3. **Representation matrices** follow from the u-matrices (`util_representation_matrix`); all harmonics of an irrep
   give the same matrices (tested up to rank 11).
4. **Binary (`harmonics.pkl`).** `{PG_id: {PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q"): (expression (dim),
   u_matrix (2l+1, dim))}}` for the 47 point groups, rank <= 11, in the order Gamma, l, Q/G, n (`util_harmonics`,
   `create_binary/create_harmonics.py`). The other binaries keyed by `PGMultipoleType` do not hold the expression or
   the u-matrix of the harmonics; `util_harmonics.harmonics_expression(pg, key)` and `harmonics_u_matrix(pg, key)`
   look them up by `(X, l, Gamma, n)` of any key, with T → Q and M → G (for a space group, use its point
   group).

---

## 3. Multipolar harmonics (`multipolar_harmonics.pkl`)

Code: `util_harmonics_multipole.create_multipolar_harmonics`, `util_multipolar_harmonics` (Olm),
`create_binary/create_multipolar_harmonics.py`. Range: s = 1, 2, 3; l <= 3; 47 point groups.

1. **Multipolar spherical harmonics.** `O_{lm}^{(s)}(k)` (`util_multipolar_harmonics.Olm`) couples the spherical
   harmonics of orbital rank `l+k` with an internal basis of rank `s`:
   `O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} sqrt(2l+1) sum_n (l+k, l, s; m-n, -m, n) O_{l+k, m-n} e_n`, where
   `O_{lm} = sqrt(4 pi/(2l+1)) r^l Y_lm` and `e_n` are the internal basis vectors, expressed in the chosen internal
   basis (`InternalBasis`).
2. **Internal basis.** The orbitals of `data_basis.cubic_orbital_basis[s]` (top Oh) or `hexagonal_orbital_basis[s]`
   (top D6h), with their u-matrices from `info["harmonics"]["basis_function"]`.
3. **Projection.** For every harmonics `(X, l, Gamma, n)` of the group with u-matrix `U` (section 2) and every k,
   the internal coefficient `U^T O_l^{(s)}(k)`, a `(dim, 2s+1)` array of polynomials in x, y, z (each element
   factorized): the coefficient of each internal basis function. Vanishing ones are dropped (`w3j = 0`).
4. **Internal type.** `x = q` if (`X = Q` and `s+k` even) or (`X = G` and `s+k` odd), otherwise `x = g`.
5. **Expression.** multipole expression `= internal coefficient @ symbols`, with the internal symbols `Q/G + orbital_latex_name`
   (`info["harmonics"]["internal_basis"][crystal][(s, x)]`).
6. **Key.** `PGMultipoleType(X, l, Gamma, n, p, s, k, x)`: `n` is the harmonics multiplicity (not renumbered); `p`
   numbers the entries with the same `(x, l, X, Gamma, n, k)` (`-1` if unique).
7. **Order.** `x` (q, g), `s`, `Gamma`, `l`, `X` (Q, G), `k`, `n`, `p`.
8. **Value.** `(internal_coefficient (dim, 2s+1), multipole_expression (dim))`. The expression of the harmonics
   `(X, l, Gamma, n)` itself is in `harmonics.pkl`.

---

## 4. Cluster bases of the top groups

Code: `util_cluster_basis`, `util_cluster_vector_basis`, `util_bond_cluster_basis` and the matching
`create_binary` modules. Data: `data_root_cluster`, `info["root_cluster"]`.

### 4.1 Clusters

- The Oh cluster has 48 points and the D6h cluster 24 points, written with `x, y, z` placeholders in cartesian
  coordinates (`data_root_cluster`). Point `i` is `g_i site0`, with `g_i` the `i`-th symmetry operation of Oh or D6h
  (in the order of `group.pkl`).
- `site0` (Oh) is `(1, 5, 6)` (not normalized, length `sqrt(62)`); the D6h one, `(2 sqrt(2), sqrt(6), 4 sqrt(3))`, is
  the Oh one with the principal axis rotated to `[111]`. `vector0` (bond direction) is the unit vector of
  `(6, 1, 5)` (rotated in the same way for D6h); `(5, 1, 6)` is not used, because it is the mirror image of `site0` on
  the plane `x + y = z` and made some bond functions vanish.
- The general point is a general position of both Oh and D6h, so the cluster has as many points as the order of the
  group, and the cluster functions form the regular representation.

### 4.2 Scalar cluster basis (`cluster_scalar_basis.pkl`)

1. For every real top harmonics `T` (rank 0..11), the vector of its values at the cluster points.
2. Gram-Schmidt within each (irrep, component) of the top group, from lower rank, then lower multiplicity; linearly
   dependent ones are dropped (whole multiplets; recorded in `"dropped"`).
3. Result: an orthonormal basis of 48 (Oh) or 24 (D6h) vectors, the regular representation.
4. Key: `PGMultipoleType("Q", l, Gamma, n, -1, 0, 0, "q")` with `(l, Gamma, n)` of the top harmonics.

### 4.3 Vector cluster basis (`cluster_vector_basis.pkl`, s = 1)

An intermediate basis for the bond cluster basis.

1. For every real top harmonics `T` and `k = -1, 0, 1`, the vector harmonics `sum_m U_T[m] O_{lm}^{(1)}(k)` with
   the internal basis `[px, py, pz]`, evaluated at the cluster points (`(N, 3)`).
2. Type: for `x = q`, `X = Q` if `s + k` is even, `G` otherwise; for `x = g`, interchanged. A polar one keeps the
   irrep/component of `T`; an axial one takes those of `data_harmonics_mapping[top]["G"]` with the coefficient `c`.
3. Gram-Schmidt within each (irrep, component) in the order rank `l`, then `k = -1, 0, 1`, then multiplicity;
   dependent ones are dropped (found numerically, the kept ones orthonormalized exactly).
4. Result: `3N` = 144 (Oh) or 72 (D6h) vectors for each `x`.
5. Key: `PGMultipoleType(X, l, Gamma, n, -1, 1, k, x)` (`p = -1`; multiplets are distinguished by `k`).

### 4.4 Bond cluster basis (`cluster_bond_basis.pkl`)

1. Bond vectors: `b_i` is point `i` of the cluster with `x, y, z` replaced by `vector0`, i.e. `g_i vector0`.
2. For each vector cluster basis `v_key` (`(N, 3)`), the bond function `f_key(i) = v_key(i) . b_i` (same key);
   `3N` vectors in the `N`-dimensional space.
3. Gram-Schmidt within each (irrep, component) in the order rank, `k = -1, 0, 1`, multiplicity; dependent and
   vanishing ones are dropped. Result: `N` vectors for each `x`.
4. `f_key` transforms as a scalar function on the cluster, with the sign `det R` for `x = g`.
5. Key: `PGMultipoleType(X, l, Gamma, n, p, 0, 0, x)`; `(X, l, Gamma, n)` from the
   source vector key; `p = 1, 2, ...` for multiplets sharing `(X, l, Gamma, n)` (in `k` order), `-1` if unique; `x`
   of the source. `"source"` maps each key to its vector key (with `k`).

---

## 5. Cluster SAMB of each Wyckoff site and bond (`cluster_samb.pkl`)

Code: `util_cluster_samb`, `create_binary/create_cluster_samb.py`. Range: 37 real point groups and 230 space groups;
`"site"`, `"bond_s"` (symmetric bond), `"bond_a"` (anti-symmetric bond). Complex point groups are not included.

1. **General point.** The general-point cluster of a group `H` is a subset of the top cluster: point `i` of `H` is
   the top point listed at position `i` of coset 0 of `data_op_mapping.op_mapping_d6h_oh_cluster[PG]`.
   - The top basis restricted to these points gives the general-point basis of `H`: `cluster_scalar_basis` for
     `"site"` and `"bond_s"`, and `cluster_bond_basis` (`x = q` only, a complete basis by itself) for `"bond_a"`.
     On the general bond, symmetric and anti-symmetric bonds span the same space.
   - The key of the top group is kept; the irrep and component are replaced by those of `H` through the top
     harmonics: key of the top group `= c_top T`, key of `H` `= c T`, factor `c c_top`.
2. **Space groups.** The operations of a space group (without the plus set) have the point parts of its point group
   in the same order. Each operation is mapped through the point group. (`op_mapping_parent_sg` maps to a parent
   space group with changes of axes and is not used.)
3. **Special points.** For a Wyckoff site (bond) with the mapping `M_j` (operations moving the representative to site
   or bond `j`, from `group.pkl`), `f_j = sum_{i in M_j} f_i`. For `"bond_a"`, an operation that reverses the bond
   (negative number in the mapping) contributes `-f_i`.
4. **Orthonormalization.** Gram-Schmidt within each irrep of `H` (component by component; components of different
   multiplets of a real 2D irrep of complex type are not orthogonal), in the order `[Gamma, l, n, p]` for `"site"`
   and `"bond_s"`, and `[Gamma, l, k, Q/G, n, p]` for `"bond_a"` (`k` of the source vector basis).
   Dependent and vanishing ones are dropped (selected numerically, then orthonormalized exactly). Each component
   has norm 1 on the `N` sites (bonds).
5. **Key.** `"site"`, `"bond_s"`: `PGMultipoleType("Q", l, Gamma, n, p, 0, 0, "q")`. `"bond_a"`: keys of
   `cluster_bond_basis` (`X = Q/G`, `x = q`). `n` is the harmonics multiplicity of `H`; `p` that of the top key.
6. **Value.** `samb (dim, N)`. The expression of the harmonics of the point group is in `harmonics.pkl`.

---

## 6. Atomic multipoles

Code: `util_atomic_multipole`, `create_binary/create_atomic_multipole.py`,
`create_binary/create_atomic_multipole_group.py`.

### 6.1 Atomic multipole matrices (`atomic_multipole.pkl`)

1. Bases (`data_basis`): spinful LMs (`atomic_lms_basis`, `(m, spin)`) and JML (`atomic_jm_basis`, `(J, M)`) of
   `L = 0..3`; spinless LM is the up-spin part.
2. Matrix elements `<b1|X_{lm}^{(s)}(k)|b2>` for `X = Q, G, T, M`, `s = 0, 1`, `l <= 7`, by the Wigner-Eckart
   theorem (`multipie.util.util_atomic_multipole`). Only elements allowed by Jz conservation are evaluated
   (`m1 + s1 = m2 + s2 + m`, `M1 = M2 + m`), and blocks forbidden by parity or the triangle rule are skipped.
3. Each bra-ket block `(L1, L2)` of `data_basis.atomic_rank_block` (upper triangle) is kept if nonzero.
4. Key: `SphericalMultipoleType(X, l, s, k, "q")`; order `Q, G, T, M`, `s`, `k`, `l`. Value: `(2l+1, n1, n2)` for
   `m = l..-l`. `"lm"` is the `s = 0`, up-spin part of `"lms"`.

### 6.2 Point-group SAMB (`atomic_multipole_group.pkl`)

1. **Top SAMB.** For each component of a top harmonics `T` (u-matrix `U_T`, `m = l..-l`) and each atomic multipole
   with matching parity (`Q, T` for polar harmonics, `G, M` for axial), `sum_m U_T[m] X_{lm}^{(s)}(k)`, normalized
   to the Frobenius norm 1 (diagonal block) or `1/sqrt(2)` (off-diagonal block). Computed once per top group.
2. The components are orthogonal (the spherical tensors are orthogonal and the `U_T` columns orthonormal), so this
   normalization is the Gram-Schmidt orthonormalization within each (`X`, irrep). In an off-diagonal block, Q-T (and
   G-M) of the same irrep overlap only in the imaginary part (they become orthogonal with the hermitian conjugate
   block).
3. **Point group.** The SAMB of a harmonics `(X, l, Gamma, n)` of the group is `c` times the top SAMB of its top key
   (`data_harmonics_mapping`), with the components stacked.
4. **Bases.** `"jml"` from the JML matrices. `"lgs"` is `U_bra^dagger M U_ket` from the LMs matrices, with `U` the
   orbitals of `data_basis.atomic_lgs_basis` (`atomic_orbital_basis`, basis-function u-matrices, spin expanded).
   `"lg"` is the `s = 0`, up-spin part of `"lgs"` times `sqrt(2)`.
5. **Key.** `PGMultipoleType(X, l, Gamma, n, -1, s, k, "q")`; `n` is the harmonics multiplicity (not renumbered).
   Order `Q, G, T, M`, `s`, `k`, `l`, `Gamma`, `n`.
6. **Value.** `matrix (dim, n1, n2)` for each block (the expression of the harmonics is in `harmonics.pkl`). All 47
   point groups.

---

## 7. Exact arithmetic and simplification

- All stored values are sympy. Saved values go through `util_symmetry_operation.simplify_value`
  (`expand(radsimp(expand(x)))`; rationals unchanged), which rationalizes denominators (`sqrt(3)/3` instead of
  `1/sqrt(3)`).
- Linear dependence is decided numerically (relative tolerance `1e-8`); the kept vectors are then orthonormalized
  exactly.
- sympy re-evaluates some products when a pickle is loaded (`(a+b)/2` → `a/2 + b/2`), so "binary is up to date"
  tests compare after a pickle round trip.
