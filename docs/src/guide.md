# User Guide

A Japanese version is [here](guide_ja.md).

This guide describes the quantities provided by MultiPieGen: their definitions (the first part) and how they are
generated (the second part).

## Definitions

This part defines the quantities provided by MultiPieGen. How to get the actual function forms and $U$
matrices from the data is described in a later section. The harmonics and the SAMB are labeled by the multipole key
`PGMultipoleType` and written with its LaTeX symbol, e.g. $\mathbb{Q}_{1}(T_{1u})$; both are explained in the
last section of this part, [Multipole key](#multipole-key).

### 1. Point-group harmonics

**Spherical harmonics.** The harmonics of rank `l` are built from

$$
O_{lm}(\boldsymbol{r}) = \sqrt{\frac{4\pi}{2l+1}}\, r^l\, Y_{lm}(\theta, \phi), \qquad m = l, l-1, \dots, -l,
$$

homogeneous polynomials of degree `l` (e.g. $O_{10} = z$, $O_{1,\pm1} = \mp(x \pm iy)/\sqrt{2}$). Their real
combinations are the tesseral harmonics ($m > 0$)

$$
C_{lm} = \frac{(-1)^m O_{lm} + O_{l,-m}}{\sqrt{2}}, \qquad
S_{lm} = \frac{(-1)^m O_{lm} - O_{l,-m}}{\sqrt{2}\, i}, \qquad C_{l0} = O_{l0},
$$

e.g. $C_{11} = x$, $S_{11} = y$, $C_{22} = \sqrt{3}(x^2 - y^2)/2$, $S_{22} = \sqrt{3}xy$.

**Point-group harmonics.** A point-group harmonics `(X, l, Gamma, n)` of a point group is a set of rank-`l` functions that
forms a basis of the irrep `Gamma`. Its component `gamma` is

$$
X_{l\Gamma n}^{(\gamma)}(\boldsymbol{r}) = \sum_{m} U_{m\gamma}\, O_{lm}(\boldsymbol{r}),
$$

with $U$ (`(2l+1) x dim`, rows `m = l, ..., -l`), whose columns are orthonormal.

*Example (Oh).* The dipole $\mathbb{Q}_{1}(T_{1u}) = (x, y, z)$ has

$$
U = \begin{pmatrix} -\frac{1}{\sqrt{2}} & \frac{i}{\sqrt{2}} & 0 \\ 0 & 0 & 1 \\ \frac{1}{\sqrt{2}} & \frac{i}{\sqrt{2}} & 0 \end{pmatrix}
\quad (m = 1, 0, -1),
$$

e.g. $x = (-O_{11} + O_{1,-1})/\sqrt{2}$, and the quadrupole $\mathbb{Q}_{2}(E_{g}) = \left(z^2 - \frac{x^2+y^2}{2},
\frac{\sqrt{3}}{2}(x^2 - y^2)\right)$.

- **Top groups.** The harmonics are defined once for the two top groups, Oh (cubic) and D6h (trigonal and
  hexagonal), up to rank 11. Every one of the 47 point groups (32 crystallographic, 5 other settings `-1`, 10
  complex `-c`) takes its harmonics from them: each component is a component of a top harmonics times a sign,
  $X_{\rm PG} = c\, X_{\rm top}$ with $c = \pm 1$.
- The representation matrices of the irreps follow from $U$, and all harmonics of an irrep give the same matrices.
- The irreps follow the order of the character table. Complex irreps (in the `-c` groups) end in `a`/`b` and have a
  single component.
- *Multiplicity.* When an irrep appears more than once at the same rank, `n = 1, 2, ...` distinguishes the
  harmonics, e.g. C4h: $\mathbb{Q}_{2}(B_{g},1) = \frac{\sqrt{3}}{2}(x^2 - y^2)$ and $\mathbb{Q}_{2}(B_{g},2) = \sqrt{3}xy$.
- *Complex group.* E.g. C4-c: $\mathbb{Q}_{1}(E^{(a)}) = (x + iy)/\sqrt{2}$ and $\mathbb{Q}_{1}(E^{(b)}) = (x - iy)/\sqrt{2}$.
- Key: `PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q")` (see [Multipole key](#multipole-key)).

**Polar (`Q`) and axial (`G`) harmonics.** An axial harmonics transforms with an extra factor `det R` under an
improper operation `R` (a mirror, the inversion, or a rotoinversion). Both come from the same top harmonics, but in a
point group a polar and an axial harmonics of the same function can belong to different irreps, and in the same irrep
their components can be in another order or have other signs (the coefficient `c` and the component assigned differ
between `Q` and `G`). So, for the same `(l, Gamma, n)`, the functions of `Q` and `G` differ in general:

- In the top groups Oh and D6h, the axial harmonics is the polar one with `g` and `u` interchanged in the irrep: the
  same function, components, and order ($c = 1$). E.g. Oh: $\mathbb{Q}_{1}(T_{1u}) = \mathbb{G}_{1}(T_{1g}) = (x, y, z)$.
- In a group of proper rotations only (C1, C2, D2, C4, D4, C3, D3, D3-1, C6, D6, T, O), there is no improper
  operation, and `Q` and `G` are identical. E.g. C4: $\mathbb{Q}_{1}(E) = \mathbb{G}_{1}(E) = (x, y)$.
- In the other groups, they differ. Examples:

| group | polar | axial | difference |
|---|---|---|---|
| C4v | $\mathbb{Q}_{1}(A_{1}) = z$ | $\mathbb{G}_{1}(A_{2}) = z$ | same function, another irrep |
| C4v | $\mathbb{Q}_{1}(E) = (x, y)$ | $\mathbb{G}_{1}(E) = (-y, x)$ | components interchanged, with a sign |
| C4h | $\mathbb{Q}_{1}(E_{u}) = (x, y)$ | $\mathbb{G}_{1}(E_{g}) = (y, x)$ | components interchanged |
| D2d | $\mathbb{Q}_{1}(E) = (x, y)$ | $\mathbb{G}_{1}(E) = (x, -y)$ | sign of a component |
| C4v | $\mathbb{Q}_{2}(B_{1}) = \frac{\sqrt{3}}{2}(x^2 - y^2)$ | $\mathbb{G}_{2}(B_{1}) = \sqrt{3}xy$ | another function for the same irrep |
| D2d | $\mathbb{Q}_{2}(A_{1}) = z^2 - \frac{x^2+y^2}{2}$ | $\mathbb{G}_{2}(A_{1}) = \frac{\sqrt{3}}{2}(x^2 - y^2)$ | another function for the same irrep |

### 2. PG multipolar harmonics

A PG multipolar harmonics combines a point-group harmonics with an internal degree of freedom.

**Internal degree of freedom.** First, the internal degree of freedom is fixed by its rank `s` (`s = 1`: a vector;
`s = 2, 3`: higher-rank tensors) and its type `x`: polar (`x = q`, e.g. a displacement or an electric polarization)
or axial (`x = g`, e.g. a spin or a rotation). Its basis $e_\nu$ is the set of rank-`s` orbitals of that type, of the
cubic family (Oh) or of the hexagonal family (D6h), written with the symbols $Q_\nu$ (`x = q`) or $G_\nu$ (`x = g`),
e.g. $(Q_x, Q_y, Q_z)$ or $(G_x, G_y, G_z)$ for `s = 1`.

**Multipolar spherical harmonics.** The multipolar spherical harmonics of total rank `l` couple the orbital part of
rank `l + k` with the internal basis:

$$
O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} \sqrt{2l+1} \sum_{\nu}
\begin{pmatrix} l+k & l & s \\ m-\nu & -m & \nu \end{pmatrix} O_{l+k, m-\nu}\, e_\nu,
\qquad k = -s, \dots, s,
$$

**Type of the multipole.** For a fixed internal type `x`, the type `X` of the multipole is determined by the parity
of `s + k`:

| internal type | `s + k` even | `s + k` odd |
|---|---|---|
| `x = q` (polar) | `Q` | `G` |
| `x = g` (axial) | `G` | `Q` |

That is, `X` has the same polar or axial character as the internal basis if `s + k` is even, and the opposite one if
`s + k` is odd. The PG multipolar harmonics of a point group is $\sum_m U_{m\gamma} O_{lm}^{(s)}(k)$ with $U$ of the
point-group harmonics `(X, l, Gamma, n)`.

- Range: `s = 1, 2, 3`, `l <= 3`, all 47 point groups. Key: `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` (see
  [Multipole key](#multipole-key)).
- The value is the coefficient of each internal basis function (polynomials in x, y, z) and the resulting expression
  in terms of the internal basis symbols.

*Examples (Oh, `s = 1`).* With the polar vector $\boldsymbol{Q} = (Q_x, Q_y, Q_z)$ (`x = q`) or the axial vector
$\boldsymbol{G} = (G_x, G_y, G_z)$ (`x = g`) as the internal basis:

| symbol | expression | meaning |
|---|---|---|
| $\vec{\mathbb{Q}}_{0}^{(1,1)}[q](A_{1g})$ | $(xQ_x + yQ_y + zQ_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{Q}$, electric monopole |
| $\vec{\mathbb{Q}}_{1}^{(1,-1)}[q](T_{1u})$ | $(Q_x, Q_y, Q_z)$ | the polar vector itself |
| $\vec{\mathbb{G}}_{1}^{(1,0)}[q](T_{1g})$ | $(zQ_y - yQ_z, xQ_z - zQ_x, yQ_x - xQ_y)/\sqrt{2}$ | $\boldsymbol{Q} \times \boldsymbol{r}$, electric toroidal dipole |
| $\vec{\mathbb{G}}_{0}^{(1,1)}[g](A_{1u})$ | $(xG_x + yG_y + zG_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{G}$, electric toroidal monopole |

### 3. Atomic orbitals

The atomic orbitals of rank `L = 0, 1, 2, 3` (s, p, d, f) are the basis of the atomic SAMB. Three bases are used.

**`lg` (real orbitals, spinless).** The tesseral harmonics, in the order $C_{LL}, S_{LL}, C_{L,L-1}, S_{L,L-1},
\dots, C_{L0}$:

| `L` | orbitals | functions |
|---|---|---|
| 0 | `s` | $1$ |
| 1 | `px`, `py`, `pz` | $x$, $y$, $z$ |
| 2 | `dv`, `dxy`, `dxz`, `dyz`, `du` | $\frac{\sqrt{3}}{2}(x^2-y^2)$, $\sqrt{3}xy$, $\sqrt{3}xz$, $\sqrt{3}yz$, $z^2 - \frac{x^2+y^2}{2}$ |
| 3 | `fc`, `fs`, `fbz`, `f3`, `fx`, `fy`, `faz` | $\frac{\sqrt{10}}{4}x(x^2-3y^2)$, $\frac{\sqrt{10}}{4}y(3x^2-y^2)$, $\frac{\sqrt{15}}{2}z(x^2-y^2)$, $\sqrt{15}xyz$, $\frac{\sqrt{6}}{4}x(4z^2-x^2-y^2)$, $\frac{\sqrt{6}}{4}y(4z^2-x^2-y^2)$, $\frac{1}{2}z(2z^2-3x^2-3y^2)$ |

**`lgs` (real orbitals with spin).** `lg` times the spin `u` (up) and `d` (down), with the spin fastest:
`(px, u), (px, d), (py, u), ...`.

**`jm` (total angular momentum).** $|J, M\rangle$ with $J = L \mp 1/2$, in the order `J` ascending and `M`
descending, e.g. for `L = 1`: $|1/2, 1/2\rangle, |1/2, -1/2\rangle, |3/2, 3/2\rangle, \dots, |3/2, -3/2\rangle$
(`L = 0` has $J = 1/2$ only). In the data this basis is named `jml`.

**Blocks.** A matrix in these bases is given for each block of bra rank `L1` and ket rank `L2`, for
`(L1, L2) = (0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)` (`L1 <= L2`).
E.g. the block `(1, 2)` (p-d) is `3 x 5` in `lg`, and `6 x 10` in `lgs` and `jml`. In `lgs`, the p orbitals are
`(px, u), (px, d), (py, u), (py, d), (pz, u), (pz, d)`.

### 4. Atomic SAMB

The atomic SAMB (symmetry-adapted multipole basis) are the matrices of the atomic multipole operators in the orbital
basis, adapted to the point group.

1. **Atomic multipoles.** The spherical multipole operators $X_{lm}^{(s)}(k)$ of `X = Q, G, T, M` (`Q` and `G`
   even, `T` and `M` odd under time reversal), with `s = 0` (orbital) and `s = 1` (spin-dependent), `k = -s, ..., s`.
   Their matrix elements between orbitals of ranks `L1` and `L2` follow from the Wigner-Eckart theorem.
   - For `s = 1`, the internal degree of freedom is fixed to the spin, an axial vector: the case `x = g` of the
     [PG multipolar harmonics](#pg-multipolar-harmonics). $X_{lm}^{(1)}(k)$ couples an orbital multipole of rank `l + k`
     with the spin. In the PG multipolar harmonics the orbital part $O_{l+k,m}$ is always a polar function, while here
     the orbital multipole can be polar (`Q`, `T`) or axial (`G`, `M`, e.g. the orbital angular momentum). So the rule
     there for `x = g` is taken relative to the orbital multipole: `X` has the polar or axial character opposite to
     that of the orbital multipole if `s + k` is even, and the same one if `s + k` is odd. In addition, the spin is
     odd under time reversal, so `X` has the opposite time-reversal parity to the orbital multipole.
   - E.g. the spin magnetic dipole $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ is the orbital monopole `Q` (rank 0, polar,
     even) with the spin: `s + k = 0` even gives axial and the spin gives odd, i.e. `M`. The spin-orbit coupling
     $\boldsymbol{l} \cdot \boldsymbol{s}$ is $\sqrt{3}\, \mathbb{Q}_{0}^{(1,1;a)}(A_{1g})$ for the p orbitals: the orbital
     angular momentum `M` (rank 1, axial, odd) with the spin, `s + k = 2` even gives polar and the spin gives even,
     i.e. `Q`.
2. **Symmetry adaptation.** With $U$ of the point-group harmonics of the type of the whole multipole (polar for
   `Q, T`, axial for `G, M`), $X_{l\Gamma n}^{(s)}(k)[\gamma] = \sum_m U_{m\gamma} X_{lm}^{(s)}(k)$.
3. **Normalization.** Each component has the Frobenius norm 1 in a diagonal block (`L1 = L2`), and `1/sqrt(2)` in
   an off-diagonal block (norm 1 together with its hermitian conjugate block `(L2, L1)`).

- Bases: `lg` (`s = 0` only), `lgs`, and `jml`.
- Key: `PGMultipoleType(X, l, Gamma, n, -1, s, k, "q")` (see [Multipole key](#multipole-key)). Value: `(dim, n1, n2)`, the matrix of each component for
  the `n1` orbitals of rank `L1` and the `n2` orbitals of rank `L2`.
- Order: `Q, G, T, M`, then `s`, `k`, `l`, `Gamma`, `n`. All 47 point groups.

*Examples (Oh).* For the p orbitals (`lg`, block `(1, 1)`, basis `px, py, pz`):

- $\mathbb{Q}_{0}^{(a)}(A_{1g}) = \frac{1}{\sqrt{3}}\,\mathrm{diag}(1, 1, 1)$, the charge;
- $\mathbb{Q}_{2}^{(a)}(E_{g}) = \left(\frac{1}{\sqrt{6}}\,\mathrm{diag}(-1, -1, 2),\ \frac{1}{\sqrt{2}}\,\mathrm{diag}(1, -1, 0)\right)$, the
  orbital quadrupole;
- $\mathbb{M}_{1}^{(a)}(T_{1g})$: $(l_x, l_y, l_z)/\sqrt{2}$, the orbital angular momentum, e.g. the `x` component is
  $l_x/\sqrt{2}$, with the matrix element $\langle p_y | l_x | p_z \rangle / \sqrt{2} = -i/\sqrt{2}$.

For the s orbital with spin (`lgs`, block `(0, 0)`, basis `(s, u), (s, d)`): $\mathbb{Q}_{0}^{(a)}(A_{1g})$ is
$\sigma_0/\sqrt{2}$, and $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ is $(\sigma_x, \sigma_y, \sigma_z)/\sqrt{2}$, the
spin.

### 5. Cluster SAMB

The cluster SAMB of a Wyckoff position are functions on its `N` equivalent sites (or bonds), which form bases of the
irreps of the group. They are given for the 37 real point groups and the 230 space groups.

- **`site`.** Functions on the `N` sites of a site Wyckoff position (e.g. `4b`). They are real, have `X = Q`, and
  form a complete orthonormal set (`N` functions).
- **`bond`.** Functions on the `N` bonds of a bond Wyckoff position. The label `2b@1o` is a bond Wyckoff position
  `2b` (2 bonds) whose centers are at the site Wyckoff position `1o`. Two kinds are put together:
  - symmetric bond (the same value for the two directions of a bond, a real hopping): `X = Q`, real values, a
    complete orthonormal set (`N` functions);
  - anti-symmetric bond (opposite signs for the two directions, an imaginary hopping): `i` times the anti-symmetric
    function, with `X = T` (polar) or `M` (axial), pure imaginary values.

  Q and T/M are orthonormal separately; they are orthogonal to each other as hermitian bond matrices (different
  time-reversal parity).
- Key: `PGMultipoleType(X, l, Gamma, n, p, 0, 0, "q")` (see [Multipole key](#multipole-key)). `(X, l, Gamma, n)` names the point-group harmonics that the
  function corresponds to (`T` and `M` with those of `Q` and `G`), and `p` distinguishes functions sharing it.
- Value: `(dim, N)`, the values on the `N` sites (bonds) for each component; each component has norm 1.
- Order: `Gamma` (character table), and in `bond`, `Q` before `T/M` in each `Gamma`.

*Examples (D4h).*

- Site `2a`: two sites at $(0, 0, \pm z)$. $\mathbb{Q}_{0}^{(s)}(A_{1g}) = (1, 1)/\sqrt{2}$ (uniform) and $\mathbb{Q}_{1}^{(s)}(A_{2u}) = (1, -1)/\sqrt{2}$
  (a dipole along `z`).
- Bond `2b@1o`: two bonds through the origin, along `x` and along `y`.
  $\mathbb{Q}_{0}^{(b)}(A_{1g}) = (1, 1)/\sqrt{2}$, $\mathbb{Q}_{2}^{(b)}(B_{1g}) = (-1, 1)/\sqrt{2}$ (real hoppings), and
  $\mathbb{T}_{1}^{(b)}(E_{u},a) = \left((i, 0), (0, i)\right)$ (imaginary hoppings): the `x` component lives on the bond along `x`, the `y`
  component on the bond along `y`.

### 6. Response tensors

A response tensor of rank 0 to 4 relates two physical quantities, e.g. a current and an electric field. For each of
the 157 magnetic point groups, MultiPieGen gives the tensors allowed by symmetry.

- Type `X` by the parity under the spatial inversion and the time reversal: `Q` (polar, time-reversal even), `G`
  (axial, even), `T` (polar, odd), `M` (axial, odd).
- `opt` specifies the part of the tensor by the symmetry of its indices: `s` symmetric and `a` anti-symmetric in a
  pair of indices, as in the table below (for a polar, time-reversal-even tensor `X = Q`; the other types have the same
  relations). The indices before `;` form the pair (or triple) of the rows of the matrix, and those after it the
  columns. The index sets are `c = [x, y, z]`, `v = [xx, yy, zz, yz, zx, xy]` (symmetric pairs), `a = [yz, zx, xy]`
  (anti-symmetric pairs), `u = [xx, yy, zz, yz, zx, xy, zy, xz, yx]` (all pairs), and
  `t = [xxx, yyy, zzz, yyz, zzx, xxy, yzz, zxx, xyy, xyz]` (symmetric triples). "independent" is the number of
  independent parameters of the part.

| rank | `opt` | relation among the components | independent | matrix (rows x columns) |
|---|---|---|---|---|
| 0 | (none) | - | 1 | 1 x 1 |
| 1 | (none) | - | 3 | 1 x 3 (`c`) |
| 2 | `s` | $X_{ab} = X_{ba}$ | 6 | 3 x 3 (`c` x `c`) |
| 2 | `a` | $X_{ab} = -X_{ba}$ | 3 | 3 x 3 (`c` x `c`) |
| 3 | `ss` | $X_{abc}$ fully symmetric in $a, b, c$ | 10 | 6 x 3 (`v` x `c`) |
| 3 | `sa` | $X_{ab;c} = X_{ba;c} = -(X_{ac;b} + X_{bc;a})$ | 8 | 6 x 3 (`v` x `c`) |
| 3 | `a` | $X_{ab;c} = -X_{ba;c}$ | 9 | 3 x 3 (`a` x `c`) |
| 4 | `sss` | $X_{ab;cd} = X_{ba;cd} = X_{ab;dc} = X_{cd;ab}$ | 21 | 6 x 6 (`v` x `v`) |
| 4 | `ssa` | $X_{ab;cd} = X_{ba;cd} = X_{ab;dc} = -X_{cd;ab}$ | 15 | 6 x 6 (`v` x `v`) |
| 4 | `aas` | $X_{ab;cd} = -X_{ba;cd} = -X_{ab;dc} = X_{cd;ab}$ | 6 | 3 x 3 (`a` x `a`) |
| 4 | `aaa` | $X_{ab;cd} = -X_{ba;cd} = -X_{ab;dc} = -X_{cd;ab}$ | 3 | 3 x 3 (`a` x `a`) |
| 4 | `sa` | $X_{ab;cd} = X_{ba;cd} = -X_{ab;dc}$ | 18 | 6 x 3 (`v` x `a`) |
| 4 | `as` | $X_{ab;cd} = -X_{ba;cd} = X_{ab;dc}$ | 18 | 3 x 6 (`a` x `v`) |
| 4 | `s` | $X_{ab;c;d} = X_{ba;c;d}$ (`sss` + `ssa` + `sa`) | 54 | 6 x 9 (`v` x `u`) |
| 4 | `a` | $X_{ab;c;d} = -X_{ba;c;d}$ (`aas` + `aaa` + `as`) | 27 | 3 x 9 (`a` x `u`) |
| 4 | `t` | $X_{abc;d}$ fully symmetric in $a, b, c$ | 30 | 10 x 3 (`t` x `c`) |

- Each tensor component is decomposed into multipoles (irreducible cartesian parts). Only the active multipoles of
  the group, those invariant under all its operations including the time reversal, can be nonzero. This gives the
  nonvanishing components and the relations among them.
- Value for `(X, rank, opt)`: the matrix of the components (symbols such as $Q_{xyz}$, `0` for vanishing ones;
  rows and columns as in the table above), and the expression of each component
  in terms of the active multipoles.

*Examples.*

- 4/mmm (`MPG:15.1.53`), `(Q, 2, s)` (e.g. the electric conductivity): $\mathrm{diag}(Q_{xx}, Q_{xx}, Q_{zz})$, with
  $Q_{xx} = \mathbb{Q}_0^{(1)} - \mathbb{Q}_u^{(1)}$ and $Q_{zz} = \mathbb{Q}_0^{(1)} + 2\mathbb{Q}_u^{(1)}$.
- 4/mm'm' (`MPG:15.6.58`, a ferromagnet along `z`), `(T, 2, a)` (e.g. the anomalous Hall effect): the only nonzero
  components are $T_{xy} = -T_{yx} = \mathbb{M}_z^{(1)}$, given by the active magnetic dipole $M_z$.

### 7. Cosets in the parents

A group `H` (a real point group or a space group) with a parent `P` (`data_parent_pg/sg`) is split into the left
cosets `g H` of `P`. A representative `g` maps the structure of `H` to another one, which has the group
`g H g^-1`; on a structure invariant under `P` (the parent structure), `g` is a symmetry. MultiPieGen gives, for every
group with a parent (not the complex point groups):

- **Parents and cosets.** The chain of parents up to the top group (Oh, D6h, SG 221-230, 191-194), the direct parent
  first, and the left cosets in each parent (the operations of the parent in each coset; for SG with the translation in
  the coordinates of `H` modulo its lattice). `H` itself and the cosets of translations only are left out.
- **First coset.** The first coset other than `H` in the direct parent, its representative `g` (the first operation),
  the order `n` of `g` (the smallest `n` with `g^n` in `H`), and the kind of `g`: a rotation, or the inversion (with a
  translation for SG). The inversion appears only for a chiral `H` (77 groups: 12 point groups and 65 space groups).
  - Rotation: the image structure is the structure of `H` written in the frame moved by `g`, and `g Z_j = +Z'_j` for
    every SAMB (the SAMB `Z'_j` of the image with the same label).
  - Inversion: the image group `H' = I H I^-1` (an enantiomorph, or `H` itself), the map `M` from the coordinates of
    `H` to those of `H'`, the image of each Wyckoff site and bond (Wyckoff position, order of its sites or bonds, and
    the reversal of a bond), and the image of the cluster SAMB: `I Z_j = eps_j Z'_j` label by label (`eps = +-1`), or
    `I Z_a = sum_b U_ab Z'_b` within an irrep for some Wyckoff positions mapped onto themselves (38 Wyckoff positions,
    e.g. 48f of SG 210).
- **Parent structure.** For each Wyckoff position with free parameters, the values of the parameters at which the
  orbit of `H` is invariant under the direct parent with no coinciding sites (the largest such sets), the Wyckoff
  position of the parent occupied there, and the cosets that keep it.

*Example (Te, P3_121 = No.152).* The direct parent is P-3m1 (No.164), the first coset has `g = -1` (`n = 2`), and its
image is P3_221 (No.154) with `x -> -x`. Every site and bond SAMB of `3a` is mapped with `eps = +1`. The parent
structure of `3a` is `x = 0` (the Wyckoff position `1a` of P-3m1, whose cell is a third of that of `H`).

(multipole-key)=
### 8. Multipole key (`PGMultipoleType`)

This section explains the key of the multipoles, which is needed to use the harmonics and the SAMB.

Every multipole of a point group (point-group harmonics, PG multipolar harmonics, atomic SAMB, cluster SAMB) is labeled by
the named tuple `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` of the `multipie` package:

| field | meaning | values |
|---|---|---|
| `X` | type, by the parity under the spatial inversion ($\mathcal{P}$) and the time reversal ($\mathcal{T}$) | `Q` electric ($\mathcal{P} = (-1)^l$, $\mathcal{T}$ even), `G` electric toroidal ($(-1)^{l+1}$, even), `T` magnetic toroidal ($(-1)^l$, odd), `M` magnetic ($(-1)^{l+1}$, odd) |
| `l` | rank | 0, 1, 2, ... |
| `Gamma` | irrep of the point group | name in the character table, e.g. `A1g`, `Eu`, `T2g` |
| `n` | multiplicity of the point-group harmonics `(X, l, Gamma)`: which of the harmonics with the same `(X, l, Gamma)` | 1, 2, ...; `-1` if unique |
| `p` | extra multiplicity, numbering entries that share all the other fields | 1, 2, ...; `-1` if unique |
| `s` | internal rank (rank of the internal degree of freedom) | 0 (scalar), 1, 2, 3 |
| `k` | internal component: the orbital part has rank `l + k` | `-s, ..., s` |
| `x` | internal type: parity of the internal basis | `q` (polar), `g` (axial) |

- `(X, l, Gamma, n)` always identifies a point-group harmonics (its xyz polynomial and $U$); `T` and `M` use the
  harmonics of `Q` and `G`. `n` is never renumbered, so different data with the same `(X, l, Gamma, n)` refer to the
  same harmonics.
- `p` appears only where several entries share `(X, l, Gamma, n, s, k, x)`, e.g. two cluster SAMB of the same
  harmonics on a Wyckoff position.
- Scalar quantities (point-group harmonics, cluster SAMB) have `s = 0, k = 0, x = "q"`.
- A key names a whole multiplet. Its value holds all components together, in the order of the components of the
  irrep (rows of the array).

#### LaTeX symbol

In the documentation and the PDFs, a multipole is written with the LaTeX symbol made from its key
(`multipie.util.util_tag.TagMultipole.latex`, with the tag below):

$$
\mathbb{X}_{l,\gamma}^{(s,k;\,\mathrm{tag})}[x](\Gamma, n\,p)
$$

| part | from the key | written |
|---|---|---|
| $\mathbb{X}$ | `X` | $\mathbb{Q}, \mathbb{G}, \mathbb{T}, \mathbb{M}$ |
| $l$ | `l` | always |
| $\gamma$ | the component (1, 2, ...) | only for one component of a multiplet; omitted for the whole multiplet and for 1D irreps |
| $(s, k)$ | `s`, `k` | only if `s > 0` |
| tag | the kind of data | none for the harmonics; `a` for the atomic SAMB; for the cluster SAMB, `s` for `site` and `b` for `bond` |
| $[x]$, $\vec{\ }$ | `x`, `s` | for the PG multipolar harmonics (`s > 0`): an arrow on $\mathbb{X}$ and the internal type $[q]$ or $[g]$ |
| $\Gamma$ | `Gamma` | always, e.g. $T_{1u}$; a complex irrep `Ea` is written $E^{(a)}$ |
| $n$ | `n` | only if `n != -1` |
| $p$ | `p` | only if `p != -1`, as a letter `1 -> a, 2 -> b, ...` after $n$ |

So the superscript and the last parentheses tell the values of `s, k` and `n, p`; the fields not shown have their
default values (`n = p = -1`, `s = k = 0`, `x = q`).

Examples:

| key | symbol | meaning |
|---|---|---|
| `("Q", 1, "T1u", -1, -1, 0, 0, "q")` | $\mathbb{Q}_{1}(T_{1u})$ | electric dipole of Oh, components $(x, y, z)$; its `x` component is $\mathbb{Q}_{1,1}(T_{1u})$ |
| `("Q", 2, "Bg", 2, -1, 0, 0, "q")` | $\mathbb{Q}_{2}(B_{g},2)$ | the second ($n = 2$) quadrupole of $B_g$ in C4h, $\sqrt{3}xy$ |
| `("G", 0, "A1u", -1, -1, 1, 1, "g")` | $\vec{\mathbb{G}}_{0}^{(1,1)}[g](A_{1u})$ | PG multipolar harmonics: electric toroidal monopole from the rank-1 orbital part and an axial vector |
| `("M", 1, "T1g", -1, -1, 1, -1, "q")` | $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ | atomic SAMB: spin magnetic dipole |
| `("T", 1, "Eu", -1, 1, 0, 0, "q")` | $\mathbb{T}_{1}^{(b)}(E_{u},a)$ | bond cluster SAMB: an anti-symmetric bond function of D4h with the symmetry of $(x, y)$ ($p = 1$) |
| `("Q", 2, "E", 1, 2, 0, 0, "q")` | $\mathbb{Q}_{2}(E,1b)$ | $n = 1$ and $p = 2$ |

## Generation rules: algorithms, keys, and orbital basis

This part collects, in one place, how MultiPieGen builds its multipole-keyed data, and the rules for keys and orbital
bases. It describes the algorithms of this package only. Each section names the code that does the work. The changes
from MultiPie Ver.2 are collected separately in "Changes from MultiPie Ver.2" of the Details chapter.

Contents:

1. Common conventions (keys, irreps, top groups, orbital names)
2. Harmonics of the point groups
3. Multipolar harmonics (`multipolar_harmonics.pkl`)
4. Cluster bases of the top groups (`cluster_scalar_basis.pkl`, `cluster_vector_basis.pkl`, `cluster_bond_basis.pkl`)
5. Cluster SAMB of each Wyckoff site and bond (`cluster_samb.pkl`)
6. Atomic multipoles (`atomic_multipole.pkl`) and their point-group SAMB (`atomic_multipole_group.pkl`)
7. Exact arithmetic and simplification
8. Binary data for MultiPie (`binary_data`)

---

### 1. Common conventions

#### 1.1 Multipole key

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

#### 1.2 Irreps

- The order of irreps in `data_character_table` is authoritative for sorting: `g` before `u`, and the `'` block
  before the `''` block (C3h, D3h, D3h-1, C3h-c are in this order).
- Complex irreps (the `-c` point groups) end in `a`/`b`, e.g. `Ega`, `T1ua`. Real 2D irreps of complex type (e.g. `E`
  of C3, C4, C6, S4, C3i, C4h, C6h, T, Th) are used in the real point groups.

#### 1.3 Top groups

- All harmonics are defined once for the two top groups, Oh and D6h, and mapped to every one of the 47 point groups
  (32 crystallographic, 5 other settings `-1`, 10 complex `-c`).
- Cubic point groups use Oh. Trigonal and hexagonal point groups use D6h. Lower-symmetry groups use the top group
  that their parent chain reaches (`data_harmonics_mapping[pg]["top"]`).

#### 1.4 Orbital names

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

### 2. Harmonics of the point groups

Data: `data_top_harmonics` (top harmonics), `data_compatibility_relation`.
Code: `util_top_harmonics` (expressions and u-matrices of the top harmonics), `util_harmonics_mapping`
(`data_harmonics_mapping`, the mapping to each point group), `util_harmonics_multipole.pg_harmonics`.

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

### 3. Multipolar harmonics (`multipolar_harmonics.pkl`)

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

### 4. Cluster bases of the top groups

Code: `util_cluster_basis`, `util_cluster_vector_basis`, `util_bond_cluster_basis` and the matching
`create_binary` modules. Data: `data_root_cluster`, `info["root_cluster"]`.

#### 4.1 Clusters

- The Oh cluster has 48 points and the D6h cluster 24 points, written with `x, y, z` placeholders in cartesian
  coordinates (`data_root_cluster`). Point `i` is `g_i site0`, with `g_i` the `i`-th symmetry operation of Oh or D6h
  (in the order of `group.pkl`).
- `site0` (Oh) is `(1, 5, 6)` (not normalized, length `sqrt(62)`); the D6h one, `(2 sqrt(2), sqrt(6), 4 sqrt(3))`, is
  the Oh one with the principal axis rotated to `[111]`. `vector0` (bond direction) is the unit vector of
  `(6, 1, 5)` (rotated in the same way for D6h); `(5, 1, 6)` is not used, because it is the mirror image of `site0` on
  the plane `x + y = z` and made some bond functions vanish.
- The general point is a general position of both Oh and D6h, so the cluster has as many points as the order of the
  group, and the cluster functions form the regular representation.

#### 4.2 Scalar cluster basis (`cluster_scalar_basis.pkl`)

1. For every real top harmonics `T` (rank 0..11), the vector of its values at the cluster points.
2. Gram-Schmidt within each (irrep, component) of the top group, from lower rank, then lower multiplicity; linearly
   dependent ones are dropped (whole multiplets; recorded in `"dropped"`).
3. Result: an orthonormal basis of 48 (Oh) or 24 (D6h) vectors, the regular representation.
4. Key: `PGMultipoleType("Q", l, Gamma, n, -1, 0, 0, "q")` with `(l, Gamma, n)` of the top harmonics.

#### 4.3 Vector cluster basis (`cluster_vector_basis.pkl`, s = 1)

An intermediate basis for the bond cluster basis.

1. For every real top harmonics `T` and `k = -1, 0, 1`, the vector harmonics `sum_m U_T[m] O_{lm}^{(1)}(k)` with
   the internal basis `[px, py, pz]`, evaluated at the cluster points (`(N, 3)`).
2. Type: for `x = q`, `X = Q` if `s + k` is even, `G` otherwise; for `x = g`, interchanged. A polar one keeps the
   irrep/component of `T`; an axial one takes those of `data_harmonics_mapping[top]["G"]` with the coefficient `c`.
3. Gram-Schmidt within each (irrep, component) in the order rank `l`, then `k = -1, 0, 1`, then multiplicity;
   dependent ones are dropped (found numerically, the kept ones orthonormalized exactly).
4. Result: `3N` = 144 (Oh) or 72 (D6h) vectors for each `x`.
5. Key: `PGMultipoleType(X, l, Gamma, n, -1, 1, k, x)` (`p = -1`; multiplets are distinguished by `k`).

#### 4.4 Bond cluster basis (`cluster_bond_basis.pkl`)

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

### 5. Cluster SAMB of each Wyckoff site and bond (`cluster_samb.pkl`)

Code: `util_cluster_samb`, `create_binary/create_cluster_samb.py`. Range: 37 real point groups and 230 space groups;
`"site"` and `"bond"`. The bond SAMB is made from the symmetric bond SAMB (`"bond_s"` below) and the
anti-symmetric bond SAMB (`"bond_a"` below), see 7. Complex point groups are not included.

1. **General point.** The general-point cluster of a group `H` is a subset of the top cluster: point `i` of `H` is
   the top point listed at position `i` of coset 0 of `util_op_mapping.op_mapping_d6h_oh_cluster[PG]`.
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
7. **Bond SAMB.** `"bond" = "bond_s" + i "bond_a"`: the anti-symmetric bond SAMB is multiplied by `i` (imaginary
   hopping) and its type is changed from `Q` to `T` and from `G` to `M` (time-reversal odd). The keys of `"bond_s"`
   (`Q`) and of `"bond_a"` (`T/M`) never coincide, so they are put together. SAMB of different time-reversal parity
   are orthogonal (as hermitian bond matrices: real symmetric and imaginary anti-symmetric), so `Q` and `T/M` are
   orthonormalized separately (4.). Order: `Gamma`, and in each `Gamma`, `Q` (`[l, n, p]`) before `T/M`
   (`[l, k, T/M, n, p]`). Values are real for `Q` and pure imaginary for `T, M`; the harmonics of `T` (`M`) are those
   of `Q` (`G`) (`util_harmonics.harmonics_type`).
8. **Reference sites for the inversion coset.** The SAMB of a special point depends on the representative (row 0)
   from which it is averaged. For a group whose first coset is the inversion (see 9.), the image of the
   representative of a Wyckoff position `W` need not be the representative of its image `W'`, and the SAMB of an irrep
   can then be mixed under the inversion. For a pair `W != W'` (or `H != H'`) whose SAMB are not mapped label by label
   already, the side with the larger (group, Wyckoff order) is averaged from the image of the representative of the
   other side (`util_cluster_samb_image.cluster_samb_reference`, `create_cluster_samb(..., reference=...)`). This
   applies to 20 Wyckoff positions of SG 98, 208, 211, 213, and 214; a Wyckoff position mapped onto itself keeps its
   data (its matrix `U` is given in `coset.pkl`).

---

### 6. Atomic multipoles

Code: `util_atomic_multipole`, `create_binary/create_atomic_multipole.py`,
`create_binary/create_atomic_multipole_group.py`.

#### 6.1 Atomic multipole matrices (`atomic_multipole.pkl`)

1. Bases (`data_basis`): spinful LMs (`atomic_lms_basis`, `(m, spin)`) and JML (`atomic_jm_basis`, `(J, M)`) of
   `L = 0..3`; spinless LM is the up-spin part.
2. Matrix elements `<b1|X_{lm}^{(s)}(k)|b2>` for `X = Q, G, T, M`, `s = 0, 1`, `l <= 7`, by the Wigner-Eckart
   theorem (`multipie.util.util_atomic_multipole`). Only elements allowed by Jz conservation are evaluated
   (`m1 + s1 = m2 + s2 + m`, `M1 = M2 + m`), and blocks forbidden by parity or the triangle rule are skipped.
3. Each bra-ket block `(L1, L2)` of `data_basis.atomic_rank_block` (upper triangle) is kept if nonzero.
4. Key: `SphericalMultipoleType(X, l, s, k, "q")`; order `Q, G, T, M`, `s`, `k`, `l`. Value: `(2l+1, n1, n2)` for
   `m = l..-l`. `"lm"` is the `s = 0`, up-spin part of `"lms"`.

#### 6.2 Point-group SAMB (`atomic_multipole_group.pkl`)

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

### 7. Exact arithmetic and simplification

- All stored values are sympy. Saved values go through `util_symmetry_operation.simplify_value`
  (`expand(radsimp(expand(x)))`; rationals unchanged), which rationalizes denominators (`sqrt(3)/3` instead of
  `1/sqrt(3)`).
- Linear dependence is decided numerically (relative tolerance `1e-8`); the kept vectors are then orthonormalized
  exactly.
- sympy re-evaluates some products when a pickle is loaded (`(a+b)/2` → `a/2 + b/2`), so "binary is up to date"
  tests compare after a pickle round trip.

---

### 8. Binary data for MultiPie (`binary_data`)

Code: `util_binary_data`, `create_binary/create_binary_data.py`.

1. **Layout.** As MultiPie's `multipie/binary_data`: `info.pkl` (the same data as `binary/info.pkl`, with the layout
   added to the header) and one file for each group, `{PG,SG,MPG,MSG}/{id}.pkl` with `:` of the id replaced by `_`, e.g. `PG/PG_32.pkl` (47, 230, 157, 1651 files; git on Windows cannot check out a file name with `:`). No
   separate `_opt` files: a PG file is at most about 260 KB.
2. **Contents.** Regrouped from the binaries in `multipiegen/binary`; nothing is recomputed except the representation
   matrices of PG (`util_representation_matrix`, below 1 KB each).
   - PG: `info`, `symmetry_operation`, `wyckoff`, `character`, `representation_matrix`, `harmonics`,
     `multipolar_harmonics`, `atomic_samb`, `cluster_samb` (real PG only), `coset` (real PG with a parent only;
     `"first"` of `coset.pkl`).
   - SG: `info`, `symmetry_operation`, `wyckoff`, `cluster_samb`, `coset` (SG with a parent only).
   - `info.pkl` has also `"parent"` and `"coset"` (`"parent"` and `"coset"` of `coset.pkl`).
   - MPG: `info`, `symmetry_operation`, `wyckoff`, `active_multipole`, `response_tensor`.
   - MSG: `info`, `symmetry_operation`, `wyckoff`.
3. **Not included.** The intermediate cluster bases (`cluster_scalar_basis`, `cluster_vector_basis`,
   `cluster_bond_basis`) and the spherical atomic multipoles (`atomic_multipole`).
4. **Harmonics.** Not repeated: the SAMB and the multipolar harmonics are looked up by `(X, l, Gamma, n)` in
   `"harmonics"` of the PG file (T -> Q, M -> G); for SG and MPG, the PG file is that of `info.PG`.
5. **Header.** Each file has its own header (`str()`), made from the headers of the source binaries, so a file can be
   read with `BinaryManager` alone.

---

### 9. Cosets (`coset.pkl`)

Code: `util_coset`, `util_cluster_samb_image`, `create_binary/create_coset.py`. Range: 35 real point groups and 216
space groups with a parent. Groups and operations are taken from the tables of this repository only
(`data_parent_pg/sg`, `data_wyckoff_site`, `group_symmetry_operation`); no table is stored.

1. **Cosets.** The parent `P` is split into left cosets `g H` (for SG modulo the lattice of `H`, so that the
   translations of `P` that are not lattice translations of `H` give separate cosets), in the order: rotation part
   the identity first, then by the length of the translation. The chain of parents follows
   `data_representative_parent_pg/sg`, with the composed transformations `x_P = R x_H`.
2. **Images.** For each coset of the direct parent, the image `g H g^-1` is looked for among the children of `P` with
   the same axes (for SG up to an origin shift `k`; among the shifts, the one that maps the largest number of Wyckoff
   positions onto themselves, then the shortest). The map `M = R'^-1 k^-1 g R` maps each Wyckoff position of `H` to one
   of the image, with the map of its parameters.
3. **First coset.** The second coset of `data_parent_pg/sg[H][P]` (the first one other than `H`), with its
   representative (first operation) and the order `n`. Its kind is `"rotation"` (det > 0) or `"inversion"`.
4. **Images of the sites, bonds, and cluster SAMB (inversion).** Each Wyckoff position `W` of `H` is paired with the
   Wyckoff position `W'` of `H'` onto which `M` maps its orbit (same orbit size and number of free parameters; lattice
   and centering translations and the reversal of a bond allowed): row `i` of `W` goes to row `perm[i]` of `W'`. The
   cluster SAMB moved in this way (the value on site `i` goes to site `perm[i]`, with the sign of the bond reversal for
   `T/M`) are expanded in the cluster SAMB of `W'` within each irrep (`Q` and `T/M` separately for bonds): `"eps"`
   (`+-1` for each component) when this is diagonal, otherwise `"U"` (exact `+-sqrt(rational)` where possible).
   Among the equivalent correspondences (`perm`, `sign`; they differ by an operation of `H`), the first one mapping
   the SAMB label by label is taken; for a Wyckoff position mapped onto itself without such one, the one whose `U`
   are involutions (`U^2 = 1`, eigenvalues `+-1`) in every block, with the fewest blocks other than `+-1`.
   Results for the 77 groups (1847 Wyckoff positions): 1809 label by label (20 of them after the change of the
   reference sites, 5.8), 38 with `U` (92 blocks other than `+-1`). All `U` are involutions except 3 blocks of SG 80
   (I4_1; `E` of `4a@4a` (T/M) and of `8c@4a` (Q, T/M)): there `U` rotates the two components of the real 2D irrep
   `E` of complex type of C4 by 90 degrees (`U^2 = -1`, eigenvalues `+-i`), whatever correspondence is taken.
5. **Parent structure.** For each Wyckoff position with free parameters, the largest affine subspaces of the
   parameters at which `c p` lies on the orbit of `H` for every representative `c` of the right cosets `H c` (then the
   orbit is invariant under `P`), with the orbit of the multiplicity of `W` (no coinciding sites). The limits found with
   the first search window come first; others are added unless their orbit of `H` is that of one listed before.
6. **Header and keys.** See the header of `coset.pkl` (`str(load_coset())`); `binary_data` has `"first"` as `coset`
   of each group file, and `"parent"` and `"coset"` in `info.pkl`.
