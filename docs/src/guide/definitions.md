---
orphan: true
---

# Definitions

This page defines the quantities provided by MultiPieGen. How to get the actual function forms and $U$
matrices from the data is described in a later section. The harmonics and the SAMB are labeled by the multipole key
`PGMultipoleType` and written with its LaTeX symbol, e.g. $\mathbb{Q}_{1}(T_{1u})$; both are explained in the
last section, [Multipole key](#multipole-key). A Japanese version is [here](definitions_ja.md).

## 1. Point-group harmonics

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

## 2. PG multipolar harmonics

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

## 3. Atomic orbitals

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

## 4. Atomic SAMB

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

## 5. Cluster SAMB

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

## 6. Response tensors

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

(multipole-key)=
## 7. Multipole key (`PGMultipoleType`)

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

### LaTeX symbol

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
