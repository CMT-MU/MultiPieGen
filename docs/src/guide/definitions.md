---
orphan: true
---

# Definitions

This page defines the quantities provided by MultiPieGen. How to get the actual function forms, $U$, and
matrices from the data is described in a later section. The harmonics and the SAMB are labeled by the multipole key
`PGMultipoleType`, explained in [Multipole key](multipole_key.md). A Japanese version is [here](definitions_ja.md).

## 1. Crystal harmonics

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

**Crystal harmonics.** A crystal harmonics `(X, l, Gamma, n)` of a point group is a set of rank-`l` functions that
forms a basis of the irrep `Gamma`. Its component `gamma` is

$$
X_{l\Gamma n}^{(\gamma)}(\boldsymbol{r}) = \sum_{m} U_{m\gamma}\, O_{lm}(\boldsymbol{r}),
$$

with $U$ (`(2l+1) x dim`, rows `m = l, ..., -l`), whose columns are orthonormal.

*Example (Oh).* The dipole $Q_1(T_{1u}) = (x, y, z)$ has

$$
U = \begin{pmatrix} -\frac{1}{\sqrt{2}} & \frac{i}{\sqrt{2}} & 0 \\ 0 & 0 & 1 \\ \frac{1}{\sqrt{2}} & \frac{i}{\sqrt{2}} & 0 \end{pmatrix}
\quad (m = 1, 0, -1),
$$

e.g. $x = (-O_{11} + O_{1,-1})/\sqrt{2}$, and the quadrupole $Q_2(E_g) = \left(z^2 - \frac{x^2+y^2}{2},
\frac{\sqrt{3}}{2}(x^2 - y^2)\right)$.

- **Top groups.** The harmonics are defined once for the two top groups, Oh (cubic) and D6h (trigonal and
  hexagonal), up to rank 11. Every one of the 47 point groups (32 crystallographic, 5 other settings `-1`, 10
  complex `-c`) takes its harmonics from them: each component is a component of a top harmonics times a sign,
  $X_{\rm PG} = c\, X_{\rm top}$ with $c = \pm 1$.
- The representation matrices of the irreps follow from $U$, and all harmonics of an irrep give the same matrices.
- The irreps follow the order of the character table. Complex irreps (in the `-c` groups) end in `a`/`b` and have a
  single component.
- *Multiplicity.* When an irrep appears more than once at the same rank, `n = 1, 2, ...` distinguishes the
  harmonics, e.g. C4h: $Q_2(B_g)$ with `n = 1`: $\frac{\sqrt{3}}{2}(x^2 - y^2)$ and `n = 2`: $\sqrt{3}xy$.
- *Complex group.* E.g. C4-c: $Q_1(E_a) = (x + iy)/\sqrt{2}$ and $Q_1(E_b) = (x - iy)/\sqrt{2}$.
- Key: `PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q")` (see [Multipole key](multipole_key.md)).

**Polar (`Q`) and axial (`G`) harmonics.** An axial harmonics transforms with an extra factor `det R` under an
improper operation `R` (a mirror, the inversion, or a rotoinversion). Both come from the same top harmonics, but in a
point group a polar and an axial harmonics of the same function can belong to different irreps, and in the same irrep
their components can be in another order or have other signs (the coefficient `c` and the component assigned differ
between `Q` and `G`). So, for the same `(l, Gamma, n)`, the functions of `Q` and `G` differ in general:

- In the top groups Oh and D6h, the axial harmonics is the polar one with `g` and `u` interchanged in the irrep: the
  same function, components, and order ($c = 1$). E.g. Oh: $Q_1(T_{1u}) = G_1(T_{1g}) = (x, y, z)$.
- In a group of proper rotations only (C1, C2, D2, C4, D4, C3, D3, D3-1, C6, D6, T, O), there is no improper
  operation, and `Q` and `G` are identical. E.g. C4: $Q_1(E) = G_1(E) = (x, y)$.
- In the other groups, they differ. Examples:

| group | polar | axial | difference |
|---|---|---|---|
| C4v | $Q_1(A_1) = z$ | $G_1(A_2) = z$ | same function, another irrep |
| C4v | $Q_1(E) = (x, y)$ | $G_1(E) = (-y, x)$ | components interchanged, with a sign |
| C4h | $Q_1(E_u) = (x, y)$ | $G_1(E_g) = (y, x)$ | components interchanged |
| D2d | $Q_1(E) = (x, y)$ | $G_1(E) = (x, -y)$ | sign of a component |
| C4v | $Q_2(B_1) = \frac{\sqrt{3}}{2}(x^2 - y^2)$ | $G_2(B_1) = \sqrt{3}xy$ | another function for the same irrep |
| D2d | $Q_2(A_1) = z^2 - \frac{x^2+y^2}{2}$ | $G_2(A_1) = \frac{\sqrt{3}}{2}(x^2 - y^2)$ | another function for the same irrep |

## 2. Multipolar harmonics

A multipolar harmonics combines a crystal harmonics with an internal degree of freedom of rank `s` (`s = 1`: a
vector such as a displacement or a spin; `s = 2, 3`: higher-rank tensors). The multipolar spherical harmonics of
total rank `l` couple the orbital part of rank `l + k` with the internal basis $e_\nu$ of rank `s`:

$$
O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} \sqrt{2l+1} \sum_{\nu}
\begin{pmatrix} l+k & l & s \\ m-\nu & -m & \nu \end{pmatrix} O_{l+k, m-\nu}\, e_\nu,
\qquad k = -s, \dots, s,
$$

and the multipolar harmonics of a point group is $\sum_m U_{m\gamma} O_{lm}^{(s)}(k)$ with $U$ of the
crystal harmonics `(X, l, Gamma, n)`.

- The internal basis $e_\nu$ is a set of rank-`s` orbitals (e.g. `px, py, pz` for `s = 1`), of the cubic family
  (Oh) or of the hexagonal family (D6h).
- Internal type: `x = q` if the internal basis is polar, `x = g` if axial. A multipole has `x = q` if
  (`X = Q` and `s + k` even) or (`X = G` and `s + k` odd), and `x = g` otherwise.
- Range: `s = 1, 2, 3`, `l <= 3`, all 47 point groups. Key: `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` (see
  [Multipole key](multipole_key.md)).
- The value is the coefficient of each internal basis function (polynomials in x, y, z) and the resulting expression
  in terms of the internal basis symbols.

*Examples (Oh, `s = 1`).* With a polar vector $\boldsymbol{Q} = (Q_x, Q_y, Q_z)$ or an axial vector
$\boldsymbol{G} = (G_x, G_y, G_z)$ as the internal basis:

| key `(X, l, Gamma, n, p, s, k, x)` | expression | meaning |
|---|---|---|
| `(Q, 0, A1g, -1, -1, 1, 1, q)` | $(xQ_x + yQ_y + zQ_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{Q}$, electric monopole |
| `(Q, 1, T1u, -1, -1, 1, -1, q)` | $(Q_x, Q_y, Q_z)$ | the polar vector itself |
| `(G, 1, T1g, -1, -1, 1, 0, q)` | $(zQ_y - yQ_z, xQ_z - zQ_x, yQ_x - xQ_y)/\sqrt{2}$ | $\boldsymbol{Q} \times \boldsymbol{r}$, electric toroidal dipole |
| `(G, 0, A1u, -1, -1, 1, 1, g)` | $(xG_x + yG_y + zG_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{G}$, electric toroidal monopole |

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

1. **Atomic multipoles.** The spherical multipole operators $X_{lm}^{(s)}(k)$ of `X = Q, G, T, M`: `s = 0`
   (orbital) and `s = 1` (spin-dependent), with `k = -s, ..., s`. Their matrix elements between orbitals of ranks
   `L1` and `L2` follow from the Wigner-Eckart theorem. `Q` and `G` are even and `T` and `M` are odd under time
   reversal.
2. **Symmetry adaptation.** For a crystal harmonics with $U$,
   $X_{l\Gamma n}^{(s)}(k)[\gamma] = \sum_m U_{m\gamma} X_{lm}^{(s)}(k)$, using polar harmonics for `Q, T` and
   axial harmonics for `G, M`.
3. **Normalization.** Each component has the Frobenius norm 1 in a diagonal block (`L1 = L2`), and `1/sqrt(2)` in
   an off-diagonal block (norm 1 together with its hermitian conjugate block `(L2, L1)`).

- Bases: `lg` (`s = 0` only), `lgs`, and `jml`.
- Key: `PGMultipoleType(X, l, Gamma, n, -1, s, k, "q")` (see [Multipole key](multipole_key.md)). Value: `(dim, n1, n2)`, the matrix of each component for
  the `n1` orbitals of rank `L1` and the `n2` orbitals of rank `L2`.
- Order: `Q, G, T, M`, then `s`, `k`, `l`, `Gamma`, `n`. All 47 point groups.

*Examples (Oh).* For the p orbitals (`lg`, block `(1, 1)`, basis `px, py, pz`):

- $Q_0(A_{1g}) = \frac{1}{\sqrt{3}}\,\mathrm{diag}(1, 1, 1)$, the charge;
- $Q_2(E_g) = \left(\frac{1}{\sqrt{6}}\,\mathrm{diag}(-1, -1, 2),\ \frac{1}{\sqrt{2}}\,\mathrm{diag}(1, -1, 0)\right)$, the
  orbital quadrupole;
- $M_1(T_{1g})$, `s = 0`: $(l_x, l_y, l_z)/\sqrt{2}$, the orbital angular momentum, e.g. the `x` component has
  $\langle p_y | \cdot | p_z \rangle = -i/\sqrt{2}$.

For the s orbital with spin (`lgs`, block `(0, 0)`, basis `(s, u), (s, d)`): $Q_0(A_{1g})$, `s = 0`, is
$\frac{1}{\sqrt{2}}\,\mathbb{1}$, and $M_1(T_{1g})$, `s = 1`, `k = -1`, is $(\sigma_x, \sigma_y, \sigma_z)/\sqrt{2}$, the
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
- Key: `PGMultipoleType(X, l, Gamma, n, p, 0, 0, "q")` (see [Multipole key](multipole_key.md)). `(X, l, Gamma, n)` names the crystal harmonics that the
  function corresponds to (`T` and `M` with those of `Q` and `G`), and `p` distinguishes functions sharing it.
- Value: `(dim, N)`, the values on the `N` sites (bonds) for each component; each component has norm 1.
- Order: `Gamma` (character table), and in `bond`, `Q` before `T/M` in each `Gamma`.

*Examples (D4h).*

- Site `2a`: two sites at $(0, 0, \pm z)$. $Q_0(A_{1g}) = (1, 1)/\sqrt{2}$ (uniform) and $Q_1(A_{2u}) = (1, -1)/\sqrt{2}$
  (a dipole along `z`).
- Bond `2b@1o`: two bonds through the origin, along `x` and along `y`.
  $Q_0(A_{1g}) = (1, 1)/\sqrt{2}$, $Q_2(B_{1g}) = (-1, 1)/\sqrt{2}$ (real hoppings), and
  $T_1(E_u) = \left((i, 0), (0, i)\right)$ (imaginary hoppings): the `x` component lives on the bond along `x`, the `y`
  component on the bond along `y`.

## 6. Response tensors

A response tensor of rank 0 to 4 relates two physical quantities, e.g. a current and an electric field. For each of
the 157 magnetic point groups, MultiPieGen gives the tensors allowed by symmetry.

- Type `X` by the parity under the spatial inversion and the time reversal: `Q` (polar, time-reversal even), `G`
  (axial, even), `T` (polar, odd), `M` (axial, odd).
- `opt` specifies the part of the tensor by the symmetry of its indices: `s` symmetric and `a` anti-symmetric
  in a pair of indices, e.g. `s`/`a` for rank 2, `ss`/`sa`/`a` for rank 3, and `sss`/`ssa`/`aas`/`aaa`/`sa`/`as`/
  `s`/`a`/`t` for rank 4.
- Each tensor component is decomposed into multipoles (irreducible cartesian parts). Only the active multipoles of
  the group, those invariant under all its operations including the time reversal, can be nonzero. This gives the
  nonvanishing components and the relations among them.
- Value for `(X, rank, opt)`: the matrix of the components (symbols such as $Q_{xyz}$, `0` for vanishing ones;
  rows and columns are index pairs or triples, e.g. `[xx, yy, zz, yz, zx, xy]`), and the expression of each component
  in terms of the active multipoles.

*Examples.*

- 4/mmm (`MPG:15.1.53`), `(Q, 2, s)` (e.g. the electric conductivity): $\mathrm{diag}(Q_{xx}, Q_{xx}, Q_{zz})$, with
  $Q_{xx} = \mathbb{Q}_0^{(1)} - \mathbb{Q}_u^{(1)}$ and $Q_{zz} = \mathbb{Q}_0^{(1)} + 2\mathbb{Q}_u^{(1)}$.
- 4/mm'm' (`MPG:15.6.58`, a ferromagnet along `z`), `(T, 2, a)` (e.g. the anomalous Hall effect): the only nonzero
  components are $T_{xy} = -T_{yx} = \mathbb{M}_z^{(1)}$, given by the active magnetic dipole $M_z$.
