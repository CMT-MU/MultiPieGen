---
orphan: true
---

# Multipole key (`PGMultipoleType`)

This page explains the key of the multipoles, which is needed to use the harmonics and the SAMB. The quantities themselves are defined in [Definitions](definitions.md). A Japanese version is [here](multipole_key_ja.md).

Every multipole of a point group (crystal harmonics, multipolar harmonics, atomic SAMB, cluster SAMB) is labeled by
the named tuple `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` of the `multipie` package:

| field | meaning | values |
|---|---|---|
| `X` | type, by the parity under the spatial inversion ($\mathcal{P}$) and the time reversal ($\mathcal{T}$) | `Q` electric ($\mathcal{P} = (-1)^l$, $\mathcal{T}$ even), `G` electric toroidal ($(-1)^{l+1}$, even), `T` magnetic toroidal ($(-1)^l$, odd), `M` magnetic ($(-1)^{l+1}$, odd) |
| `l` | rank | 0, 1, 2, ... |
| `Gamma` | irrep of the point group | name in the character table, e.g. `A1g`, `Eu`, `T2g` |
| `n` | multiplicity of the crystal harmonics `(X, l, Gamma)`: which of the harmonics with the same `(X, l, Gamma)` | 1, 2, ...; `-1` if unique |
| `p` | extra multiplicity, numbering entries that share all the other fields | 1, 2, ...; `-1` if unique |
| `s` | internal rank (rank of the internal degree of freedom) | 0 (scalar), 1, 2, 3 |
| `k` | internal component: the orbital part has rank `l + k` | `-s, ..., s` |
| `x` | internal type: parity of the internal basis | `q` (polar), `g` (axial) |

- `(X, l, Gamma, n)` always identifies a crystal harmonics (its xyz polynomial and $U$); `T` and `M` use the
  harmonics of `Q` and `G`. `n` is never renumbered, so different data with the same `(X, l, Gamma, n)` refer to the
  same harmonics.
- `p` appears only where several entries share `(X, l, Gamma, n, s, k, x)`, e.g. two cluster SAMB of the same
  harmonics on a Wyckoff position.
- Scalar quantities (crystal harmonics, cluster SAMB) have `s = 0, k = 0, x = "q"`.
- A key names a whole multiplet. Its value holds all components together, in the order of the components of the
  irrep (rows of the array).

Examples:

| key | meaning |
|---|---|
| `("Q", 1, "T1u", -1, -1, 0, 0, "q")` | electric dipole of Oh, components $(x, y, z)$ |
| `("G", 0, "A1u", -1, -1, 1, 1, "g")` (s = 1) | electric toroidal monopole made of the rank-1 orbital part and an axial vector |
| `("T", 1, "Eu", -1, 1, 0, 0, "q")` | an anti-symmetric bond cluster SAMB of D4h with the symmetry of $(x, y)$ |
