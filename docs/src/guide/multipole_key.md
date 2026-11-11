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

## LaTeX symbol

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
| $[x]$, $\vec{\ }$ | `x`, `s` | for the multipolar harmonics (`s > 0`): an arrow on $\mathbb{X}$ and the internal type $[q]$ or $[g]$ |
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
| `("G", 0, "A1u", -1, -1, 1, 1, "g")` | $\vec{\mathbb{G}}_{0}^{(1,1)}[g](A_{1u})$ | multipolar harmonics: electric toroidal monopole from the rank-1 orbital part and an axial vector |
| `("M", 1, "T1g", -1, -1, 1, -1, "q")` | $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ | atomic SAMB: spin magnetic dipole |
| `("T", 1, "Eu", -1, 1, 0, 0, "q")` | $\mathbb{T}_{1}^{(b)}(E_{u},a)$ | bond cluster SAMB: an anti-symmetric bond function of D4h with the symmetry of $(x, y)$ ($p = 1$) |
| `("Q", 2, "E", 1, 2, 0, 0, "q")` | $\mathbb{Q}_{2}(E,1b)$ | $n = 1$ and $p = 2$ |
