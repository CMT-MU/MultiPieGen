---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Examples

Sample code for each feature, in `docs/examples`. Run each file from anywhere after `pip install -e .`, e.g.
`python docs/examples/group.py`. The outputs below are made when the documentation is built.

## cluster_basis.py

Load the orthonormal (scalar) Oh/D6h cluster basis (cluster_scalar_basis.pkl), keyed by PGMultipoleType, expand a function on the cluster.

```{code-cell} ipython3
import numpy as np
from multipiegen.util.util_cluster_basis import load_cluster_basis, cluster_basis_matrix
```

```{code-cell} ipython3
# 1. load the cluster basis of Oh and D6h.
data = load_cluster_basis()
print("groups:", list(data.keys()))

oh = data["Oh"]
print("site0 :", oh["site0"])
print("number of cluster points:", len(oh["cluster"]))
print("first 3 cluster points (exact):")
for p in oh["cluster"][:3]:
    print("   ", p)
```

```{code-cell} ipython3
# 2. basis multiplets are keyed by PGMultipoleType(X="Q", l, Gamma, n, p=-1, s=0, k=0, x="q") (MultiPie Ver.2 site SAMB);
#    the value has all components, (dim, N).
print("\nbasis keys of Oh:")
for key in oh["basis"].keys():
    print("   ", key, oh["basis"][key].shape)

key = [k for k in oh["basis"] if k.l == 1 and k.Gamma == "T1u"][0]
print(f"\n{key}, x component, first 4 values (exact):", oh["basis"][key][0][:4])

# harmonics dropped by Gram-Schmidt (linearly dependent on lower ranks).
print("number of dropped harmonics:", len(oh["dropped"]), ", e.g.", oh["dropped"][:3])
```

```{code-cell} ipython3
# 3. basis as a float matrix (each row is a component, columns follow oh["cluster"]); keys are (PGMultipoleType, comp).
m, keys = cluster_basis_matrix(oh)
print("\nbasis matrix shape:", m.shape)
print("orthonormal ?", np.allclose(m @ m.T, np.eye(len(keys))))
```

```{code-cell} ipython3
# 4. expand a function on the cluster in the basis.
#    f is a weight on the cluster points, e.g. f(r) = x + 0.5 x y z at each point.
r = np.array(oh["cluster"], dtype=float)
f = r[:, 0] + 0.5 * r[:, 0] * r[:, 1] * r[:, 2]

c = m @ f  # coefficients, f = sum_k c_k e_k.
print("\nexpansion of f = x + 0.5 xyz:")
for k, ck in zip(keys, c):
    if abs(ck) > 1e-8:
        print(f"    {k}: {ck:.6f}")
print("reconstructed ?", np.allclose(c @ m, f))

# weight of each irrep (sum of |c_k|^2 over the irrep), total equals |f|^2.
weight = {}
for k, ck in zip(keys, c):
    weight[k[0].Gamma] = weight.get(k[0].Gamma, 0.0) + ck**2
print("irrep weights:", {ir: round(float(w), 6) for ir, w in weight.items() if w > 1e-8})
print("|f|^2 =", round(float(f @ f), 6))
```

```{code-cell} ipython3
# 5. the same for D6h (24 points).
d6h = data["D6h"]
m6, keys6 = cluster_basis_matrix(d6h)
print("\nD6h basis matrix shape:", m6.shape, ", orthonormal ?", np.allclose(m6 @ m6.T, np.eye(len(keys6))))
print("D6h irreps:", sorted({k[0].Gamma for k in keys6}))
```

## cluster_vector_basis.py

Load the orthonormal vector cluster basis (internal rank s=1, polar/axial internal basis) of the Oh/D6h clusters, keyed by PGMultipoleType, expand a vector field.

```{code-cell} ipython3
import numpy as np

from multipiegen.util.util_cluster_vector_basis import cluster_vector_basis_matrix, load_cluster_vector_basis

data = load_cluster_vector_basis()

for group in ["Oh", "D6h"]:
    d = data[group]
    print(f"{group}: {len(d['cluster'])} points, internal basis {d['internal_basis']}")
    for x in ["q", "g"]:
        m, keys = cluster_vector_basis_matrix(d, x)
        n_mult = len([k for k in d["basis"] if k.x == x])
        n_drop = len([k for k in d["dropped"] if k.x == x])
        ok = np.allclose(m @ m.T, np.eye(len(m)))
        print(f"  x={x}: {len(keys)} vectors in {n_mult} multiplets, orthonormal: {ok}, dropped multiplets: {n_drop}")
```

```{code-cell} ipython3
# keys: PGMultipoleType(X, l, Gamma, n, p=-1, s=1, k, x) as MultiPie Ver.2 vector SAMB; value: (dim, N, 3) at the cluster points.
d = data["Oh"]
for key in list(d["basis"])[:4]:
    print(key, "at the first point:", d["basis"][key][:, 0].tolist())

# a vector field on the cluster expanded in the basis: coefficients = <basis|field>.
P = np.array(d["cluster"], dtype=float)
field = np.cross(P, [0, 0, 1])  # rotation around z (axial dipole like).
m, keys = cluster_vector_basis_matrix(d, "q")
c = m @ field.ravel()
print("nonzero components of r x e_z:", [(keys[i], round(c[i], 6)) for i in range(len(c)) if abs(c[i]) > 1e-9])
```

## bond_cluster_basis.py

Load the orthonormal bond cluster basis of the Oh/D6h clusters (vector cluster basis contracted with the bond vectors and re-orthonormalized), keyed by PGMultipoleType as MultiPie Ver.2's bond_a.

```{code-cell} ipython3
import numpy as np

from multipiegen.util.util_bond_cluster_basis import bond_cluster_basis_matrix, load_bond_cluster_basis

data = load_bond_cluster_basis()

for group in ["Oh", "D6h"]:
    d = data[group]
    print(f"{group}: {len(d['cluster'])} points, vector0 = {d['vector0']}, bond b_0 = {d['bond'][0]}")
    for x in ["q", "g"]:
        m, keys = bond_cluster_basis_matrix(d, x)
        ok = np.allclose(m @ m.T, np.eye(len(m)))
        n_mult = len([k for k in d["basis"] if k.x == x])
        n_drop = len([k for k in d["dropped"] if k.x == x])
        print(f"  x={x}: {len(keys)} vectors in {n_mult} multiplets, orthonormal: {ok}, dropped vector multiplets: {n_drop}")
```

```{code-cell} ipython3
# keys: PGMultipoleType(X, l, Gamma, n, p, s=0, k=0, x) as MultiPie Ver.2 bond_a; value: (dim, N) values on the bonds.
# source: the multiplet of the vector cluster basis (PGMultipoleType with s=1).
d = data["Oh"]
for key in list(d["basis"])[:4]:
    print(key, d["basis"][key].shape, "<-", d["source"][key])
print("with add_multiplicity:", [k for k in d["basis"] if k.p != -1][:2])
print("dropped vector cluster basis multiplets (first 2):", d["dropped"][:2])
```

## cluster_samb.py

Cluster SAMB of each Wyckoff site and bond (site; bond = symmetric Q and anti-symmetric T/M) of the real point groups and the space groups (cluster_samb.pkl).

```{code-cell} ipython3
import numpy as np

from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_harmonics import harmonics_expression
from multipiegen.util.util_info import load_info

info = load_info()
data = load_cluster_samb()
```

```{code-cell} ipython3
# {PG_id/SG_id: {"site"/"bond": {wyckoff: {PGMultipoleType: samb (dim, N)}}}}.
# "bond": symmetric bond (X = Q, real) and i times anti-symmetric bond (X = T/M, pure imaginary).
# the expression of the harmonics (X, l, Gamma, n) is looked up by the key in harmonics.pkl (T, M as Q, G).
d = data[info["id"]["D4h"]]
for kind in ["site", "bond"]:
    w = list(d[kind].keys())[1]
    print(f"* {kind} {w}")
    for key, samb in d[kind][w].items():
        print("  ", key, "ex =", list(harmonics_expression("D4h", key)))
        print("  ", samb.tolist())

# orthonormality on the N sites (bonds) of each Wyckoff position, for Q and T/M separately (different time-reversal
# parity).
d = data["SG:225"]
for kind in ["site", "bond"]:
    for X in ["Q", "TM"]:
        parts = [[np.array(s, dtype=complex) for k, s in samb.items() if k.X in X] for samb in d[kind].values()]
        parts = [np.vstack(p) for p in parts if p]
        if not parts:
            continue
        ok = all(np.allclose(v.conj() @ v.T, np.eye(len(v))) for v in parts)
        print(f"SG:225 {kind} {X}: {len(parts)} Wyckoff positions, orthonormal: {ok}")
```

## group.py

Data of each group (group.pkl): symmetry operations, product, mapping to Oh/D6h, character, Wyckoff site and bond, active multipoles.

```{code-cell} ipython3
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

info = load_info()
group = load_group()
print("number of groups:", len(group.keys()), "(PG, SG, MPG, MSG)\n")
```

```{code-cell} ipython3
# 1. symmetry operations of a point group (D4h).
d4h = group[info["id"]["D4h"]]
so = d4h["symmetry_operation"]
print("D4h:", d4h["info"].international, ", SO:", so["tag"])
print("  matrix of 4+[001] (cartesian):", so["cartesian"][so["tag"].index("4+[001]")].tolist())
print("  product 4+[001] * 2[100] =", so["product"][("4+[001]", "2[100]")])
print("  dipole matrix of 4+[001], axial (G, s=1):", so["so_matrix"][("G", 1)][2].tolist())
```

```{code-cell} ipython3
# 2. mapping to the Oh cluster: mapping[0] is D4h itself, mapping[k] is a left coset in Oh.
top = so["mapping"]["top"]
oh_tag = group[info["id"][top]]["symmetry_operation"]["tag"]
print(f"\nmapping to {top}:")
for k, row in enumerate(so["mapping"]["mapping"]):
    print(f"  coset {k}:", [oh_tag[i - 1] for i in row][:6], "...")
```

```{code-cell} ipython3
# 3. character table (complex PG has its own table, other data are those of the real PG).
c3h = group[info["id"]["C3h-c"]]
print("\nC3h-c irreps:", list(c3h["character"]["table"]))
print("  character of Ea' for all SO:", c3h["character"]["table_full"]["Ea'"].tolist())
```

```{code-cell} ipython3
# 4. Wyckoff site and bond of a space group (P6_3/mmc).
sg = group["SG:194"]
site = sg["wyckoff"]["site"]["2a"]
print("\nSG:194", sg["info"].international, "Wyckoff 2a:", site["conventional"].tolist(), ", symmetry:", site["symmetry"])
print("  mapping:", site["mapping"])
print("  bonds at 2a:", site["bond"][:4])
bw = site["bond"][0]
bond = sg["wyckoff"]["bond"][bw]
print(f"  bond {bw}: symmetry {bond['symmetry']}, (vector, center) =", bond["conventional"].tolist())
```

```{code-cell} ipython3
# 5. magnetic point group: active multipoles.
mpg = group["MPG:24.5.149"]
print("\nMPG", mpg["info"].international, "time-reversal sign:", mpg["symmetry_operation"]["tr_sign"].tolist())
print("  active multipoles:", mpg["active_multipole"])
```

## representation_matrix.py

Representation matrices of a point group from the harmonics up to rank 11, and the check that all harmonics of an irrep give the same matrices.

```{code-cell} ipython3
from multipiegen.util.util_representation_matrix import (
    check_representation_matrix,
    harmonics_blocks,
    representation_matrix,
)

pg = "D3h"
```

```{code-cell} ipython3
# 1. harmonics blocks (X, rank, irrep, multiplicity) -> [(coeff, top key)] of each component.
blocks = harmonics_blocks(pg)
print(f"{pg}: {len(blocks)} harmonics blocks up to rank 11, e.g.")
for b in list(blocks)[:4]:
    print("  ", b, blocks[b])
```

```{code-cell} ipython3
# 2. all harmonics of the same irrep have the same representation matrices.
print("\nirrep: number of harmonics blocks, max deviation from the first one, max closure residual")
for irrep, d in check_representation_matrix(pg).items():
    print(f"  {irrep:4s}: {len(d['blocks']):3d}, {d['deviation']:.1e}, {d['closure']:.1e}")
```

```{code-cell} ipython3
# 3. exact representation matrices (in the order of symmetry operations).
rm = representation_matrix(pg)
D = rm["matrix"]["E'"]
for tag, d in list(zip(rm["tag"], D))[:3]:
    print(f"\nE' {tag}:\n{d}")
```

## response_tensor.py

Response tensors (rank 0-4) of a magnetic point group in terms of active multipoles (response_tensor.pkl).

```{code-cell} ipython3
from multipiegen.util.util_group import load_group
from multipiegen.util.util_response_tensor import load_response_tensor

rt = load_response_tensor()
group = load_group()

mpg = "MPG:24.5.149"
print(mpg, group[mpg]["info"].international, ", active multipoles:", group[mpg]["active_multipole"])
```

```{code-cell} ipython3
# 1. keys: (X, rank, opt), X = Q (polar), G (axial), T (magnetic polar), M (magnetic axial).
d = rt[mpg]
print("number of tensors:", len(d), ", e.g.", list(d)[:4])
```

```{code-cell} ipython3
# 2. non-vanishing tensors and their expressions in terms of active multipoles.
for (X, rank, opt), (m, ex) in d.items():
    if ex and rank <= 3:
        print(f"\n{X} rank {rank} {opt}:\n{m}")
        for c, e in ex.items():
            print(f"  {c} = {e}")
```

## harmonics.py

Expression and u-matrix of the harmonics of a point group (harmonics.pkl), looked up by a multipole key of any binary.

```{code-cell} ipython3
from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_harmonics import harmonics_expression, harmonics_u_matrix, load_harmonics
from multipiegen.util.util_info import load_info

info = load_info()
```

```{code-cell} ipython3
# 1. {PG_id: {PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q"): (expression (dim), u_matrix (2l+1, dim))}}.
d = load_harmonics()[info["id"]["D4h"]]
print("D4h:", len(d), "harmonics (rank <= 11)")
for key, (ex, U) in list(d.items())[:4]:
    print(" ", key, list(ex), U.shape)
```

```{code-cell} ipython3
# 2. look up the harmonics of a key of another binary (only (X, l, Gamma, n) of the key is used).
samb = load_cluster_samb()[info["id"]["D4h"]]["site"]["4b"]
for key in list(samb)[:3]:
    print(key, "->", list(harmonics_expression("D4h", key)))
key = list(samb)[1]
print("u-matrix of", key, ":\n", harmonics_u_matrix("PG:15", key))
```

## multipolar_harmonics.py

Multipolar harmonics (internal rank s=1,2,3, rank l<=3) with polar/axial internal basis of the data_basis orbitals (multipolar_harmonics.pkl).

```{code-cell} ipython3
import collections

from multipiegen.util.util_harmonics import harmonics_expression
from multipiegen.util.util_harmonics_multipole import load_multipolar_harmonics
from multipiegen.util.util_info import load_info

info = load_info()
hm = load_multipolar_harmonics()
```

```{code-cell} ipython3
# 1. keys: PGMultipoleType(X, l, Gamma, n, p, s, k, x); value: (internal_coefficient, multipole_expression).
#    the expression of the harmonics (X, l, Gamma, n) is looked up by the key in harmonics.pkl.
for pg in ["Oh", "D6h"]:
    d = hm[info["id"][pg]]
    print(pg, len(d), "multipolar harmonics:", dict(collections.Counter((k.x, k.s) for k in d)))
    crystal = "cubic" if pg == "Oh" else "hexagonal"
    print("  internal basis (s=3, q):", info["harmonics"]["internal_basis"][crystal][(3, "q")])
```

```{code-cell} ipython3
# 2. an example: a dipole-type (s=1) multipole of rank 2 in Oh.
d = hm[info["id"]["Oh"]]
key = [k for k in d if k.s == 1 and k.l == 2 and k.Gamma == "Eg" and k.k == 1][0]
coeff, expr = d[key]
print("\n", key)
print("  harmonics:", harmonics_expression("Oh", key).tolist())
print("  internal coefficient of [Q_x, Q_y, Q_z] (first component):", coeff[0].tolist())
print("  multipole expression (first component):", expr[0])
```

```{code-cell} ipython3
# 3. axial internal basis (g) of rank s=2 in D6h.
d = hm[info["id"]["D6h"]]
key = [k for k in d if k.x == "g" and k.s == 2][0]
print("\n", key, "\n  ", d[key][1].tolist())
```

## atomic_multipole.py

Atomic multipole matrices (atomic_multipole.pkl) and atomic multipole SAMB of a point group in the jml/lgs/lg basis, keyed by PGMultipoleType (atomic_multipole_group.pkl).

```{code-cell} ipython3
import numpy as np

from multipiegen.data.data_basis import atomic_lgs_basis
from multipiegen.util.util_atomic_multipole import load_atomic_multipole, load_atomic_multipole_group
from multipiegen.util.util_harmonics import harmonics_expression
from multipiegen.util.util_info import load_info
```

```{code-cell} ipython3
# atomic multipole: {"lms"/"jml"/"lm": {(L1, L2): {SphericalMultipoleType(X, l, s, k, "q"): ndarray(2l+1, n1, n2)}}}.
am = load_atomic_multipole()
for key, mat in am["lms"][(1, 1)].items():
    print(key, mat.shape)
```

```{code-cell} ipython3
# atomic multipole SAMB of Oh: {"jml"/"lgs"/"lg": {(L1, L2): {PGMultipoleType: matrix (dim, n1, n2)}}}.
# the expression of the harmonics (X, l, Gamma, n) is looked up by the key in harmonics.pkl.
info = load_info()
pg = load_atomic_multipole_group()[info["id"]["Oh"]]
print("lgs basis of p:", atomic_lgs_basis[1])
for key, mat in list(pg["lgs"][(1, 1)].items())[:3]:
    print(key, "symmetry:", list(harmonics_expression("Oh", key)))
    print(mat[0])

# orthonormality: Re tr(A^dagger B) = delta (diagonal block) or delta/2 (off-diagonal block).
for block in [(1, 1), (1, 2)]:
    v = np.array([np.asarray(c, dtype=complex).ravel() for m in pg["lgs"][block].values() for c in m])
    g = (v.conj() @ v.T).real
    print(block, len(v), "matrices, orthonormal:", np.allclose(g, np.eye(len(v)) * (1 if block[0] == block[1] else 0.5)))
```

## harmonics_symbol.py

Polar and axial symbols of harmonics side by side for all 47 point groups (as in the harmonics PDF).

```{code-cell} ipython3
from multipiegen.util.util_harmonics_symbol import harmonics_symbol_table, polar_to_axial, multipole_symbol
```

```{code-cell} ipython3
# 1. table of Oh harmonics: (key, polar symbol, axial symbol, expression in LaTeX).
tbl = harmonics_symbol_table("Oh")
print("number of Oh harmonics (rank 0-11):", len(tbl))
for key, p_sym, a_sym, expr in tbl[:9]:
    print(f"    {str(key):22s} {p_sym:28s} {a_sym:30s} {expr}")
```

```{code-cell} ipython3
# 2. for Oh and D6h, the axial irrep is obtained by g <-> u (same component, no sign).
for key, p_sym, a_sym, expr in harmonics_symbol_table("D6h"):
    if key[0] == 3:
        print(f"    D6h: {p_sym:28s} -> {a_sym}")
```

```{code-cell} ipython3
# 3. polar -> axial for a single key (any of the 47 point groups).
a_key, sign = polar_to_axial("Oh", (2, "Eg", -1, 1))
print("\nOh  Q(Eg,2) ->", a_key, sign, multipole_symbol("G", a_key, sign))

#   Td: components are interchanged and a sign appears.
for key, p_sym, a_sym, expr in harmonics_symbol_table("Td")[4:9]:
    print(f"    Td: {p_sym:24s} -> {a_sym:26s} {expr}")

#   complex group: complex irreps E^(a), E^(b) and complex expressions.
for key, p_sym, a_sym, expr in harmonics_symbol_table("C3h-c")[1:4]:
    print(f"    C3h-c: {p_sym:30s} -> {a_sym:30s} {expr}")
```

```{code-cell} ipython3
# 4. create the PDFs of all 47 point groups (needs ptex2pdf):
#    python -m multipiegen.create_pdf.create_pdf_harmonics
```

## info.py

Overall information of groups (info.pkl): group id/tag/name, orbital basis and basis functions, root cluster.

```{code-cell} ipython3
from multipiegen.util.util_info import load_info
```

```{code-cell} ipython3
# 1. load info. str(info) is the description of the contents.
info = load_info()
print("keys:", list(info.keys()))
print(str(info)[:400], "...\n")
```

```{code-cell} ipython3
# 2. group id, tag and name.
print("PG:32 ->", info["tag"]["PG:32"], ", Oh ->", info["id"]["Oh"])
print("number of SG:", len(info["id_set"]["SG"]["all"]))
print("name of PG:32:", info["group_name"]["PG"]["PG:32"])
```

```{code-cell} ipython3
# 3. orbital basis (renamed tags) and basis functions.
h = info["harmonics"]
print("\nvariables: site", h["variable_site"], ", vector", h["variable_vector"])
print("\nregular atomic orbitals (lg):", h["atomic_basis"]["spinless"]["lg"])
print("cubic orbitals:", h["atomic_basis"]["spinless"]["cubic"][3])
print("hexagonal orbitals:", h["atomic_basis"]["spinless"]["hexagonal"][3])
for name in ["fc", "fx", "gc4"]:
    expr, u = h["basis_function"][name]
    print(f"  {name:4s}: key={h['harmonics_key'][name]}, expression={expr}, <m|gamma>={u}")
print("tesseral tag of fc:", h["tesseral"]["fc"], ", Wannier90 fxz2 ->", h["wannier90"]["fxz2"])
print("internal basis (cubic, s=3, polar):", h["internal_basis"]["cubic"][(3, "q")])
```

```{code-cell} ipython3
# 4. root cluster and compatibility relation.
rc = info["root_cluster"]
print("\nroot cluster:", rc["point_group"], ", scale:", rc["scale"])
print("  Oh  site0 (cart):", rc["oh_site0_cart"], ", vector0 (cart):", rc["oh_vector0_cart"])
print("  D6h site0 (cart):", rc["d6h_site0_cart"], ", site0 (hex. fractional):", rc["d6h_site0_frac"])
print("compatibility Oh -> Td (first 4):", info["character"]["compatibility_relation"][("Oh", "Td")][:4])
```

## binary_data.py

Binary data for MultiPie (binary_data): info.pkl and one file for each group, read with BinaryManager alone, and the look-up of the harmonics of a SAMB key in the PG file.

```{code-cell} ipython3
from multipie.util.util_binary import BinaryManager

from multipiegen.util.util_binary_data import binary_data_dir, group_file
```

```{code-cell} ipython3
# 1. info.pkl: information of all groups (the same as binary/info.pkl).
info = BinaryManager("info", topdir=binary_data_dir)
print("info:", list(info.keys()))
print("PG:32 ->", info["tag"]["PG:32"], ", name:", info["group_name"]["PG"]["PG:32"]["name"])
```

```{code-cell} ipython3
# 2. a PG file. str() shows the description of its keys.
pg = BinaryManager("PG/PG_32", topdir=binary_data_dir)
print("\nPG/PG_32:", list(pg.keys()))
print(str(pg)[:300], "...")
```

```{code-cell} ipython3
# 3. the SAMB hold no expression: look up the harmonics by (X, l, Gamma, n) of the key, with T -> Q and M -> G.
harmonics_type = {"Q": "Q", "G": "G", "T": "Q", "M": "G"}


def expression(harmonics, key):
    ex, _ = harmonics[(harmonics_type[key.X], key.l, key.Gamma, key.n, -1, 0, 0, "q")]
    return list(ex)


for key, samb in list(pg["cluster_samb"]["bond"]["3a@1o"].items())[:3]:
    print(" ", key, "ex =", expression(pg["harmonics"], key), "samb =", samb.tolist())
```

```{code-cell} ipython3
# 4. an SG file: its harmonics are those of its PG (info.PG).
sg = BinaryManager("SG/SG_221", topdir=binary_data_dir)
pg_of_sg = BinaryManager(group_file(sg["info"].PG), topdir=binary_data_dir)
print("\nSG/SG_221:", list(sg.keys()), ", PG:", sg["info"].PG)
w = list(sg["cluster_samb"]["site"].keys())[1]
for key, samb in list(sg["cluster_samb"]["site"][w].items())[:2]:
    print(" ", w, key, "ex =", expression(pg_of_sg["harmonics"], key))
```

```{code-cell} ipython3
# 5. an MPG file: active multipoles and response tensors.
mpg = BinaryManager("MPG/MPG_15.6.58", topdir=binary_data_dir)
print("\nMPG/MPG_15.6.58:", mpg["info"].international, list(mpg.keys()))
matrix, ex = mpg["response_tensor"][("T", 2, "a")]
print("  (T, 2, a):", matrix.tolist(), ex)
```
