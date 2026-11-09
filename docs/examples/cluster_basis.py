# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: orthonormal cluster basis of Oh and D6h.

    python docs/examples/cluster_basis.py

The binary data (multipiegen/binary/cluster_scalar_basis.pkl) is created by
multipiegen/create_binary/create_cluster_scalar_basis.py
(python -m multipiegen.create_binary.create_cluster_scalar_basis).
"""

import numpy as np
from multipiegen.util.util_cluster_basis import load_cluster_basis, cluster_basis_matrix

# ==================================================
# 1. load the cluster basis of Oh and D6h.
data = load_cluster_basis()
print("groups:", list(data.keys()))

oh = data["Oh"]
print("site0 :", oh["site0"])
print("number of cluster points:", len(oh["cluster"]))
print("first 3 cluster points (exact):")
for p in oh["cluster"][:3]:
    print("   ", p)

# ==================================================
# 2. basis multiplets are keyed by PGMultipoleType(X="Q", l, Gamma, n, p=-1, s=0, k=0, x="q") (MultiPie site SAMB);
#    the value has all components, (dim, N).
print("\nbasis keys of Oh:")
for key in oh["basis"].keys():
    print("   ", key, oh["basis"][key].shape)

key = [k for k in oh["basis"] if k.l == 1 and k.Gamma == "T1u"][0]
print(f"\n{key}, x component, first 4 values (exact):", oh["basis"][key][0][:4])

# harmonics dropped by Gram-Schmidt (linearly dependent on lower ranks).
print("number of dropped harmonics:", len(oh["dropped"]), ", e.g.", oh["dropped"][:3])

# ==================================================
# 3. basis as a float matrix (each row is a component, columns follow oh["cluster"]); keys are (PGMultipoleType, comp).
m, keys = cluster_basis_matrix(oh)
print("\nbasis matrix shape:", m.shape)
print("orthonormal ?", np.allclose(m @ m.T, np.eye(len(keys))))

# ==================================================
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

# ==================================================
# 5. the same for D6h (24 points).
d6h = data["D6h"]
m6, keys6 = cluster_basis_matrix(d6h)
print("\nD6h basis matrix shape:", m6.shape, ", orthonormal ?", np.allclose(m6 @ m6.T, np.eye(len(keys6))))
print("D6h irreps:", sorted({k[0].Gamma for k in keys6}))
