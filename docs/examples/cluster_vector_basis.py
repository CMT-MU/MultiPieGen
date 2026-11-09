# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: orthonormal vector cluster basis (internal rank s=1) of the Oh and D6h clusters.

    python docs/examples/cluster_vector_basis.py

The binary data (multipiegen/binary/cluster_vector_basis.pkl) is created by
multipiegen/create_binary/create_cluster_vector_basis.py (python -m multipiegen.create_binary.create_cluster_vector_basis).
"""

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

# ==================================================
# keys: PGMultipoleType(X, l, Gamma, n, p=-1, s=1, k, x) as MultiPie vector SAMB; value: (dim, N, 3) at the cluster points.
d = data["Oh"]
for key in list(d["basis"])[:4]:
    print(key, "at the first point:", d["basis"][key][:, 0].tolist())

# a vector field on the cluster expanded in the basis: coefficients = <basis|field>.
P = np.array(d["cluster"], dtype=float)
field = np.cross(P, [0, 0, 1])  # rotation around z (axial dipole like).
m, keys = cluster_vector_basis_matrix(d, "q")
c = m @ field.ravel()
print("nonzero components of r x e_z:", [(keys[i], round(c[i], 6)) for i in range(len(c)) if abs(c[i]) > 1e-9])
