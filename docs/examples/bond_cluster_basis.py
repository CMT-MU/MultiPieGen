# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: orthonormal bond cluster basis of the Oh and D6h clusters.

    python docs/examples/bond_cluster_basis.py

The binary data (multipiegen/binary/cluster_bond_basis.pkl) is created by
multipiegen/create_binary/create_cluster_bond_basis.py (python -m multipiegen.create_binary.create_cluster_bond_basis)
from the vector cluster basis (cluster_vector_basis.pkl).
"""

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

# ==================================================
# keys: PGMultipoleType(X, l, Gamma, n, p, s=0, k=0, x) as MultiPie Ver.2 bond_a; value: (dim, N) values on the bonds.
# source: the multiplet of the vector cluster basis (PGMultipoleType with s=1).
d = data["Oh"]
for key in list(d["basis"])[:4]:
    print(key, d["basis"][key].shape, "<-", d["source"][key])
print("with add_multiplicity:", [k for k in d["basis"] if k.p != -1][:2])
print("dropped vector cluster basis multiplets (first 2):", d["dropped"][:2])
