# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: orthonormal bond cluster basis of the Oh and D6h clusters.

    python docs/examples/bond_cluster_basis.py

The binary data (multipiegen/binary/bond_cluster_basis.pkl) is created by
multipiegen/create_binary/create_bond_cluster_basis.py (python -m multipiegen.create_binary.create_bond_cluster_basis)
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
        print(
            f"  x={x}: {len(keys)} vectors, orthonormal: {ok}, dropped: {len(d['dropped'][x])}, max rank: {max(k[1] for k in keys)}"
        )

# ==================================================
# keys: (X, rank, irrep, multiplicity, component, k) of the vector cluster basis; value: (N,) values on the bonds.
d = data["Oh"]
for key in list(d["basis"]["q"])[:5]:
    print(key, d["basis"]["q"][key][:3], "...")
print("dropped (first 5):", d["dropped"]["q"][:5])
