# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: cluster SAMB of the Wyckoff sites and bonds of a point group and a space group.

    python docs/examples/cluster_samb.py

The binary data (multipiegen/binary/cluster_samb.pkl) is created by
python -m multipiegen.create_binary.create_cluster_samb from group.pkl, cluster_scalar_basis.pkl, and
cluster_bond_basis.pkl.
"""

import numpy as np

from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_info import load_info

info = load_info()
data = load_cluster_samb()

# ==================================================
# {PG_id/SG_id: {"site"/"bond_s"/"bond_a": {wyckoff: {PGMultipoleType: (samb (dim, N), ex (dim))}}}}.
d = data[info["id"]["D4h"]]
for kind in ["site", "bond_s", "bond_a"]:
    w = list(d[kind].keys())[1]
    print(f"* {kind} {w}")
    for key, (samb, ex) in d[kind][w].items():
        print("  ", key, "ex =", list(ex))
        print("  ", samb.tolist())

# orthonormality on the N sites (bonds) of each Wyckoff position.
d = data["SG:225"]
for kind in ["site", "bond_s", "bond_a"]:
    ok = all(
        np.allclose((v := np.vstack([np.array(s, dtype=float) for s, _ in samb.values()])) @ v.T, np.eye(len(v)))
        for samb in d[kind].values()
    )
    print(f"SG:225 {kind}: {len(d[kind])} Wyckoff positions, orthonormal: {ok}")
