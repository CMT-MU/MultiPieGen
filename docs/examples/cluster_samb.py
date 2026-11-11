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
from multipiegen.util.util_harmonics import harmonics_expression
from multipiegen.util.util_info import load_info

info = load_info()
data = load_cluster_samb()

# ==================================================
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
