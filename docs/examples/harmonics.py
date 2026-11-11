# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: expression and u-matrix of the harmonics of a point group, looked up by a multipole key.

    python docs/examples/harmonics.py

The binary data (multipiegen/binary/harmonics.pkl) is created by python -m multipiegen.create_binary.create_harmonics.
The other binaries keyed by PGMultipoleType do not hold the expression of the harmonics; look it up here.
"""

from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_harmonics import harmonics_expression, harmonics_u_matrix, load_harmonics
from multipiegen.util.util_info import load_info

info = load_info()

# ==================================================
# 1. {PG_id: {PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q"): (expression (dim), u_matrix (2l+1, dim))}}.
d = load_harmonics()[info["id"]["D4h"]]
print("D4h:", len(d), "harmonics (rank <= 11)")
for key, (ex, U) in list(d.items())[:4]:
    print(" ", key, list(ex), U.shape)

# ==================================================
# 2. look up the harmonics of a key of another binary (only (X, l, Gamma, n) of the key is used).
samb = load_cluster_samb()[info["id"]["D4h"]]["site"]["4b"]
for key in list(samb)[:3]:
    print(key, "->", list(harmonics_expression("D4h", key)))
key = list(samb)[1]
print("u-matrix of", key, ":\n", harmonics_u_matrix("PG:15", key))
