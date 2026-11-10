# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: multipolar harmonics (internal rank s = 1, 2, 3, up to rank 3) with polar (q) and axial (g) internal basis.

    python docs/examples/harmonics_multipole.py

The binary data (multipiegen/binary/harmonics_multipole.pkl) is created by
multipiegen/create_binary/create_harmonics_multipole.py (python -m multipiegen.create_binary.create_harmonics_multipole).
"""

import collections

from multipiegen.util.util_harmonics_multipole import load_harmonics_multipole
from multipiegen.util.util_info import load_info

info = load_info()
hm = load_harmonics_multipole()

# ==================================================
# 1. keys: PGMultipoleType(X, l, Gamma, n, p, s, k, x); value: (cartesian ex., coefficient of internal basis, multipole ex.)
for pg in ["Oh", "D6h"]:
    d = hm[info["id"][pg]]
    print(pg, len(d), "multipolar harmonics:", dict(collections.Counter((k.x, k.s) for k in d)))
    crystal = "cubic" if pg == "Oh" else "hexagonal"
    print("  internal basis (s=3, q):", info["harmonics"]["internal_basis"][crystal][(3, "q")])

# ==================================================
# 2. an example: a dipole-type (s=1) multipole of rank 2 in Oh.
d = hm[info["id"]["Oh"]]
key = [k for k in d if k.s == 1 and k.l == 2 and k.Gamma == "Eg" and k.k == 1][0]
ex, mharm, m_ex = d[key]
print("\n", key)
print("  harmonics:", ex.tolist())
print("  coefficient of [Q_x, Q_y, Q_z] (first component):", mharm[0].tolist())
print("  multipole expression (first component):", m_ex[0])

# ==================================================
# 3. axial internal basis (g) of rank s=2 in D6h.
d = hm[info["id"]["D6h"]]
key = [k for k in d if k.x == "g" and k.s == 2][0]
print("\n", key, "\n  ", d[key][2].tolist())
