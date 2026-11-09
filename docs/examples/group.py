# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: data of each group (group.pkl).

    python docs/examples/group.py

The binary data (multipiegen/binary_data/group.pkl) is created by
multipiegen/create_binary/create_group.py (python -m multipiegen.create_binary.create_group).
"""

from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

info = load_info()
group = load_group()
print("number of groups:", len(group.keys()), "(PG, SG, MPG, MSG)\n")

# ==================================================
# 1. symmetry operations of a point group (D4h).
d4h = group[info["id"]["D4h"]]
so = d4h["symmetry_operation"]
print("D4h:", d4h["info"].international, ", SO:", so["tag"])
print("  matrix of 4+[001] (cartesian):", so["cartesian"][so["tag"].index("4+[001]")].tolist())
print("  product 4+[001] * 2[100] =", so["product"][("4+[001]", "2[100]")])
print("  dipole matrix of 4+[001], axial (G, s=1):", so["so_matrix"][("G", 1)][2].tolist())

# ==================================================
# 2. mapping to the Oh cluster: mapping[0] is D4h itself, mapping[k] is a left coset in Oh.
top = so["mapping"]["top"]
oh_tag = group[info["id"][top]]["symmetry_operation"]["tag"]
print(f"\nmapping to {top}:")
for k, row in enumerate(so["mapping"]["mapping"]):
    print(f"  coset {k}:", [oh_tag[i - 1] for i in row][:6], "...")

# ==================================================
# 3. character table (complex PG has its own table, other data are those of the real PG).
c3h = group[info["id"]["C3h-c"]]
print("\nC3h-c irreps:", list(c3h["character"]["table"]))
print("  character of Ea' for all SO:", c3h["character"]["table_full"]["Ea'"].tolist())

# ==================================================
# 4. Wyckoff site and bond of a space group (P6_3/mmc).
sg = group["SG:194"]
site = sg["wyckoff"]["site"]["2a"]
print("\nSG:194", sg["info"].international, "Wyckoff 2a:", site["conventional"].tolist(), ", symmetry:", site["symmetry"])
print("  mapping:", site["mapping"])
print("  bonds at 2a:", site["bond"][:4])
bw = site["bond"][0]
bond = sg["wyckoff"]["bond"][bw]
print(f"  bond {bw}: symmetry {bond['symmetry']}, (vector, center) =", bond["conventional"].tolist())

# ==================================================
# 5. magnetic point group: active multipoles.
mpg = group["MPG:24.5.149"]
print("\nMPG", mpg["info"].international, "time-reversal sign:", mpg["symmetry_operation"]["tr_sign"].tolist())
print("  active multipoles:", mpg["active_multipole"])
