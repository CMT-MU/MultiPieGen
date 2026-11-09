# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: overall information of groups (info.pkl).

    python docs/examples/info.py

The binary data (multipiegen/binary/info.pkl) is created by
multipiegen/create_binary/create_info.py (python -m multipiegen.create_binary.create_info).
"""

from multipiegen.util.util_info import load_info

# ==================================================
# 1. load info. str(info) is the description of the contents.
info = load_info()
print("keys:", list(info.keys()))
print(str(info)[:400], "...\n")

# ==================================================
# 2. group id, tag and name.
print("PG:32 ->", info["tag"]["PG:32"], ", Oh ->", info["id"]["Oh"])
print("number of SG:", len(info["id_set"]["SG"]["all"]))
print("name of PG:32:", info["group_name"]["PG"]["PG:32"])

# ==================================================
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

# ==================================================
# 4. root cluster and compatibility relation.
rc = info["root_cluster"]
print("\nroot cluster:", rc["point_group"], ", scale:", rc["scale"])
print("  Oh  site0 (cart):", rc["oh_site0_cart"], ", vector0 (cart):", rc["oh_vector0_cart"])
print("  D6h site0 (cart):", rc["d6h_site0_cart"], ", site0 (hex. fractional):", rc["d6h_site0_frac"])
print("compatibility Oh -> Td (first 4):", info["character"]["compatibility_relation"][("Oh", "Td")][:4])
