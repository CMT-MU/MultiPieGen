# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: binary data for MultiPie (multipiegen/binary_data), one file for each group.

    python docs/examples/binary_data.py

The binary data are created by python -m multipiegen.create_binary.create_binary_data from the binaries in
multipiegen/binary. They are read with MultiPie's BinaryManager alone, without code from MultiPieGen.
"""

from multipie.util.util_binary import BinaryManager

from multipiegen.util.util_binary_data import binary_data_dir, group_file

# ==================================================
# 1. info.pkl: information of all groups (the same as binary/info.pkl).
info = BinaryManager("info", topdir=binary_data_dir)
print("info:", list(info.keys()))
print("PG:32 ->", info["tag"]["PG:32"], ", name:", info["group_name"]["PG"]["PG:32"]["name"])

# ==================================================
# 2. a PG file. str() shows the description of its keys.
pg = BinaryManager("PG/PG_32", topdir=binary_data_dir)
print("\nPG/PG_32:", list(pg.keys()))
print(str(pg)[:300], "...")

# ==================================================
# 3. the SAMB hold no expression: look up the harmonics by (X, l, Gamma, n) of the key, with T -> Q and M -> G.
harmonics_type = {"Q": "Q", "G": "G", "T": "Q", "M": "G"}


def expression(harmonics, key):
    ex, _ = harmonics[(harmonics_type[key.X], key.l, key.Gamma, key.n, -1, 0, 0, "q")]
    return list(ex)


for key, samb in list(pg["cluster_samb"]["bond"]["3a@1o"].items())[:3]:
    print(" ", key, "ex =", expression(pg["harmonics"], key), "samb =", samb.tolist())

# ==================================================
# 4. an SG file: its harmonics are those of its PG (info.PG).
sg = BinaryManager("SG/SG_221", topdir=binary_data_dir)
pg_of_sg = BinaryManager(group_file(sg["info"].PG), topdir=binary_data_dir)
print("\nSG/SG_221:", list(sg.keys()), ", PG:", sg["info"].PG)
w = list(sg["cluster_samb"]["site"].keys())[1]
for key, samb in list(sg["cluster_samb"]["site"][w].items())[:2]:
    print(" ", w, key, "ex =", expression(pg_of_sg["harmonics"], key))

# ==================================================
# 5. an MPG file: active multipoles and response tensors.
mpg = BinaryManager("MPG/MPG_15.6.58", topdir=binary_data_dir)
print("\nMPG/MPG_15.6.58:", mpg["info"].international, list(mpg.keys()))
matrix, ex = mpg["response_tensor"][("T", 2, "a")]
print("  (T, 2, a):", matrix.tolist(), ex)
