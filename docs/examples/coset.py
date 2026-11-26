# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: cosets of a group in its parents (coset.pkl): the first coset in the direct parent, the image of the inversion
coset with the correspondence of the cluster SAMB, the parent structures, and the cosets of all parents.

    python docs/examples/coset.py

The binary data (multipiegen/binary/coset.pkl) is created by python -m multipiegen.create_binary.create_coset from
group.pkl and cluster_samb.pkl. In binary_data, "first" is "coset" of each group file, and "parent" and "coset" are in
info.pkl.
"""

from multipiegen.util.util_coset import load_coset
from multipiegen.util.util_info import load_info

info = load_info()
coset = load_coset()


def show(gid):
    d = coset[gid]
    f = d["first"]
    print(f"{gid} ({info['tag'][gid] if gid.startswith('PG') else gid}): parents {d['parent']}")
    print(f"  direct parent {f['parent']}, index {f['index']}, representative g = {f['rep']}, order {f['order']}, {f['kind']}")
    print(f"  g in the coordinates of H: {f['rep_matrix'].tolist()}")
    if f["kind"] == "rotation":
        print(f"  frame of the image system: {f['frame'].tolist()}")
    else:
        inv = f["inversion"]
        print(f"  image {inv['image']}, map {inv['map'].tolist()}")
        for kind in ["site", "bond"]:
            for w, t in list(inv[kind].items())[:3]:
                if t["diag"]:
                    eps = {f"{k.X}{k.l}({k.Gamma})": e for k, e in t["eps"].items()}
                    print(f"  {kind} {w} -> {t['image'][1]}, perm {t['perm']}, eps {eps}")
                else:
                    print(f"  {kind} {w} -> {t['image'][1]}, U in blocks {list(t['U'])}")
    for w, lims in f["parent_limit"].items():
        print(f"  parent structure of {w}: " + ", ".join(f"{l['param']} -> {l['parent_wyckoff']}" for l in lims))
    for p, cos in d["coset"].items():
        print(f"  {len(cos)} cosets in {p}, first: {cos[0][:4]} ...")
    print()


# ==================================================
# Te (P3_121, No.152): the inversion coset of P-3m1 maps it onto P3_221 (No.154).
show("SG:152")
# D2h in D4h: g = 4+[001] (rotation); the image system is written in the frame rotated by g.
show(info["id"]["D2h"])
# SG 210 (F4_132): Wyckoff positions mapped onto themselves with the transport matrix U (first non-diagonal block).
import numpy as np

inv = coset["SG:210"]["first"]["inversion"]
for kind in ["site", "bond"]:
    for w, t in inv[kind].items():
        if t["diag"]:
            continue
        for blk, (rows, rows2, U) in t["U"].items():
            if not np.allclose(np.abs(np.array(U, dtype=complex)), np.eye(len(rows))):
                print(f"SG:210 {kind} {w} -> {t['image'][1]}, block {blk}, {len(rows)} rows: U =\n{U}")
                break
        break
