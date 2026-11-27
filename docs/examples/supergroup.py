# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: direct (maximal) supergroups of a group with the embeddings x_P = A x_H (supergroup.pkl), and the chains of
them up to the representative parents of coset.pkl.

    python docs/examples/supergroup.py

The binary data (multipiegen/binary/supergroup.pkl) is created by python -m multipiegen.create_binary.create_supergroup
from group.pkl and info.pkl. In binary_data, it is "supergroup" of info.pkl.
"""

from multipiegen.util.util_supergroup import load_supergroup

data = load_supergroup()


def show(gid):
    print(f"{gid}: {len(data[gid])} direct supergroups")
    for d in data[gid]:
        kind = "t" if d["index"] == d["index_parent"] else ("k" if d["index_parent"] == 1 else "t+k")
        print(
            f"  {d['parent']:7s} index {d['index']} (modulo the lattice of P: {d['index_parent']}, {kind}), A = {d['A'].tolist()}"
        )


# ==================================================
# P3_121 (Te, No.152): R-3m (No.166) is not a direct supergroup, but a direct supergroup of R32 (No.155).
show("SG:152")
show("SG:155")


# ==================================================
# the chain of direct supergroups from SG 152 to SG 166, and to the representative parent P-3m1 (No.164).
def chains(h, p, path=None):
    path = [h] if path is None else path
    if h == p:
        yield path
        return
    for d in data[h]:
        if d["parent"] not in path:
            yield from chains(d["parent"], p, path + [d["parent"]])


for p in ["SG:166", "SG:164"]:
    print(p, "<-", min(chains("SG:152", p), key=len))

# ==================================================
# point groups: D2h in D4h (4-fold axis along z) and in Th; no supergroup for the top groups.
show("PG:8")
print([g for g, v in data.items() if not v])
