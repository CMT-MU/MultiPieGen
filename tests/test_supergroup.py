# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest
from multipiegen.util.util_coset import load_coset
from multipiegen.util.util_info import load_info
from multipiegen.util.util_supergroup import _wkey, direct_supergroups, exact_matrix, group_data, load_supergroup

info = load_info()
data = load_supergroup()


def _parents(gid):
    return sorted({d["parent"] for d in data[gid]}, key=lambda p: int(p.split(":")[1]))


# ==================================================
def test_groups():
    assert list(data) == info["id_set"]["PG"]["all"] + info["id_set"]["SG"]["all"]
    # no supergroup: the top groups D6h, Oh, P6/mmm (Pm-3m is a k-subgroup of Im-3m and Fm-3m).
    assert [g for g, v in data.items() if not v] == ["PG:27", "PG:32", "SG:191"]


# ==================================================
@pytest.mark.parametrize("gid", ["PG:8", "PG:1", "SG:1", "SG:152", "SG:198", "SG:155"])
def test_binary_is_up_to_date(gid):
    new = [(d["parent"], exact_matrix(d["A"]), d["index"], d["index_parent"]) for d in direct_supergroups(gid)]
    old = [(d["parent"], d["A"], d["index"], d["index_parent"]) for d in data[gid]]
    assert new == old


# ==================================================
def test_examples():
    # P3_121: k (154 with 2c, 150 with c/3, 155 R32) and t (178, 181; 151 with a hexagonal supercell).
    assert _parents("SG:152") == ["SG:150", "SG:151", "SG:154", "SG:155", "SG:178", "SG:181"]
    assert _parents("SG:155") == ["SG:149", "SG:166", "SG:167"] + [f"SG:{n}" for n in range(207, 215)]
    assert _parents("SG:198") == ["SG:196", "SG:199", "SG:205", "SG:212", "SG:213"]
    assert _parents("PG:8") == ["PG:15", "PG:29"]
    d = next(d for d in data["SG:152"] if d["parent"] == "SG:155")
    assert (d["index"], d["index_parent"]) == (3, 1)
    d = next(d for d in data["SG:198"] if d["parent"] == "SG:205")
    assert (d["index"], d["index_parent"], d["A"].tolist()) == (2, 2, np.eye(4).astype(int).tolist())


# ==================================================
def test_embedding():
    # A maps every operation of H to one of P, R is right-handed.
    for h, lst in data.items():
        Hd = group_data(h)
        for d in lst:
            A = np.array(d["A"].evalf(), dtype=float)
            assert np.linalg.det(A[:3, :3]) > 0
            Pd = group_data(d["parent"])
            Ai = np.linalg.inv(A)
            for W in Hd["full"]:
                M = A @ W @ Ai
                t = Pd["point"].get(_wkey(M[:3, :3]))
                assert t is not None, (h, d["parent"])
                v = M[:3, 3] - t[1]
                if Pd["L"] is not None:
                    v = np.linalg.solve(Pd["L"], v)
                    v = v - np.rint(v)
                assert np.allclose(v, 0, atol=1e-6), (h, d["parent"])
            n = len(Pd["point"]) // len(Hd["point"])
            assert d["index_parent"] == n and d["index"] % n == 0


# ==================================================
def test_parent_chain_reachable():
    # every parent of the chain of representative parents (coset.pkl) is reached by the direct supergroups.
    for h, v in load_coset().items():
        reach, frontier = set(), [h]
        while frontier:
            new = []
            for q in frontier:
                for d in data[q]:
                    if d["parent"] not in reach:
                        reach.add(d["parent"])
                        new.append(d["parent"])
            frontier = new
        assert set(v["parent"]) <= reach, h
