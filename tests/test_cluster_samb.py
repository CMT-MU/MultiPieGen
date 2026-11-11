# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pickle

import numpy as np
import pytest
from multipiegen.data.data_character_table import character
from multipiegen.util.util_op_mapping import op_mapping_d6h_oh_cluster
from multipiegen.util.util_cluster_samb import create_cluster_samb, load_cluster_samb
from multipiegen.util.util_group import load_group
from multipiegen.util.util_harmonics import harmonics_expression
from multipiegen.util.util_info import load_info

info = load_info()
group = load_group()
data = load_cluster_samb()
SAMPLE = ["PG:32", "PG:27", "PG:17", "PG:13", "SG:4", "SG:18", "SG:148", "SG:194", "SG:225"]


def _args(gid):
    pg_id = gid if gid.startswith("PG") else group[gid]["info"].PG
    return info["tag"][pg_id], group[pg_id]["symmetry_operation"]["tag"], group[gid]["symmetry_operation"]["tag"]


# ==================================================
def test_all_groups():
    assert list(data.keys()) == info["id_set"]["PG"]["all"] + info["id_set"]["SG"]["all"]


# ==================================================
def test_operation_mapping():
    # coset 0 of a point group lists its operations in the top group; a space group has the same point parts.
    for gid in info["id_set"]["PG"]["all"]:
        top, cosets = op_mapping_d6h_oh_cluster[info["tag"][gid]]
        top_ops = group[info["id"][top]]["symmetry_operation"]["tag"]
        assert [top_ops[i - 1] for i in cosets[0]] == group[gid]["symmetry_operation"]["tag"]
    for gid in info["id_set"]["SG"]["all"]:
        pg_ops = group[group[gid]["info"].PG]["symmetry_operation"]["tag"]
        assert sorted(t.split(":")[0] for t in group[gid]["symmetry_operation"]["tag"]) == sorted(pg_ops)


# ==================================================
@pytest.mark.parametrize("gid", ["PG:13", "PG:17", "SG:4", "SG:148"])
def test_binary_is_up_to_date(gid):
    new = pickle.loads(pickle.dumps(create_cluster_samb(*_args(gid), group[gid]["wyckoff"])))
    assert repr(new) == repr(data[gid])


# ==================================================
@pytest.mark.parametrize("gid", SAMPLE)
def test_orthonormal_and_complete(gid):
    pg = _args(gid)[0]
    irreps = list(character[pg][1])
    wyckoff = group[gid]["wyckoff"]
    assert list(data[gid].keys()) == ["site", "bond"]
    for kind in ["site", "bond"]:
        assert list(data[gid][kind].keys()) == list(wyckoff[kind].keys())
        for w, samb in data[gid][kind].items():
            n = len(wyckoff[kind][w]["mapping"])
            assert all(k.Gamma in irreps and k.s == 0 and k.k == 0 and k.x == "q" for k in samb)
            assert all(len(s) == len(harmonics_expression(pg, k)) for k, s in samb.items())
            # order: Gamma (character table), and Q before T/M in each Gamma.
            order = [(irreps.index(k.Gamma), k.X != "Q") for k in samb]
            assert order == sorted(order)
            # Q: real; T/M: pure imaginary (i times the anti-symmetric bond SAMB). Q (symmetric functions on the sites
            # or bonds) is complete. Q and T/M are orthonormal separately (different time-reversal parity).
            assert all(k.X in (["Q"] if kind == "site" else ["Q", "T", "M"]) for k in samb)
            for X, f in [("Q", 1), ("TM", -1j)]:
                part = [np.array(s, dtype=complex) * f for k, s in samb.items() if k.X in X]
                if not part:
                    continue
                v = np.vstack(part)
                assert np.allclose(v.imag, 0, atol=1e-12)
                v = v.real
                assert v.shape[1] == n
                assert np.allclose(v @ v.T, np.eye(len(v)), atol=1e-10)
                if X == "Q":
                    assert len(v) == n
