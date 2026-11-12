# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os
import pickle

import pytest
from multipiegen.util.util_binary_data import (
    binary_data_dir,
    create_group_data,
    create_info_data,
    group_ids,
    group_types,
    load_binary_data,
)
from multipiegen.util.util_harmonics import harmonics_type
from multipiegen.util.util_info import load_info

info = load_info()
KEYS = {
    "PG": [
        "info",
        "symmetry_operation",
        "wyckoff",
        "character",
        "representation_matrix",
        "harmonics",
        "multipolar_harmonics",
        "atomic_samb",
        "cluster_samb",
    ],
    "SG": ["info", "symmetry_operation", "wyckoff", "cluster_samb"],
    "MPG": ["info", "symmetry_operation", "wyckoff", "active_multipole", "response_tensor"],
    "MSG": ["info", "symmetry_operation", "wyckoff"],
}
SAMPLE = ["PG:32", "PG:27", "PG:15", "PG:41", "SG:221", "SG:194", "MPG:15.6.58", "MSG:221.92"]


def _round_trip(obj):
    # sympy re-evaluates some products when a pickle is loaded; compare after a round trip.
    return pickle.loads(pickle.dumps(obj))


# ==================================================
def test_files():
    files = {f"{tp}/{f[:-4]}" for tp in group_types for f in os.listdir(os.path.join(binary_data_dir, tp))}
    assert files == {f"{gid.split(':')[0]}/{gid}" for gid in group_ids()}
    assert len(group_ids()) == len(info["tag"]) == 47 + 230 + 157 + 1651


# ==================================================
def test_info_is_up_to_date():
    saved = load_binary_data("info")
    new = create_info_data()
    assert str(saved) == str(new)
    assert repr(saved) == repr(_round_trip(new))
    assert "group_name" in saved


# ==================================================
@pytest.mark.parametrize("gid", SAMPLE)
def test_group_is_up_to_date(gid):
    saved = load_binary_data(gid)
    new = create_group_data(gid)
    assert str(saved) == str(new)
    assert repr(saved) == repr(_round_trip(new))


# ==================================================
@pytest.mark.parametrize("tp", group_types)
def test_keys(tp):
    for gid in [g for g in group_ids() if g.startswith(tp + ":")][:: 10 if tp == "MSG" else 1]:
        keys = KEYS[tp] if not (tp == "PG" and info["tag"][gid].endswith("-c")) else KEYS[tp][:-1]
        assert list(load_binary_data(gid).keys()) == keys, gid


# ==================================================
@pytest.mark.parametrize("gid", ["PG:32", "PG:15", "SG:221", "SG:194"])
def test_harmonics_of_samb_keys_in_pg_file(gid):
    d = load_binary_data(gid)
    pg = load_binary_data(d["info"].PG) if gid.startswith("SG") else d
    harmonics = pg["harmonics"]
    find = lambda k: (harmonics_type[k.X], k.l, k.Gamma, k.n, -1, 0, 0, "q") in harmonics
    for kind in d["cluster_samb"].values():
        for samb in kind.values():
            assert all(find(k) for k in samb)
    if gid.startswith("PG"):
        assert all(find(k) for k in d["multipolar_harmonics"])
        assert all(find(k) for b in d["atomic_samb"].values() for blk in b.values() for k in blk)
