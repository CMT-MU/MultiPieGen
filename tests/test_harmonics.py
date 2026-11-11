# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pickle

import numpy as np
import pytest
import sympy as sp
from multipiegen.util.util_harmonics_mapping import data_harmonics_mapping
from multipiegen.util.util_atomic_multipole import load_atomic_multipole_group
from multipiegen.util.util_bond_cluster_basis import load_bond_cluster_basis
from multipiegen.util.util_cluster_basis import load_cluster_basis
from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_cluster_vector_basis import load_cluster_vector_basis
from multipiegen.util.util_group import load_group
from multipiegen.util.util_harmonics import (
    create_harmonics,
    harmonics_expression,
    harmonics_type,
    harmonics_u_matrix,
    load_harmonics,
    max_rank,
)
from multipiegen.util.util_harmonics_multipole import load_multipolar_harmonics
from multipiegen.util.util_info import load_info

info = load_info()
data = load_harmonics()
SAMPLE = ["Oh", "D6h", "C3h-c", "D2d-1", "C1"]


# ==================================================
def test_all_groups():
    assert list(data.keys()) == info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]


# ==================================================
@pytest.mark.parametrize("pg", SAMPLE)
def test_binary_is_up_to_date(pg):
    new = pickle.loads(pickle.dumps(create_harmonics(pg)))
    assert repr(new) == repr(data[info["id"][pg]])


# ==================================================
@pytest.mark.parametrize("pg", SAMPLE)
def test_keys_and_mapping(pg):
    d = data[info["id"][pg]]
    m = data_harmonics_mapping[pg]
    keys = {(X, l, irrep, n) for X in ["Q", "G"] for (l, irrep, n, _) in m[X] if l <= max_rank}
    assert {(k.X, k.l, k.Gamma, k.n) for k in d} == keys
    for k, (ex, U) in d.items():
        assert (k.p, k.s, k.k, k.x) == (-1, 0, 0, "q")
        assert U.shape == (2 * k.l + 1, len(ex))
        assert np.allclose(np.array(U, dtype=complex).conj().T @ np.array(U, dtype=complex), np.eye(len(ex)))  # orthonormal


# ==================================================
def test_lookup_any_key():
    # a key with any (p, s, k, x), and the PG given by tag or id, gives the harmonics of its (X, l, Gamma, n).
    k = list(load_multipolar_harmonics()[info["id"]["Oh"]])[5]
    ex, U = data["PG:32"][(k.X, k.l, k.Gamma, k.n, -1, 0, 0, "q")]
    assert repr(harmonics_expression("Oh", k)) == repr(ex) == repr(harmonics_expression("PG:32", tuple(k)))
    assert repr(harmonics_u_matrix("Oh", k)) == repr(np.asarray(U))


# ==================================================
def test_keys_of_other_binaries_are_found():
    group = load_group()
    mh = load_multipolar_harmonics()
    amg = load_atomic_multipole_group()
    pgs = info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]
    for gid in pgs:
        d = data[gid]
        find = lambda k: (harmonics_type[k.X], k.l, k.Gamma, k.n, -1, 0, 0, "q") in d
        assert all(find(k) for k in mh[gid])
        assert all(find(k) for blocks in amg[gid].values() for b in blocks.values() for k in b)
    samb = load_cluster_samb()
    for gid in samb.keys():
        pg_id = gid if gid.startswith("PG") else group[gid]["info"].PG
        d = data[pg_id]
        for kind in samb[gid].values():
            for w in kind.values():
                for k, s in w.items():
                    key = (harmonics_type[k.X], k.l, k.Gamma, k.n, -1, 0, 0, "q")
                    assert key in d and len(d[key][0]) == len(s)
    bases = [load_cluster_basis(), load_cluster_vector_basis(), load_bond_cluster_basis()]
    for top, gid in [("Oh", "PG:32"), ("D6h", "PG:27")]:
        d = data[gid]
        for basis in bases:
            assert all((k.X, k.l, k.Gamma, k.n, -1, 0, 0, "q") in d for k in basis[top]["basis"])
