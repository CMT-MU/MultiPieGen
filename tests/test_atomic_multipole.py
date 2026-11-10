# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pickle

import numpy as np
import pytest
from multipiegen.data.data_basis import atomic_jm_basis, atomic_lms_basis, atomic_rank_block
from multipiegen.data.data_character_table import character
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.util.util_atomic_multipole import (
    atomic_multipole_matrix,
    create_atomic_multipole_pg,
    create_top_samb,
    heads,
    load_atomic_multipole,
    load_atomic_multipole_group,
    max_rank,
)
from multipiegen.util.util_harmonics_multipole import pg_harmonics
from multipiegen.util.util_info import load_info

info = load_info()
am = load_atomic_multipole()
amg = load_atomic_multipole_group()
SAMPLE = ["Oh", "D6h", "D3d", "C3h-c", "T-c", "D2d-1"]


_cache = {}


def _c(v):
    if v not in _cache:
        _cache[v] = complex(v)
    return _cache[v]


def _num(m):
    return np.array([_c(v) for v in m.ravel()], dtype=complex).reshape(m.shape)


# ==================================================
@pytest.mark.parametrize("arg", [("lms", "Q", 2, 0, 0), ("lms", "M", 3, 1, -1), ("jml", "G", 3, 1, 1), ("jml", "T", 4, 1, 0)])
def test_atomic_multipole_is_up_to_date(arg):
    bt, X, l, s, k = arg
    new = atomic_multipole_matrix(*arg)
    old = {b: d[key] for b, d in am[bt].items() for key in d if (key.X, key.l, key.s, key.k) == (X, l, s, k)}
    assert repr(pickle.loads(pickle.dumps(new))) == repr(old)


# ==================================================
def test_atomic_multipole_shape():
    for bt, basis in [("lms", atomic_lms_basis), ("jml", atomic_jm_basis), ("lm", atomic_lms_basis)]:
        assert list(am[bt].keys()) == atomic_rank_block
        for (L1, L2), d in am[bt].items():
            f = 2 if bt == "lm" else 1
            for key, mat in d.items():
                assert mat.shape == (2 * key.l + 1, len(basis[L1]) // f, len(basis[L2]) // f)
                assert not (mat == 0).all()


# ==================================================
@pytest.mark.parametrize("pg", SAMPLE)
def test_atomic_multipole_group_is_up_to_date(pg):
    top = data_harmonics_mapping[pg]["top"]
    ex = {key: e for key, (e, _) in pg_harmonics(pg, max_rank).items()}
    blocks = [(0, 1), (1, 1)]
    top_data = {bt: {b: create_top_samb(am, top, b, bt) for b in blocks} for bt in ["jml", "lgs"]}
    for bt in top_data:
        for b in atomic_rank_block:
            top_data[bt].setdefault(b, {})
    new = pickle.loads(pickle.dumps(create_atomic_multipole_pg(pg, top_data, ex)))
    old = amg[info["id"][pg]]
    for bt in ["jml", "lgs", "lg"]:
        for b in blocks:
            assert repr(new[bt][b]) == repr(old[bt][b])


# ==================================================
@pytest.mark.parametrize("pg", SAMPLE)
def test_atomic_multipole_group_keys_and_orthonormality(pg):
    irreps = list(character[pg][1])
    d = amg[info["id"][pg]]
    assert list(d.keys()) == ["jml", "lgs", "lg"]
    for bt, blocks in d.items():
        for (L1, L2), samb in blocks.items():
            keys = list(samb.keys())
            order = [(heads.index(k.X), k.s, k.k, k.l, irreps.index(k.Gamma), k.n) for k in keys]
            assert order == sorted(order) and len(set(order)) == len(order)
            assert all(k.p == -1 and k.x == "q" and k.l <= max_rank for k in keys)
            if not samb:
                continue
            v = np.array([c.ravel() for m, ex in samb.values() for c in _num(m)])
            g = v.conj() @ v.T
            # Q-T and G-M of the same irrep overlap only in the imaginary part in an off-diagonal block.
            assert np.allclose(g.real, np.eye(len(v)) * (1 if L1 == L2 else 0.5), atol=1e-10)
            if L1 == L2:
                assert np.allclose(g, np.eye(len(v)), atol=1e-10)
