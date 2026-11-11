# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import collections
import functools

import numpy as np
import pytest
from multipiegen.data.data_character_table import character
from multipiegen.util.util_bond_cluster_basis import (
    bond_cluster_basis_matrix,
    create_bond_cluster_basis,
    load_bond_cluster_basis,
)
from multipiegen.util.util_cluster_vector_basis import load_cluster_vector_basis
from multipiegen.util.util_representation_matrix import _group, representation_matrix

data = load_bond_cluster_basis()
CASES = [(g, x) for g in ["Oh", "D6h"] for x in ["q", "g"]]


@functools.cache
def _vector_basis():
    return load_cluster_vector_basis()


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_orthonormal(group, x):
    n = len(data[group]["cluster"])
    m, keys = bond_cluster_basis_matrix(data[group], x)
    assert m.shape == (n, n)
    assert np.allclose(m @ m.T, np.eye(n))


# ==================================================
def _multiplets(group, x):
    return {key: np.array(v, dtype=float) for key, v in data[group]["basis"].items() if key.x == x}


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_irrep_count(group, x):
    # scalar functions on N points (regular representation): each irrep appears d^2 times.
    count = collections.Counter()
    for key, v in _multiplets(group, x).items():
        count[key.Gamma] += len(v)
    for irrep, chi in character[group][1].items():
        d = int(chi.strip("[]").split(",")[0])
        assert count[irrep] == d * d, irrep


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_keys(group, x):
    # keys as MultiPie Ver.2 bond_a: PGMultipoleType(X, l, Gamma, n, p, s=0, k=0, x);
    # source + dropped = all multiplets of the vector cluster basis.
    vkeys = [k for k in _vector_basis()[group]["basis"] if k.x == x]
    src = []
    for key in _multiplets(group, x):
        assert key.s == 0 and key.k == 0
        v = data[group]["source"][key]
        assert (v.X, v.l, v.Gamma, v.n, v.x) == (key.X, key.l, key.Gamma, key.n, x)
        src.append(v)
    dropped = [k for k in data[group]["dropped"] if k.x == x]
    assert sorted(src + dropped) == sorted(vkeys)
    # add_multiplicity: -1 if unique, otherwise 1, 2, ... in the order of k of the source.
    same = collections.defaultdict(list)
    for key in _multiplets(group, x):
        same[(key.X, key.l, key.Gamma, key.n)].append(key)
    for lst in same.values():
        if len(lst) == 1:
            assert lst[0].p == -1
        else:
            lst = sorted(lst, key=lambda k: data[group]["source"][k].k)
            assert [k.p for k in lst] == list(range(1, len(lst) + 1))


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_representation_matrix(group, x):
    # (g f)(i) = s f(j), R_g^{-1} p_i = p_j, s = det R for axial internal basis, with the matrices of its irrep.
    P = np.array(data[group]["cluster"], dtype=float)
    B = np.array(data[group]["bond"], dtype=float)
    _, R, det = _group(group)
    D = {ir: m.astype(complex).real for ir, m in representation_matrix(group)["matrix"].items()}
    perm = [[int(np.argmin(np.linalg.norm(P - q, axis=1))) for q in P @ Rg] for Rg in R]
    for Rg, pp in zip(R, perm):  # bonds move with points.
        assert np.allclose(B @ Rg, B[pp])
    for key, V in _multiplets(group, x).items():
        assert len(V) == D[key.Gamma].shape[1], key  # complete multiplet.
        for g in range(len(R)):
            s = det[g] if x == "g" else 1.0
            for c in range(len(V)):
                pred = sum(V[cp] * D[key.Gamma][g][cp, c] for cp in range(len(V)))
                assert np.allclose(s * V[c][perm[g]], pred), (key, g)


# ==================================================
@pytest.mark.parametrize("group", ["Oh", "D6h"])
def test_binary_is_up_to_date(group):
    new = create_bond_cluster_basis(group)
    old = data[group]
    assert list(new["basis"].keys()) == list(old["basis"].keys())
    assert new["dropped"] == old["dropped"]
    assert new["source"] == old["source"]
    for k in new["basis"]:
        assert repr(new["basis"][k]) == repr(old["basis"][k]), k
