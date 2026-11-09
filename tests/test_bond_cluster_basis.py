# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import collections

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


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_orthonormal(group, x):
    n = len(data[group]["cluster"])
    m, keys = bond_cluster_basis_matrix(data[group], x)
    assert m.shape == (n, n)
    assert np.allclose(m @ m.T, np.eye(n))


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_irrep_count(group, x):
    # scalar functions on N points (regular representation): each irrep appears d^2 times.
    count = collections.Counter(key[2] for key in data[group]["basis"][x])
    for irrep, chi in character[group][1].items():
        d = int(chi.strip("[]").split(",")[0])
        assert count[irrep] == d * d, irrep


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_keys_from_vector_basis(group, x):
    # keys are those of the vector cluster basis, in the same order; kept + dropped = all.
    vkeys = list(load_cluster_vector_basis()[group]["basis"][x].keys())
    keys = list(data[group]["basis"][x].keys())
    assert keys == [k for k in vkeys if k in keys]
    assert sorted(keys + data[group]["dropped"][x]) == sorted(vkeys)


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
    blocks = collections.defaultdict(dict)
    for (X, rank, irrep, mul, comp, k), v in data[group]["basis"][x].items():
        blocks[(X, rank, irrep, mul, k)][comp] = np.array(v, dtype=float)
    for (X, rank, irrep, mul, k), comps in blocks.items():
        assert len(comps) == D[irrep].shape[1], (X, rank, irrep, mul, k)  # complete multiplet.
        V = [comps[c] for c in sorted(comps)]
        for g in range(len(R)):
            s = det[g] if x == "g" else 1.0
            for c in range(len(V)):
                pred = sum(V[cp] * D[irrep][g][cp, c] for cp in range(len(V)))
                assert np.allclose(s * V[c][perm[g]], pred), (X, rank, irrep, mul, k, g)


# ==================================================
@pytest.mark.parametrize("group", ["Oh", "D6h"])
def test_binary_is_up_to_date(group):
    new = create_bond_cluster_basis(group)
    old = data[group]
    for x in ["q", "g"]:
        assert list(new["basis"][x].keys()) == list(old["basis"][x].keys())
        assert new["dropped"][x] == old["dropped"][x]
        for k in new["basis"][x]:
            assert repr(new["basis"][x][k]) == repr(old["basis"][x][k]), k
