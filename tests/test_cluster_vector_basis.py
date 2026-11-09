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
from multipiegen.util.util_cluster_vector_basis import (
    cluster_vector_basis_matrix,
    create_cluster_vector_basis,
    load_cluster_vector_basis,
)
from multipiegen.util.util_representation_matrix import _group, representation_matrix

data = load_cluster_vector_basis()
CASES = [(g, x) for g in ["Oh", "D6h"] for x in ["q", "g"]]


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_orthonormal(group, x):
    n = 3 * len(data[group]["cluster"])
    m, keys = cluster_vector_basis_matrix(data[group], x)
    assert m.shape == (n, n)
    assert np.allclose(m @ m.T, np.eye(n))


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_irrep_count(group, x):
    # vector (3 dim) x regular representation: each irrep appears 3 d^2 times (d components).
    count = collections.Counter(key[2] for key in data[group]["basis"][x])
    for irrep, chi in character[group][1].items():
        d = int(chi.strip("[]").split(",")[0])
        assert count[irrep] == 3 * d * d, irrep


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_type(group, x):
    # polar internal basis: Q for even s+k, G for odd s+k (s=1); interchanged for axial internal basis.
    for X, rank, irrep, mul, comp, k in data[group]["basis"][x]:
        even = (1 + k) % 2 == 0
        assert X == (("Q" if even else "G") if x == "q" else ("G" if even else "Q"))


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_representation_matrix(group, x):
    # each multiplet transforms with the representation matrices of its irrep:
    # (g v)(p) = s R v(R^{-1} p), s = det R for axial internal basis.
    P = np.array(data[group]["cluster"], dtype=float)
    _, R, det = _group(group)
    D = {ir: m.astype(complex).real for ir, m in representation_matrix(group)["matrix"].items()}
    perm = [[int(np.argmin(np.linalg.norm(P - q, axis=1))) for q in P @ Rg] for Rg in R]
    blocks = collections.defaultdict(dict)
    for (X, rank, irrep, mul, comp, k), v in data[group]["basis"][x].items():
        blocks[(X, rank, irrep, mul, k)][comp] = np.array(v, dtype=float)
    for (X, rank, irrep, mul, k), comps in blocks.items():
        V = [comps[c] for c in sorted(comps)]
        for g in range(len(R)):
            s = det[g] if x == "g" else 1.0
            for c in range(len(V)):
                gv = s * (V[c][perm[g]] @ R[g].T)
                pred = sum(V[cp] * D[irrep][g][cp, c] for cp in range(len(V)))
                assert np.allclose(gv, pred), (X, rank, irrep, mul, k, g)


# ==================================================
@pytest.mark.parametrize("group", ["Oh", "D6h"])
def test_binary_is_up_to_date(group):
    new = create_cluster_vector_basis(group)
    old = data[group]
    for x in ["q", "g"]:
        assert list(new["basis"][x].keys()) == list(old["basis"][x].keys())
        assert new["dropped"][x] == old["dropped"][x]
        for k in new["basis"][x]:
            assert repr(new["basis"][x][k]) == repr(old["basis"][x][k]), k
