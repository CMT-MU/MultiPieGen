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
def _multiplets(group, x):
    return {key: np.array(v, dtype=float) for key, v in data[group]["basis"].items() if key.x == x}


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_irrep_count(group, x):
    # vector (3 dim) x regular representation: each irrep appears 3 d^2 times (d components).
    count = collections.Counter()
    for key, v in _multiplets(group, x).items():
        count[key.Gamma] += len(v)
    for irrep, chi in character[group][1].items():
        d = int(chi.strip("[]").split(",")[0])
        assert count[irrep] == 3 * d * d, irrep


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_keys(group, x):
    # PGMultipoleType as MultiPie vector SAMB: p = -1, s = 1; Q for even s+k (q), interchanged for g.
    for key in _multiplets(group, x):
        assert key.p == -1 and key.s == 1
        even = (1 + key.k) % 2 == 0
        assert key.X == (("Q" if even else "G") if x == "q" else ("G" if even else "Q"))
    keys = [k for k in data[group]["basis"]] + data[group]["dropped"]
    assert len(set(keys)) == len(keys)


# ==================================================
@pytest.mark.parametrize("group, x", CASES)
def test_representation_matrix(group, x):
    # each multiplet transforms with the representation matrices of its irrep:
    # (g v)(p) = s R v(R^{-1} p), s = det R for axial internal basis.
    P = np.array(data[group]["cluster"], dtype=float)
    _, R, det = _group(group)
    D = {ir: m.astype(complex).real for ir, m in representation_matrix(group)["matrix"].items()}
    perm = [[int(np.argmin(np.linalg.norm(P - q, axis=1))) for q in P @ Rg] for Rg in R]
    for key, V in _multiplets(group, x).items():
        assert len(V) == D[key.Gamma].shape[1], key
        for g in range(len(R)):
            s = det[g] if x == "g" else 1.0
            for c in range(len(V)):
                gv = s * (V[c][perm[g]] @ R[g].T)
                pred = sum(V[cp] * D[key.Gamma][g][cp, c] for cp in range(len(V)))
                assert np.allclose(gv, pred), (key, g)


# ==================================================
@pytest.mark.parametrize("group", ["Oh", "D6h"])
def test_binary_is_up_to_date(group):
    new = create_cluster_vector_basis(group)
    old = data[group]
    assert list(new["basis"].keys()) == list(old["basis"].keys())
    assert new["dropped"] == old["dropped"]
    for k in new["basis"]:
        assert repr(new["basis"][k]) == repr(old["basis"][k]), k
