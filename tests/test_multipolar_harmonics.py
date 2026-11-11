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
from multipiegen.data.data_character_table import character
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.util.util_harmonics_multipole import (
    _crystal,
    _orbital_basis,
    create_multipolar_harmonics,
    internal_symbols,
    internal_type,
    load_multipolar_harmonics,
    max_rank,
    max_s_rank,
)
from multipiegen.util.util_harmonics import harmonics_expression
from multipiegen.util.util_info import load_info
from multipiegen.util.util_representation_matrix import _group, representation_matrix

info = load_info()
data = load_multipolar_harmonics()
SAMPLE = ["Oh", "D6h", "D3d", "C3h-c", "T-c", "D2d-1"]
_x, _y, _z = sp.symbols("x y z", real=True)
_x2, _y2, _z2 = sp.symbols("X2 Y2 Z2", real=True)


# ==================================================
def test_all_groups():
    ids = info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]
    assert list(data.keys()) == ids


# ==================================================
@pytest.mark.parametrize("pg", SAMPLE)
def test_binary_is_up_to_date(pg):
    # compare after a pickle round trip, as sympy re-evaluates some products (e.g. (a+b)/2 -> a/2+b/2) on loading.
    new = pickle.loads(pickle.dumps(create_multipolar_harmonics(pg)))
    assert repr(new) == repr(data[info["id"][pg]])


# ==================================================
@pytest.mark.parametrize("pg", SAMPLE)
def test_keys_and_values(pg):
    crystal = _crystal[data_harmonics_mapping[pg]["top"]]
    irreps = list(character[pg][1])
    keys = list(data[info["id"][pg]].keys())
    order = [(["q", "g"].index(k.x), k.s, irreps.index(k.Gamma), k.l, ["Q", "G"].index(k.X), k.k, k.n) for k in keys]
    assert order == sorted(order)
    for k in keys:
        coeff, expr = data[info["id"][pg]][k]
        ex = harmonics_expression(pg, k)
        assert 1 <= k.s <= max_s_rank and k.l <= max_rank and abs(k.k) <= k.s and k.p == -1
        assert k.x == internal_type(k.X, k.s, k.k)
        assert coeff.shape == (len(ex), 2 * k.s + 1) and expr.shape == ex.shape
        sym = internal_symbols(crystal)[(k.s, k.x)]
        assert all(sp.expand(e - sum(c * s for c, s in zip(row, sym))) == 0 for e, row in zip(expr, coeff))


# ==================================================
@pytest.mark.parametrize("pg", SAMPLE)
def test_representation_matrix(pg):
    # F_c(r, r') = sum_a coeff[c,a](r) phi_a(r') (phi: internal orbitals) transforms with the representation
    # matrices of its irrep: s F_c(R^T r, R^T r') = sum_c' D_c'c F_c'(r, r'), s = det R for axial internal basis.
    crystal = _crystal[data_harmonics_mapping[pg]["top"]]
    bf = info["harmonics"]["basis_function"]
    phi = {
        s: [bf[n][0].subs({_x: _x2, _y: _y2, _z: _z2}, simultaneous=True) for n in _orbital_basis[crystal][s]] for s in (1, 2, 3)
    }
    _, R, det = _group(pg)
    D = {ir: m.astype(complex) for ir, m in representation_matrix(pg)["matrix"].items()}
    rng = np.random.default_rng(0)
    P1, P2 = rng.normal(size=(4, 3)), rng.normal(size=(4, 3))
    for k, (coeff, expr) in data[info["id"][pg]].items():
        f = sp.lambdify((_x, _y, _z, _x2, _y2, _z2), [sum(m * p for m, p in zip(row, phi[k.s])) for row in coeff], "numpy")

        def F(A, B):
            return np.array([np.broadcast_to(np.asarray(v, dtype=complex), (len(A),)) for v in f(*A.T, *B.T)])

        V = F(P1, P2)
        for g in range(len(R)):
            s = det[g] if k.x == "g" else 1.0
            assert np.allclose(s * F(P1 @ R[g], P2 @ R[g]), D[k.Gamma][g].T @ V), (k, g)
