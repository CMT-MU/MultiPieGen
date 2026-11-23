# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import re

import pytest
import sympy as sp
from multipiegen.create_binary.create_response_tensor import create_response_tensor_group
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info
from multipiegen.util.util_response_tensor import load_response_tensor, multipole_latex, tensor_shape, tensor_type

info = load_info()
group = load_group()
rt = load_response_tensor()


# ==================================================
def test_all_groups():
    assert list(rt.keys()) == info["id_set"]["MPG"]["all"]
    keys = [(X, *k) for X in tensor_type for k in tensor_shape]
    for gid in rt.keys():
        assert list(rt[gid].keys()) == keys, gid


# ==================================================
@pytest.mark.parametrize("gid", ["MPG:1.1.1", "MPG:24.5.149", "MPG:32.5.122", "MPG:27.9.157"])
def test_binary_is_up_to_date(gid):
    assert repr(create_response_tensor_group(gid)) == repr(rt[gid])


# ==================================================
def test_only_active_multipoles():
    # the expressions use only the active multipoles of the group, and every active multipole appears.
    latex = info["harmonics"]["latex"]["orbital"]

    def base(s):  # symbol without the label (n) of the cartesian decomposition.
        return re.sub(r"\(\d+\)", "", str(s))

    for gid in rt.keys():
        allowed = {base(multipole_latex(a[0], a[1:], "", latex)) for a in group[gid]["active_multipole"]}
        used = {base(s) for m, ex in rt[gid].values() for e in ex.values() for s in e.free_symbols}
        assert used == allowed, (gid, used ^ allowed)


# ==================================================
@pytest.mark.parametrize("crystal", ["cubic", "hexagonal"])
def test_cartesian_multipole_positive(crystal):
    # the traceless part of each cartesian monomial is sum_k c_k N_k f_k with the coefficients c_k of cartesian_multipole
    # and the basis functions f_k: every N_k is positive and the same for all keys, i.e., every multipole in the
    # response tensors is a positive multiple of the harmonic of the same name.
    x, y, z = sp.symbols("x y z")
    v = {"x": x, "y": y, "z": z}
    r2 = x**2 + y**2 + z**2
    bf = info["harmonics"]["basis_function"]

    def f(name):
        e = bf[name][0]
        return sp.expand(e.subs({s: v[str(s)] for s in e.free_symbols}))

    def traceless(p, l):
        h = 0
        for k in range(l // 2 + 1):
            h += (
                sp.Integer(-1) ** k
                * sp.factorial2(2 * l - 2 * k - 1)
                / (sp.factorial2(2 * l - 1) * sp.factorial2(2 * k))
                * r2**k
                * p
            )
            p = sp.diff(p, x, 2) + sp.diff(p, y, 2) + sp.diff(p, z, 2)
        return sp.expand(h)

    N = {}
    for key, lst in info["response_tensor"]["cartesian_multipole"][crystal].items():
        h = traceless(sp.Mul(*[v[c] for c in key]), len(key)) if key != "s" else sp.S(1)
        n = sp.symbols(f"n0:{len(lst)}")
        sol = sp.solve(sp.Poly(h - sum(c * ni * f(name) for (c, name), ni in zip(lst, n)), x, y, z).coeffs(), n, dict=True)
        assert len(sol) == 1, key
        for (c, name), ni in zip(lst, n):
            N.setdefault(name, set()).add(sp.nsimplify(sol[0][ni]))
    for name, s in N.items():
        assert len(s) == 1 and s.pop() > 0, name
