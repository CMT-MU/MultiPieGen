# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import sympy as sp
from multipie.util.util import str_to_sympy
from multipiegen.create_binary.create_info import create_info
from multipiegen.data.data_basis import atomic_orbital_basis, cubic_orbital_basis, hexagonal_orbital_basis
from multipiegen.data.data_top_harmonics_info import top_harmonics_symmetry
from multipiegen.util.util_info import load_info

info = load_info()


def test_binary_is_up_to_date():
    new = create_info(verbose=False)
    assert list(new.keys()) == list(info.keys())
    for k in new.keys():
        assert repr(new[k]) == repr(info[k]), k
    assert str(new) == str(info)


def test_group():
    assert len(info["id_set"]["PG"]["all"]) + len(info["id_set"]["PG"]["complex"]) == 47
    assert len(info["id_set"]["SG"]["all"]) == 230
    for g_type in ["PG", "SG", "MPG", "MSG"]:
        for id_s, name in info["group_name"][g_type].items():
            # name list uses "_" instead of "^" in SG tags (e.g., C1_1 for C1^1).
            assert info["tag"][id_s].replace("^", "_") == name["tag"]
            assert name[g_type] == id_s


def test_tesseral_tag():
    # basis function of each regular orbital equals the tesseral harmonics given by its tesseral tag.
    h = info["harmonics"]
    for name, (l, m, cs) in h["tesseral"].items():
        expr = h["basis_function"][name][0]
        t = str_to_sympy(top_harmonics_symmetry["tesseral"][(l, m if cs == "c" else -m)])
        assert sp.expand(sp.sympify(str(expr)) - sp.sympify(str(t))) == 0, name


def test_orbitals_have_basis_function_and_latex():
    h = info["harmonics"]
    for basis in [atomic_orbital_basis, cubic_orbital_basis, hexagonal_orbital_basis]:
        for names in basis:
            for name in names:
                assert name in h["basis_function"]
                assert name in h["latex"]["orbital"]
                assert len(h["basis_function"][name][1]) == 2 * h["harmonics_key"][name][1][0] + 1
