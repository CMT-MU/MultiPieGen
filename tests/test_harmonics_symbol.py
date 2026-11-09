# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pytest
from multipiegen.data.data_group_pg import pg_id_set, pg_tag
from multipiegen.util.util_harmonics_symbol import harmonics_symbol_table, polar_to_axial, multipole_symbol, irrep_latex

all_pg = [pg_tag[i] for i in pg_id_set["all"] + pg_id_set["complex"]]


def _swap_gu(irrep):
    return irrep[:-1] + {"g": "u", "u": "g"}[irrep[-1]]


def test_number_of_groups():
    assert len(all_pg) == 47


@pytest.mark.parametrize("group", all_pg)
def test_table(group):
    tbl = harmonics_symbol_table(group)
    assert len(tbl) == 144  # sum of (2l+1) for l = 0-11.
    assert len({r[1] for r in tbl}) == 144  # polar symbols are unique.
    assert len({r[2].lstrip("-") for r in tbl}) == 144  # axial symbols are unique.
    assert [r[0][0] for r in tbl] == sorted(r[0][0] for r in tbl)  # sorted by rank.


@pytest.mark.parametrize("group", ["Oh", "D6h"])
def test_top_group_axial(group):
    # axial of top group: g <-> u, same component, no sign.
    for key, p_sym, a_sym, expr in harmonics_symbol_table(group):
        a_key, sign = polar_to_axial(group, key)
        assert a_key == (key[0], _swap_gu(key[1]), key[2], key[3])
        assert sign == 1


def test_subgroup_sign():
    # Td: polar E (comp 2) is axial -G(E) comp 1.
    a_key, sign = polar_to_axial("Td", (2, "E", -1, 1))
    assert a_key == (2, "E", -1, 0) and sign == -1


def test_symbol_format():
    assert multipole_symbol("Q", (0, "A1g", -1, -1)) == "\\mathbb{Q}_{0}(A_{1g})"
    assert multipole_symbol("Q", (1, "T1u", -1, 0)) == "\\mathbb{Q}_{1,1}(T_{1u})"
    assert multipole_symbol("G", (11, "Eu", 2, 1), -1) == "-\\mathbb{G}_{11,2}(E_{u},2)"
    assert irrep_latex("A2''") == "A_{2}^{\\prime\\prime}"
    assert irrep_latex("E1ua") == "E_{1u}^{(a)}"
    assert irrep_latex("Ea'") == "E^{\\prime(a)}"
    assert irrep_latex("Ag") == "A_{g}"
