# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for the mapping part of create_data_harmonics.py (no binary data / sympy calculation needed).

Usage:
    pytest test_create_data_harmonics.py
(place this file in the same directory as create_data_harmonics.py.)
"""

import json
import sys

import pytest

import create_data_harmonics as m

ALL_TAGS = list(m.compatibility_relation)


# ==================================================
# create_mapping.
# ==================================================
@pytest.mark.parametrize("tag", ALL_TAGS)
def test_every_top_harmonics_used_once(tag):
    """Each top harmonics is used exactly once for each of Q and G."""
    top = m.compatibility_relation[tag][0]
    mapping = m.create_mapping(tag)
    for X in m.HARMONICS_TYPES:
        used = [tk for data in mapping[X].values() for _, tk in data]
        assert len(used) == len(m.top_harmonics_tag[top])
        assert len(set(used)) == len(used)


@pytest.mark.parametrize("tag", ALL_TAGS)
def test_coefficients_are_pm1(tag):
    mapping = m.create_mapping(tag)
    for X in m.HARMONICS_TYPES:
        for data in mapping[X].values():
            assert all(c in (1, -1) for c, _ in data)


def test_oh_identity_and_axial_gu_swap():
    """Oh: polar keeps the irrep, axial interchanges g and u."""
    mapping = m.create_mapping("Oh")
    assert mapping["Q"][(2, "T2g", 1)] == [(1, ("Oh", 2, "T2g", -1, 1))]
    assert mapping["G"][(2, "T2u", 1)] == [(1, ("Oh", 2, "T2g", -1, 1))]
    assert mapping["Q"][(3, "A2u", -1)] == [(1, ("Oh", 3, "A2u", -1, -1))]
    assert mapping["G"][(3, "A2g", -1)] == [(1, ("Oh", 3, "A2u", -1, -1))]


def test_td_values():
    mapping = m.create_mapping("Td")
    assert mapping["Q"][(2, "E", 0)] == [(1, ("Oh", 2, "Eg", -1, 0))]
    assert mapping["Q"][(2, "T2", 1)] == [(1, ("Oh", 2, "T2g", -1, 1))]
    # axial: components of E are interchanged with a sign change.
    assert mapping["G"][(2, "E", 1)] == [(1, ("Oh", 2, "Eg", -1, 0))]
    assert mapping["G"][(2, "E", 0)] == [(-1, ("Oh", 2, "Eg", -1, 1))]
    assert mapping["G"][(2, "T1", 0)] == [(1, ("Oh", 2, "T2g", -1, 0))]


def test_complex_group_uses_ab_expression():
    """"-c" groups refer to the complex (a/b) top harmonics with a single component."""
    mapping = m.create_mapping("C6h-c")
    assert mapping["Q"][(1, "E1ua", -1)] == [(1, ("D6h", 1, "E1ua", -1, -1))]
    assert mapping["Q"][(1, "E1ub", -1)] == [(1, ("D6h", 1, "E1ub", -1, -1))]
    # non-degenerate irrep is not affected.
    assert mapping["Q"][(1, "Au", -1)] == [(1, ("D6h", 1, "A2u", -1, -1))]


def test_real_group_does_not_use_ab_expression():
    mapping = m.create_mapping("C6h")
    tops = [tk for data in mapping["Q"].values() for _, tk in data]
    assert all(not tk[2].endswith(("a", "b")) for tk in tops)


def test_top_group_of_d6h_subgroup():
    assert m.compatibility_relation["D3h"][0] == "D6h"
    mapping = m.create_mapping("D3h")
    assert all(tk[0] == "D6h" for data in mapping["Q"].values() for _, tk in data)


def test_unknown_group_raises():
    with pytest.raises(KeyError):
        m.create_mapping("XXX")


# ==================================================
# _top_irrep_comp.
# ==================================================
@pytest.mark.parametrize(
    "irrep, comp, sub_irrep, expected",
    [
        ("Eg", 0, "Ega", ("Ega", -1)),
        ("Eg", 1, "Egb", ("Egb", -1)),
        ("Eg", 0, "Eg", ("Eg", 0)),  # real: unchanged.
        ("T1u", 2, "Au", ("T1u", 2)),  # comp 2 is not converted.
        ("A1g", -1, "Ag", ("A1g", -1)),
    ],
)
def test_top_irrep_comp(irrep, comp, sub_irrep, expected):
    assert m._top_irrep_comp(irrep, comp, sub_irrep) == expected


# ==================================================
# renumber_multiplicity.
# ==================================================
def test_renumber_single_is_minus_one():
    renum = m.renumber_multiplicity(m.create_mapping("D4h"))
    assert renum["Q"][(2, "A1g", -1, -1)] == (1, ("Oh", 2, "Eg", -1, 0))
    assert renum["Q"][(2, "B1g", -1, -1)] == (1, ("Oh", 2, "Eg", -1, 1))


def test_renumber_multiplicity_follows_tag_order():
    """D4h l=4: A1g comes from Oh A1g and Oh Eg(0), in this order."""
    renum = m.renumber_multiplicity(m.create_mapping("D4h"))
    assert renum["Q"][(4, "A1g", 1, -1)] == (1, ("Oh", 4, "A1g", -1, -1))
    assert renum["Q"][(4, "A1g", 2, -1)] == (1, ("Oh", 4, "Eg", -1, 0))
    assert (4, "A1g", -1, -1) not in renum["Q"]


@pytest.mark.parametrize("tag", ALL_TAGS)
def test_renumber_keeps_all_entries(tag):
    mapping = m.create_mapping(tag)
    renum = m.renumber_multiplicity(mapping)
    for X in m.HARMONICS_TYPES:
        assert len(renum[X]) == sum(len(v) for v in mapping[X].values())
        # multiplicity is -1 or consecutive 1..N for each (l, irrep, comp).
        groups = {}
        for (l, irrep, n, comp) in renum[X]:
            groups.setdefault((l, irrep, comp), []).append(n)
        for ns in groups.values():
            assert ns == [-1] or sorted(ns) == list(range(1, len(ns) + 1))


# ==================================================
# create_mappings.
# ==================================================
def test_create_mappings_default_all():
    assert list(m.create_mappings()) == ALL_TAGS


def test_create_mappings_selected():
    result = m.create_mappings(["Td", "D4h"])
    assert list(result) == ["Td", "D4h"]
    assert result["Td"] == m.create_mapping("Td")


# ==================================================
# output (text / json / CLI).
# ==================================================
def test_format_mapping_text():
    text = m.format_mapping("Td", m.create_mapping("Td"))
    assert text.startswith("[Td] top=Oh")
    assert "(2, 'E', 1) <- +1 * ('Oh', 2, 'Eg', -1, 0)" in text
    assert "(2, 'E', 0) <- -1 * ('Oh', 2, 'Eg', -1, 1)" in text


def test_format_mapping_text_renumbered():
    renum = m.renumber_multiplicity(m.create_mapping("D4h"))
    text = m.format_mapping("D4h", renum, renumbered=True)
    assert "(l, irrep, n, comp)" in text
    assert "(4, 'A1g', 2, -1) <- +1 * ('Oh', 4, 'Eg', -1, 0)" in text


def test_mapping_to_jsonable_roundtrip():
    data = m.mapping_to_jsonable(m.create_mapping("Td"))
    data = json.loads(json.dumps(data))  # must be serializable.
    assert data["G"]["2,E,0"] == [{"coeff": -1, "top": ["Oh", 2, "Eg", -1, 1]}]


def test_mapping_to_jsonable_renumbered():
    renum = m.renumber_multiplicity(m.create_mapping("D4h"))
    data = m.mapping_to_jsonable(renum)
    assert data["Q"]["4,A1g,2,-1"] == [{"coeff": 1, "top": ["Oh", 4, "Eg", -1, 0]}]


def test_output_mapping_stdout(capsys):
    m.output_mapping(["Td"])
    out = capsys.readouterr().out
    assert "[Td]" in out and "[D4h]" not in out


def test_output_mapping_json_file(tmp_path):
    path = tmp_path / "map.json"
    m.output_mapping(["Td", "D4h"], renumber=True, fmt="json", output=str(path))
    data = json.loads(path.read_text(encoding="utf-8"))
    assert set(data) == {"Td", "D4h"}
    assert "4,A1g,1,-1" in data["D4h"]["Q"]


def test_output_mapping_unknown_group():
    with pytest.raises(SystemExit):
        m.output_mapping(["XXX"])


def test_main_mapping_option(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["create_data_harmonics.py", "--mapping", "-g", "Td", "--renumber"])
    m.main()
    out = capsys.readouterr().out
    assert "[Td]" in out and "(l, irrep, n, comp)" in out


def test_main_default_creates_binary(monkeypatch):
    called = []
    monkeypatch.setattr(m, "create_harmonics", lambda: called.append(True))
    monkeypatch.setattr(sys, "argv", ["create_data_harmonics.py"])
    m.main()
    assert called == [True]
