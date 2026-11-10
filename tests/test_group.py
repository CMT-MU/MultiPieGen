# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest
import sympy as sp
from multipie import Group
from multipiegen.create_binary.create_group import create_group_data, group_ids
from multipiegen.data.data_op_mapping import op_mapping_d6h_oh_cluster
from multipiegen.data.data_wyckoff_bond_sym import wyckoff_bond_data
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info
from multipiegen.util.util_symmetry_operation import simplify_value

info = load_info()
group = load_group()

# representative groups: cubic/hexagonal PG, complex PG, SG with plus set, hexagonal SG, MPG, MSG.
SAMPLE = ["PG:32", "PG:27", "PG:44", "SG:225", "SG:194", "MPG:24.5.149", "MSG:177.151"]


# ==================================================
def _norm(x):
    """Comparable form of nested data (arrays by shape and str)."""
    if isinstance(x, dict) or (hasattr(x, "items") and not hasattr(x, "_fields")):
        return {k: _norm(v) for k, v in x.items()}
    if isinstance(x, np.ndarray):
        return (x.shape, str(x.tolist()))
    if isinstance(x, list):
        return [_norm(i) for i in x]
    return repr(x)


# ==================================================
def _same(a, b):
    a, b = np.asarray(a, dtype=object), np.asarray(b, dtype=object)
    return a.shape == b.shape and all(sp.simplify(sp.S(x) - sp.S(y)) == 0 for x, y in zip(a.ravel(), b.ravel()))


# ==================================================
def test_group_ids():
    ids = [i for i, _ in group_ids(info)]
    assert list(group.keys()) == ids
    n = {k: len([i for i in ids if i.split(":")[0] == k]) for k in ["PG", "SG", "MPG", "MSG"]}
    assert n == {"PG": 47, "SG": 230, "MPG": 157, "MSG": 1651}


# ==================================================
@pytest.mark.parametrize("gid", SAMPLE)
def test_binary_is_up_to_date(gid):
    assert info["tag"]["PG:44"] == "C3h-c"
    assert _norm(create_group_data(gid)) == _norm(group[gid])


# ==================================================
@pytest.mark.parametrize("gid", ["PG:32", "PG:27", "PG:44", "PG:16", "SG:225", "SG:194", "SG:167"])
def test_agree_with_multipie(gid):
    # items that do not depend on the new root cluster agree with MultiPie.
    g, mp = group[gid], Group(info["tag"][gid])
    so, mso = g["symmetry_operation"], mp.symmetry_operation
    assert so["tag"] == mso["tag"]
    for k in ["cartesian", "fractional", "det"]:
        assert _same(so[k], mso[k]), k
    for wp, v in mp.wyckoff["site"].items():
        assert v["mapping"] == g["wyckoff"]["site"][wp]["mapping"], wp
        assert _same(v["conventional"], g["wyckoff"]["site"][wp]["conventional"]), wp
    for bw, v in mp.wyckoff.get("bond", {}).items():
        assert v["mapping"] == g["wyckoff"]["bond"][bw]["mapping"], bw
    if "character" in g:
        for irrep, ch in mp.character["table_full"].items():
            assert _same(ch, g["character"]["table_full"][irrep]), irrep


# ==================================================
def test_product_and_mapping():
    top = {"Oh": group["PG:32"]["symmetry_operation"], "D6h": group["PG:27"]["symmetry_operation"]}
    for gid in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"] + info["id_set"]["SG"]["all"]:
        so = group[gid]["symmetry_operation"]
        tags = [t.split(":")[0] for t in so["tag"]]
        if "product" in so:  # closed product table
            assert set(so["product"].values()) == set(so["tag"])
            assert len(so["product"]) == len(so["tag"]) ** 2
        mp = so["mapping"]
        t_tag = top[mp["top"]]["tag"]
        assert sorted(sum(mp["mapping"], [])) == list(range(1, len(t_tag) + 1))  # cosets cover the top group once
        assert [t_tag[i - 1] for i in mp["mapping"][0]] == [t.split(":")[0] for t in tags]  # first coset is the group
        pg_id = gid if gid.startswith("PG") else group[gid]["info"].PG
        key = info["tag"][pg_id].replace("-c", "")
        assert mp == {"top": op_mapping_d6h_oh_cluster[key][0], "mapping": op_mapping_d6h_oh_cluster[key][1]}


# ==================================================
def test_complex_pg_has_real_data():
    for gid in info["id_set"]["PG"]["complex"]:
        tag = info["tag"][gid]
        rid = info["id"][tag.replace("-c", "")]
        for k in ["symmetry_operation", "wyckoff"]:
            assert _norm(group[gid][k]) == _norm(group[rid][k]), (tag, k)
        assert list(group[gid]["character"]["table"]) != list(group[rid]["character"]["table"])


# ==================================================
def test_bond_symmetry():
    for gid in info["id_set"]["PG"]["all"] + info["id_set"]["SG"]["all"]:
        tag = info["tag"][gid]
        bond = group[gid]["wyckoff"]["bond"]
        assert list(bond) == list(wyckoff_bond_data[tag])
        for bw, v in bond.items():
            assert v["symmetry"] == wyckoff_bond_data[tag][bw][3]
            assert bw in group[gid]["wyckoff"]["site"][bw.split("@")[1]]["bond"]


# ==================================================
@pytest.mark.parametrize("gid", SAMPLE)
def test_values_are_simplest(gid):
    def walk(x):
        if isinstance(x, dict) or (hasattr(x, "items") and not hasattr(x, "_fields")):
            for k, v in x.items():
                if k != "table":  # "table" keeps wp/wm.
                    walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, np.ndarray) and x.dtype == object:
            for v in x.ravel():
                if isinstance(v, sp.Basic):
                    assert simplify_value(v) == v, (gid, v)

    walk(group[gid])


# ==================================================
def test_op_mapping_parent_sg():
    # SG -> top parent SG (221-230 or 191-194): cosets cover the parent SG once, and the top SG maps to itself.
    from multipiegen.data.data_op_mapping import op_mapping_parent_sg

    assert sorted(op_mapping_parent_sg) == list(range(1, 231))
    for no, (top, mapping) in op_mapping_parent_sg.items():
        assert top in list(range(221, 231)) + list(range(191, 195))
        n_top = len(group[f"SG:{top}"]["symmetry_operation"]["tag"])
        assert sorted(sum(mapping, [])) == list(range(1, n_top + 1))
        assert len(mapping[0]) == len(group[f"SG:{no}"]["symmetry_operation"]["tag"])
        if no == top:
            assert mapping == [list(range(1, n_top + 1))]
