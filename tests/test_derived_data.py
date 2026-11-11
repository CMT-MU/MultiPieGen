# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Data derived from the tables in multipiegen/data (computed in multipiegen, not stored)."""

import numpy as np
import pytest
from multipiegen.data.data_group_mpg import mpg_info
from multipiegen.data.data_group_msg import msg_info
from multipiegen.data.data_group_pg import pg_info
from multipiegen.data.data_group_sg import sg_info
from multipiegen.util.util_group import load_group
from multipiegen.util.util_group_data import group_name_list, msg_id_set, pg_id, pg_id_set, pg_tag, sg_id_set
from multipiegen.util.util_info import load_info
from multipiegen.util.util_representation_matrix import representation_matrix
from multipiegen.util.util_top_harmonics import top_harmonics_symmetry, top_harmonics_u_matrix

info = load_info()
group = load_group()


# ==================================================
def test_group_ids():
    assert pg_tag == {k: v.tag for k, v in pg_info.items()} and pg_id == {v: k for k, v in pg_tag.items()}
    assert sorted(pg_id_set["all"] + pg_id_set["complex"]) == sorted(pg_info)
    assert sorted(sum(pg_id_set["irrep"].values(), [])) == sorted(pg_info)
    assert sorted(sum(sg_id_set["PG"].values(), [])) == sorted(sg_info)
    assert all(sg_info[i].PG == pg for pg, lst in sg_id_set["PG"].items() for i in lst)
    assert sorted(sum(msg_id_set["type"].values(), [])) == sorted(msg_info)
    assert [list(group_name_list[t]) for t in ["PG", "SG", "MPG", "MSG"]] == [
        list(pg_info),
        list(sg_info),
        list(mpg_info),
        list(msg_info),
    ]


# ==================================================
def test_top_harmonics():
    assert list(top_harmonics_symmetry) == list(top_harmonics_u_matrix) == ["tesseral", "Oh", "D6h"]
    for top in ["Oh", "D6h"]:
        assert list(top_harmonics_symmetry[top]) == list(top_harmonics_u_matrix[top])


# ==================================================
@pytest.mark.parametrize("pg", ["C4", "S4", "C4h", "C3", "C6h", "T", "Th", "D4h", "Oh", "D6h", "C3h-c", "C4-c"])
def test_character_plus(pg):
    # products of irreps in group.pkl agree with those from the representation matrices (independent of the characters).
    ch = group[info["id"][pg]]["character"]
    rm = representation_matrix(pg)
    D = {ir: np.array(m, dtype=complex) for ir, m in rm["matrix"].items()}
    chi = {ir: np.array([np.trace(d) for d in m]) for ir, m in D.items()}

    def dec(c):
        return sorted(
            (int(round((np.vdot(x, c) / np.vdot(x, x)).real)), ir) for ir, x in chi.items() if abs(np.vdot(x, c)) > 1e-8
        )

    for (a, b), v in ch["symmetric_product"].items():
        c = np.array([(np.trace(m) ** 2 + np.trace(m @ m)) / 2 for m in D[a]]) if a == b else chi[a] * chi[b]
        assert dec(c) == sorted(v), (a, b)
    for a, v in ch["anti_symmetric_product"].items():
        assert dec(np.array([(np.trace(m) ** 2 - np.trace(m @ m)) / 2 for m in D[a]])) == sorted(v), a
    det = np.array(
        [np.linalg.det(np.array(m, dtype=float)[:3, :3]) for m in group[info["id"][pg]]["symmetry_operation"]["cartesian"]]
    )
    for a, b in ch["polar_axial_conversion"].items():
        assert np.allclose(det * chi[a], chi[b]), (a, b)
