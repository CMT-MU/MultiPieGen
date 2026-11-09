# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest
from multipiegen.data.data_group_pg import pg_id_set, pg_tag
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info
from multipiegen.util.util_representation_matrix import check_representation_matrix, representation_matrix

PG = [pg_tag[i] for i in pg_id_set["all"] + pg_id_set["complex"]]
TOL = 1e-9


# ==================================================
@pytest.mark.parametrize("pg", PG)
def test_same_matrix_in_each_irrep(pg):
    # all polar and axial harmonics (up to rank 11) of the same irrep have the same representation matrices.
    for irrep, d in check_representation_matrix(pg).items():
        assert d["closure"] < TOL, (pg, irrep, d["closure"])
        assert d["deviation"] < TOL, (pg, irrep, d["deviation"])


# ==================================================
@pytest.mark.parametrize("pg", ["Oh", "D6h", "D3h", "C3h-c", "T-c", "D3d-1"])
def test_exact_matrix(pg):
    # trace = character, and D(a) D(b) = D(ab).
    gid = load_info()["id"][pg]
    g = load_group()[gid]
    rm = representation_matrix(pg)
    idx = {t: i for i, t in enumerate(rm["tag"])}
    prod = g["symmetry_operation"]["product"]
    assert list(rm["matrix"]) == list(g["character"]["table"])
    for irrep, D in rm["matrix"].items():
        Df = D.astype(complex)
        assert np.allclose([np.trace(d) for d in Df], g["character"]["table_full"][irrep].astype(complex), atol=TOL), irrep
        for (a, b), ab in prod.items():
            assert np.allclose(Df[idx[a]] @ Df[idx[b]], Df[idx[ab]], atol=TOL), (irrep, a, b)
