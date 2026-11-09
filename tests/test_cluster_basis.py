# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest
from multipiegen.data.data_character_table import character
from multipiegen.util.util_cluster_basis import load_cluster_basis, cluster_basis_matrix, create_cluster_basis

data = load_cluster_basis()


@pytest.mark.parametrize("group, n", [("Oh", 48), ("D6h", 24)])
def test_orthonormal(group, n):
    m, keys = cluster_basis_matrix(data[group])
    assert m.shape == (n, n)
    assert np.allclose(m @ m.T, np.eye(n))


@pytest.mark.parametrize("group", ["Oh", "D6h"])
def test_regular_representation(group):
    # every irrep appears, as many times as its dimension (character of identity).
    count = {}
    for key, v in data[group]["basis"].items():
        assert (key.X, key.p, key.s, key.k, key.x) == ("Q", -1, 0, 0, "q")  # PGMultipoleType as MultiPie site SAMB.
        count[key.Gamma] = count.get(key.Gamma, 0) + len(v)
    table = character[group][1]
    assert set(count) == set(table)
    for irrep, chi in table.items():
        d = int(chi.strip("[]").split(",")[0])
        assert count[irrep] == d * d, irrep


@pytest.mark.parametrize("group", ["Oh", "D6h"])
def test_binary_is_up_to_date(group):
    new = create_cluster_basis(group)
    old = data[group]
    assert list(new["basis"].keys()) == list(old["basis"].keys())
    assert new["dropped"] == old["dropped"]
    for k in new["basis"]:
        assert repr(new["basis"][k]) == repr(old["basis"][k]), k
