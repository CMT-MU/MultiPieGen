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
