# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
from multipie.util.util import str_to_sympy
from multipie import Group

TOL = 1e-6


# ==================================================
def same(a, b):
    """3x3 operations a and b are equal."""
    return bool(np.abs(a - b).max() < TOL)


# ==================================================
def check(c_group, p_group):
    cg = Group(c_group)
    pg = Group(p_group)
    cosets, R = data[c_group][p_group]  # cosets[0] = parent operation numbers (1-based) of child operations
    R = np.array(str_to_sympy(R).astype(float))  # 3x3, x_parent = M x_child
    Ri = np.linalg.inv(R)

    p_so = pg.symmetry_operation["fractional"].astype(float)  # (n_p, 3, 3)
    c_so = cg.symmetry_operation["fractional"].astype(float)  # (n_c, 3, 3)

    err = []
    op = list(cosets[0])
    n_c, n_p = len(c_so), len(p_so)

    def find(a):  # 1-based number of the parent operation equal to a
        for k, b in enumerate(p_so):
            if same(a, b):
                return k + 1
        return None

    # --- 1. number / injectivity / M -------------------------------------------------------
    if abs(np.linalg.det(R)) < TOL:
        err.append("M is singular.")
        return err
    if len(op) != n_c:
        err.append(f"len(op)={len(op)} but the child point group has {n_c} operations.")
    if len(set(op)) != len(op):
        err.append("op is not injective.")

    # --- 2. each child operation, transformed to the parent's coordinates, is the mapped parent operation
    for no, m in enumerate(op):
        if not same(R @ c_so[no] @ Ri, p_so[m - 1]):
            err.append(f"child op {no + 1} -> parent op {m}: linear part differs.")

    # --- 3. the image is closed under products (a subgroup of the parent point group)
    image = set(op)
    for a in op:
        for b in op:
            k = find(p_so[a - 1] @ p_so[b - 1])
            if k not in image:
                err.append(f"image not closed: parent ops {a}, {b} give {k}.")

    # --- 4. coset decomposition: every parent operation exactly once, row k = g_k o (row 0)
    flat = [x for row in cosets for x in row]
    if sorted(flat) != list(range(1, n_p + 1)):
        err.append("cosets do not contain every parent operation exactly once.")
    for row in cosets[1:]:
        g = row[0]
        for a, x in zip(op, row):
            k = find(p_so[g - 1] @ p_so[a - 1])
            if k != x:
                err.append(f"coset row {row[:4]}...: g({g}) o a({a}) = {k}, not {x}.")
                break
    return err


# ==================================================
from multipiegen.data.data_parent_pg import data_parent_pg as data


def test_parent_pg_mapping():
    err = 0
    for c_group, v in data.items():
        for p_group in v.keys():
            lst = check(c_group, p_group)
            if lst:
                print(f"error in {p_group} => {c_group}.")
                for i in lst:
                    print("   ", i)
                err += 1
    assert err == 0
