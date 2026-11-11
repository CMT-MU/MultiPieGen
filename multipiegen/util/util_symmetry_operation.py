# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Symmetry operation matrices (sympy).

Ported from MultiPie Ver.2 (others/tools/utils/util_symmetry_operation.py and util_symmetry_operation_spherical.py).
The matrices of the point-group part are cached, because the same operations appear in many groups.
group_symmetry_operation gives the symmetry operations of a PG or SG from the group information tables
(data_group_pg/sg) alone, so that the table generators in others/creation do not depend on other data.
"""

from functools import lru_cache

import numpy as np
import sympy as sp
from sympy.physics.quantum.spin import Rotation
from multipie.util.util import str_to_sympy
from multipie.util.util_crystal import convert_to_fractional_hexagonal_matrix
from multipie.util.util_tag import TagSymmetryOperation

from multipiegen.data.data_group_pg import pg_info
from multipiegen.data.data_group_sg import sg_info
from multipiegen.util.util_group_data import pg_id, sg_id
from multipiegen.util.util_multipolar_harmonics import _create_unitary_matrix_tesseral


# ==================================================
def simplify_value(x):
    """
    Bring a sympy value into a simple canonical form.

    Args:
        x (sympy or int): value, a number or an expression with symbols.

    Returns:
        - (sympy) -- expanded value with rationalized denominators (e.g. sqrt(3)/3 instead of 1/sqrt(3)).
    """
    x = sp.S(x)
    if x.is_Rational:
        return x
    return sp.expand(sp.radsimp(sp.expand(x)))


# ==================================================
def simplify_array(a):
    """
    Apply simplify_value to all elements of an array.

    Args:
        a (array-like): array of sympy values.

    Returns:
        - (ndarray) -- array of simplified values (dtype=object).
    """
    a = np.asarray(a, dtype=object)
    if a.size == 0:
        return a
    return np.vectorize(simplify_value, otypes=[object])(a)


# ==================================================
def _rotation_matrix(theta, axis):
    """
    Rotation matrix (cartesian coordinate).

    Args:
        theta (sympy): rotation angle (radian).
        axis (array-like): rotation axis (cartesian coordinate, no need for normalization).

    Returns:
        - (ndarray) -- 3x3 rotation matrix, [[sympy]].
    """
    n = np.asarray([sp.S(axis[0]), sp.S(axis[1]), sp.S(axis[2])])
    norm = sp.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
    n /= norm

    K = np.array([[sp.S(0), -n[2], n[1]], [n[2], sp.S(0), -n[0]], [-n[1], n[0], sp.S(0)]])
    I = np.full((3, 3), sp.S(0))
    I[0, 0] = I[1, 1] = I[2, 2] = sp.S(1)
    R = I + sp.sin(theta) * K + (1 - sp.cos(theta)) * K @ K

    return R


# ==================================================
def _axis_to_euler_zyz(theta, axis):
    """
    Convert (theta, axis) to zyz-Euler angle, (alpha, beta, gamma).

    Args:
        theta (sympy): rotation angle (radian).
        axis (array-like): rotation axis (cartesian coordinate, no need for normalization).

    Returns:
        - (tuple) -- angles, ( alpha(sympy), beta(sympy), gamma(sympy) ).
    """
    R = _rotation_matrix(theta, axis)

    if R[2, 2] < 1 and R[2, 2] > -1:  # |Rzz| ne 1.
        alpha = sp.atan2(R[1, 2], R[0, 2])
        beta = sp.acos(R[2, 2])
        gamma = sp.atan2(R[2, 1], -R[2, 0])
    else:  # gimbal lock.
        alpha = sp.atan2(-R[0, 1], R[1, 1])
        beta = sp.S(0) if R[2, 2] == sp.S(1) else sp.pi
        gamma = sp.S(0)

    return alpha, beta, gamma


# ==================================================
def _wigner_D_matrix(j, theta, axis):
    """
    Wigner D matrix.

    Args:
        j (int): angular momentum.
        theta (sympy): rotation angle (radian).
        axis (array-like): rotation axis (cartesian coordinate, no need for normalization).

    Returns:
        - (ndarray) -- rotation matrix, <m1|D|m2>, (m1,m2=j,j-1,...,-j) [O' = D * O].
    """
    alpha, beta, gamma = _axis_to_euler_zyz(theta, axis)
    rng = [sp.S(m) / 2 for m in range(2 * j, -2 * j - 1, -2)]
    pairs = [(m1, m2) for m1 in rng for m2 in rng]
    D = np.asarray([Rotation.D(j, m1, m2, alpha, beta, gamma).simplify() for m1, m2 in pairs])
    D = D.reshape(2 * j + 1, 2 * j + 1)

    return D


# ==================================================
def representation_matrix(l, n, axis, inversion=False, mirror=False, axial=False, tesseral=False):
    """
    Representation matrix of symmetry operation for spherical harmonics, l.

    Args:
        l (int): rank.
        n (int): n-fold rotation.
        axis (list): rotation or mirror axis (cartesian coordinate, no need for normalization).
        inversion (bool, optional): with inversion ?
        mirror (bool, optional): with mirror ?
        axial (bool, optional): for axial quantity ?
        tesseral (bool, optional): tesseral basis ?

    Returns:
        - (ndarray) -- representation matrix, <lm|g|lm'> or <alpha|g|beta> (tesseral=True) [alpha, beta=cl,sl,cl-1, ..., c0].
    """
    if mirror:
        n = 2
        inversion = True

    angle = 2 * sp.pi / n
    g = _wigner_D_matrix(l, angle, axis)

    p = l
    if axial:
        p += 1
    if inversion:
        g = (-1) ** p * g

    if tesseral:
        U = _create_unitary_matrix_tesseral(l)
        g = np.vectorize(sp.expand_complex)((U.conjugate().T) @ g @ U)

    return g


# ==================================================
def point_tag(tag):
    """
    Point-group part of symmetry operation tag.

    Args:
        tag (str): symmetry operation tag, e.g. "2[001]:[0,0,1/2]" or "2'[100]".

    Returns:
        - (str) -- tag without translation and time reversal, e.g. "2[001]".
    """
    return tag.split(":")[0].replace("'", "")


# ==================================================
@lru_cache(maxsize=None)
def _point_matrix(ptag, hexagonal, axial):
    d = TagSymmetryOperation.parse(ptag, hexagonal)
    g = representation_matrix(1, d["n"], str_to_sympy(d["axis"]), d["inversion"], d["mirror"], axial, tesseral=True)
    g = simplify_array(g)
    g.flags.writeable = False
    return g


# ==================================================
def symmetry_operation_matrix(tag, hexagonal=False, axial=False):
    """
    Symmetry operation matrix (cartesian coordinate, conventional cell).

    Args:
        tag (str): symmetry operation tag (time-reversal "'" is ignored).
        hexagonal (bool, optional): hexagonal systems ?
        axial (bool, optional): SO matrix for axial quantity ?

    Returns:
        - (ndarray) -- symmetry operation matrix in cartesian basis (3x3 or 4x4 with translation).
    """
    tag = tag.replace("'", "")
    g = _point_matrix(point_tag(tag), hexagonal, axial).copy()

    d = TagSymmetryOperation.parse(tag, hexagonal)
    if not d["point_group"]:
        t = str_to_sympy(d["partial_trans"])
        g4 = np.full((4, 4), sp.S(0))
        g4[0:3, 0:3] = g
        g4[0:3, 3] = t
        g4[3, 3] = sp.S(1)
        g = g4

    return g


# ==================================================
@lru_cache(maxsize=None)
def _multipole_matrix(s, ptag, hexagonal, axial, tesseral):
    if s == 0:
        d = TagSymmetryOperation.parse(ptag, False)
        g = sp.S(-1) if axial and (d["mirror"] or d["inversion"]) else sp.S(1)
        g = np.array(g, ndmin=2)
    else:
        d = TagSymmetryOperation.parse(ptag, hexagonal)
        g = representation_matrix(s, d["n"], str_to_sympy(d["axis"]), d["inversion"], d["mirror"], axial, tesseral)
        g = np.array(g, ndmin=2)
    g = simplify_array(g)
    g.flags.writeable = False
    return g


# ==================================================
def symmetry_operation_matrix_multipole(s, tag, hexagonal=False, axial=False, tesseral=False):
    """
    Symmetry operation matrix for multipole (lm or tesseral basis).

    Args:
        s (int): internal rank.
        tag (str): symmetry operation tag.
        hexagonal (bool, optional): hexagonal systems ?
        axial (bool, optional): SO matrix for axial quantity ?
        tesseral (bool, optional): tesseral basis ?

    Returns:
        - (ndarray) -- symmetry operation matrix in lm or tesseral basis, (2s+1 x 2s+1).

    Note:
        - translation and time reversal are ignored.
    """
    return _multipole_matrix(s, point_tag(tag), hexagonal, axial, tesseral).copy()


# ==================================================
def symmetry_operation_matrices(so, hexagonal):
    """
    Cartesian and fractional symmetry operation matrices of a list of operations.

    Args:
        so (list): symmetry operation tags.
        hexagonal (bool): hexagonal systems ?

    Returns:
        - (ndarray) -- cartesian SO matrices (3x3, or 4x4 with translation).
        - (ndarray) -- fractional SO matrices (3x3, or 4x4 with translation).
        - (ndarray) -- determinants.
    """
    cart = np.asarray([symmetry_operation_matrix(i, hexagonal) for i in so])
    frac = simplify_array(convert_to_fractional_hexagonal_matrix(cart)) if hexagonal else cart.copy()
    det = simplify_array([sp.Matrix(i[0:3, 0:3]).det() for i in cart])
    return cart, frac, det


# ==================================================
def group_symmetry_operation(tag):
    """
    Symmetry operations of a PG or SG from the group information tables (data_group_pg/sg).

    Args:
        tag (str): PG tag (e.g. "Oh", "D3-1", "C3-c") or SG tag (e.g. "Oh^1").

    Returns:
        - (dict) -- "tag" (list of SO tags), "cartesian" and "fractional" (ndarray of SO matrices, 4x4 for SG),
          "lattice" ("0" for PG, centering such as "P", "F" for SG).
    """
    info = sg_info[sg_id[tag]] if "^" in tag else pg_info[pg_id[tag]]
    cart, frac, _ = symmetry_operation_matrices(info.SO, info.hexagonal_g)
    return {"tag": info.SO, "cartesian": cart, "fractional": frac, "lattice": info.lattice}
