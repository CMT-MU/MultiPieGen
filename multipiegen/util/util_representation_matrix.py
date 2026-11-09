# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Representation matrices of point groups from the harmonics (up to rank 11).

The harmonics f_c (component c) of a point group are those of data_harmonics_mapping, f_c = coeff * T (top harmonics).
A symmetry operation g (cartesian matrix R) acts as (g f)(r) = f(R^T r) for the polar (Q) harmonics, and with an
extra det(R) for the axial (G) harmonics. The representation matrix D(g) is defined by

    g f_c = sum_c' f_c' D_c'c(g),

so that D(g) = R for the dipole [x,y,z]. D(g) is obtained by least squares from the values of f_c at random points
on the unit sphere. The residual measures whether g f_c stays in the space of the components (closure).

- representation matrix (dict) = {"tag": [SO tag], "matrix": {irrep: ndarray(n_SO, dim, dim) (sympy)}}
    irreps in the order of the character table, matrices in the order of SO.
"""

from functools import lru_cache

import numpy as np
import sympy as sp
from multipie.util.util import str_to_sympy

from multipiegen.data.data_character_table import character
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.data.data_top_harmonics_info import top_harmonics_symmetry
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

TOL = 1e-9
N_POINT = 64
_x, _y, _z = sp.symbols("x y z", real=True)


# ==================================================
@lru_cache(maxsize=None)
def _top_function(top, key):
    """Numerical function of a top harmonics, f(x, y, z) for arrays."""
    e = str_to_sympy(top_harmonics_symmetry[top][key])
    e = e.subs({s: {"x": _x, "y": _y, "z": _z}[s.name] for s in e.free_symbols})
    f = sp.lambdify((_x, _y, _z), e, "numpy")
    return lambda p: np.broadcast_to(np.asarray(f(p[..., 0], p[..., 1], p[..., 2]), dtype=complex), p.shape[:-1])


# ==================================================
@lru_cache(maxsize=None)
def _points():
    rng = np.random.default_rng(1)
    p = rng.normal(size=(N_POINT, 3))
    return p / np.linalg.norm(p, axis=1)[:, None]


# ==================================================
def _group(pg):
    """(SO tags, cartesian matrices (float), det) of point group pg (tag)."""
    info = load_info()
    so = _group_data()[info["id"][pg]]["symmetry_operation"]
    R = np.asarray(so["cartesian"], dtype=float)
    det = np.asarray(so["det"], dtype=float)
    return list(so["tag"]), R, det


@lru_cache(maxsize=None)
def _group_data():
    return load_group()


# ==================================================
def harmonics_blocks(pg):
    """
    Harmonics of a point group grouped by (X, rank, irrep, multiplicity).

    Args:
        pg (str): point group tag.

    Returns:
        - (dict) -- {(X, rank, irrep, mul): [(coeff, top_key)]} with components in order.
    """
    m = data_harmonics_mapping[pg]
    blocks = {}
    for X in ["Q", "G"]:
        for key, (c, tk) in sorted(m[X].items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3])):
            blocks.setdefault((X, key[0], key[1], key[2]), []).append((key[3], c, tk))
    return {k: [(c, tk) for _, c, tk in sorted(v)] for k, v in blocks.items()}


# ==================================================
def numerical_representation_matrix(pg, block):
    """
    Representation matrices of a harmonics block (numerical).

    Args:
        pg (str): point group tag.
        block (tuple): (X, rank, irrep, mul).

    Returns:
        - (ndarray) -- D(g), (n_SO, dim, dim) complex.
        - (float) -- max residual of the least squares (closure, should be ~0).
    """
    top = data_harmonics_mapping[pg]["top"]
    comps = harmonics_blocks(pg)[block]
    _, R, det = _group(pg)
    p = _points()

    F = np.stack([complex(c) * _top_function(top, tk)(p) for c, tk in comps], axis=1)  # (N, dim)
    pg_ = np.einsum("ni,gij->gnj", p, R)  # R^T r as row vectors: r @ R.
    D, res = [], 0.0
    for g in range(len(R)):
        Fg = np.stack([complex(c) * _top_function(top, tk)(pg_[g]) for c, tk in comps], axis=1)
        if block[0] == "G":
            Fg = det[g] * Fg
        d, *_ = np.linalg.lstsq(F, Fg, rcond=None)
        res = max(res, float(np.abs(F @ d - Fg).max()))
        D.append(d)
    return np.asarray(D), res


# ==================================================
def check_representation_matrix(pg):
    """
    Check that all harmonics of the same irrep (polar and axial, all ranks and multiplicities) have the same
    representation matrices.

    Args:
        pg (str): point group tag.

    Returns:
        - (dict) -- {irrep: {"blocks": [(X, rank, mul)], "closure": max residual, "deviation": max abs(D - D_ref)}}.
    """
    by_irrep = {}
    for block in harmonics_blocks(pg):
        X, rank, irrep, mul = block
        D, res = numerical_representation_matrix(pg, block)
        d = by_irrep.setdefault(irrep, {"blocks": [], "closure": 0.0, "deviation": 0.0, "ref": D})
        d["blocks"].append((X, rank, mul))
        d["closure"] = max(d["closure"], res)
        d["deviation"] = max(d["deviation"], float(np.abs(D - d["ref"]).max()))
    for d in by_irrep.values():
        del d["ref"]
    return by_irrep


# ==================================================
def _exact(v):
    """Exact sympy value of a float (complex) number made of rationals and sqrt(3)."""
    parts = []
    for x in [v.real, v.imag]:
        e = sp.nsimplify(x, [sp.sqrt(3)], tolerance=1e-10) if abs(x) > TOL else sp.S(0)
        if abs(float(e) - x) > TOL:
            raise ValueError(f"cannot convert {x} to exact value.")
        parts.append(e)
    return sp.expand(parts[0] + sp.I * parts[1])


# ==================================================
def representation_matrix(pg):
    """
    Representation matrices of a point group (exact), from the lowest-rank harmonics of each irrep (polar first).

    Args:
        pg (str): point group tag.

    Returns:
        - (dict) -- {"tag": [SO tag], "matrix": {irrep: ndarray(n_SO, dim, dim) (sympy)}}.
    """
    tag, _, _ = _group(pg)
    blocks = harmonics_blocks(pg)
    ref = {}
    for block in sorted(blocks, key=lambda b: (b[1], b[0] != "Q", b[3])):  # rank, Q before G, multiplicity.
        ref.setdefault(block[2], block)

    matrix = {}
    for irrep in character[pg][1].keys():
        D, _ = numerical_representation_matrix(pg, ref[irrep])
        matrix[irrep] = np.vectorize(_exact, otypes=[object])(D)
    return {"tag": tag, "matrix": matrix}
