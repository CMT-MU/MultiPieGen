# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Orthonormal bond cluster basis of Oh and D6h clusters.

- Bond vectors. b_i is the general point i of the Oh or D6h cluster with x,y,z replaced by vector0 (unit vector,
  info["root_cluster"]["oh/d6h_vector0_cart"]).
- Bond cluster vectors. For each vector cluster basis v_key (util_cluster_vector_basis, (N,3)), the N-component
  vector f_key(i) = v_key(i) . b_i, with the same key. They are 3N vectors in the N-dimensional space.
- Orthonormalization. Within each (irrep, component), Gram-Schmidt in the order of rank, then k = -1, 0, 1, then
  multiplicity; linearly dependent ones (including vanishing ones) are dropped. The result is N = 48 (Oh) or
  24 (D6h) vectors for each internal basis type x = q, g. The linear dependence is first found numerically,
  then the kept vectors are orthonormalized exactly.
- Symmetry. f_key transforms as a scalar function on the cluster, (g f)(i) = s f(j) with R_g^{-1} p_i = p_j,
  s = 1 (x = q) or det R_g (x = g), with the representation matrices of its irrep.

- bond_cluster_basis (dict)::

    { "Oh/D6h": {
        "site0": [x0,y0,z0] (sympy),
        "vector0": [X0,Y0,Z0] (sympy),
        "cluster": cluster points, ndarray (N,3) (sympy),
        "bond": bond vectors b_i, ndarray (N,3) (sympy),
        "basis": { "q/g": { key: ndarray (N,) (sympy) } },  # orthonormal.
        "dropped": { "q/g": [key] },  # linearly dependent on lower ones (or vanishing).
    } }

    key = (X, rank, irrep, multiplicity [1-based], component [0-based], k) of the vector cluster basis.
"""

import numpy as np
import sympy as sp
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.data.data_root_cluster import d6h_cluster_cart, oh_cluster_cart
from multipiegen.util.util_cluster_basis import _simplify
from multipiegen.util.util_cluster_vector_basis import (
    TOL_DEPENDENT,
    _gram_schmidt_exact,
    _order,
    load_cluster_vector_basis,
)
from multipiegen.util.util_info import load_info

# ==================================================
bond_cluster_basis_file = "bond_cluster_basis"
_cluster_cart = {"Oh": oh_cluster_cart, "D6h": d6h_cluster_cart}
_vector0 = {"Oh": "oh_vector0_cart", "D6h": "d6h_vector0_cart"}
_xyz = sp.symbols("x y z")


# ==================================================
def create_bond(group):
    """
    Bond vectors by substituting vector0 into general points.

    Args:
        group (str): "Oh" or "D6h".

    Returns:
        - (ndarray) -- vector0 (3,) (sympy).
        - (ndarray) -- bond vectors (N,3) (sympy).
    """
    vector0 = load_info()["root_cluster"][_vector0[group]]
    sub = dict(zip(_xyz, vector0))
    bond = np.vectorize(lambda e: _simplify(sp.sympify(e).subs(sub)), otypes=[object])(str_to_sympy(_cluster_cart[group]))
    return vector0, bond


# ==================================================
def _numerical_selection(vectors_num):
    """
    Numerical Gram-Schmidt within each (irrep, component).

    Args:
        vectors_num (dict): { key: float ndarray (N,) }.

    Returns:
        - (list) -- kept keys.
        - (list) -- dropped keys.
        - (tuple) -- (min relative residual of kept, max relative residual of dropped (vanishing ones as 0)).
    """
    scale = max(np.linalg.norm(v) for v in vectors_num.values())
    groups = {}
    for key in vectors_num:
        groups.setdefault((key[2], key[4]), []).append(key)
    kept, dropped, rk, rd = [], [], [1.0], [0.0]
    for keys in groups.values():
        B = []
        for key in sorted(keys, key=_order):
            v = vectors_num[key].copy()
            v0 = np.linalg.norm(v)
            if v0 < TOL_DEPENDENT * scale:  # vanishing.
                dropped.append(key)
                continue
            for b in B:
                v -= (b @ v) * b
            r = np.linalg.norm(v) / v0
            if r > TOL_DEPENDENT:
                B.append(v / np.linalg.norm(v))
                kept.append(key)
                rk.append(r)
            else:
                dropped.append(key)
                rd.append(r)
    return kept, dropped, (min(rk), max(rd))


# ==================================================
def create_bond_cluster_basis(group, verbose=False):
    """
    Create orthonormal bond cluster basis from the vector cluster basis.

    Args:
        group (str): "Oh" or "D6h".
        verbose (bool, optional): print the separation of kept and dropped vectors ?

    Returns:
        - (dict) -- { "site0", "vector0", "cluster", "bond", "basis", "dropped" }.
    """
    vb = load_cluster_vector_basis()[group]
    vector0, bond = create_bond(group)
    N = len(bond)
    bond_f = np.array(bond, dtype=float)

    basis, dropped = {}, {}
    for x in ["q", "g"]:
        vbx = vb["basis"][x]
        num = {k: np.einsum("ia,ia->i", np.array(v, dtype=float), bond_f) for k, v in vbx.items()}
        kept, drop, (rk, rd) = _numerical_selection(num)
        if verbose:
            print(f"{group} {x}: kept {len(kept)}, dropped {len(drop)}, min residual kept {rk:.1e}, max dropped {rd:.1e}")
        if len(kept) != N:
            raise ValueError(f"number of basis ({len(kept)}) differs from N ({N}).")
        vec = {
            k: np.array([_simplify(sum(a * b for a, b in zip(vbx[k][i], bond[i]))) for i in range(N)], dtype=object) for k in kept
        }
        b = _gram_schmidt_exact(vec)
        basis[x] = {k: b[k] for k in vbx if k in b}
        dropped[x] = [k for k in vbx if k in drop]

    return {
        "site0": vb["site0"],
        "vector0": vector0,
        "cluster": vb["cluster"],
        "bond": bond,
        "basis": basis,
        "dropped": dropped,
    }


# ==================================================
def save_bond_cluster_basis(groups=("Oh", "D6h"), verbose=True):
    """
    Create and save bond cluster basis as binary (multipiegen/binary/bond_cluster_basis.pkl).

    Args:
        groups (tuple, optional): groups.
        verbose (bool, optional): verbose ?
    """
    bm = BinaryManager(topdir=__bin_dir__, verbose=verbose, comment=__doc__)
    for group in groups:
        bm[group] = create_bond_cluster_basis(group, verbose=verbose)
    bm.save_binary(bond_cluster_basis_file)


# ==================================================
def load_bond_cluster_basis():
    """
    Load bond cluster basis.

    Returns:
        - (dict) -- { "Oh/D6h": { "site0", "vector0", "cluster", "bond", "basis", "dropped" } }.
    """
    return dict(BinaryManager(bond_cluster_basis_file, topdir=__bin_dir__))


# ==================================================
def bond_cluster_basis_matrix(data, x, numeric=True):
    """
    Basis as matrix, each row is a basis vector.

    Args:
        data (dict): bond cluster basis of a group, load_bond_cluster_basis()[group].
        x (str): "q" or "g".
        numeric (bool, optional): float ndarray ? otherwise sympy object ndarray.

    Returns:
        - (ndarray) -- basis matrix (N, N).
        - (list) -- key of each row.
    """
    keys = list(data["basis"][x].keys())
    m = np.array([data["basis"][x][k] for k in keys], dtype=object)
    if numeric:
        m = m.astype(float)
    return m, keys
