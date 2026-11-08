# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Orthonormal cluster basis of Oh and D6h clusters.

The cluster is made by substituting the representative site (site0) into the general points
of the Oh (48 points) or D6h (24 points) cluster. The value of each real top harmonics at
each cluster point gives a cluster vector. Different irreps (and components) are orthogonal
by symmetry. Within each (irrep, component), the vectors are orthogonalized by the
Gram-Schmidt process from lower rank (then lower multiplicity), and linearly dependent ones
are dropped. The result is an orthonormal basis of 48 (Oh) or 24 (D6h) vectors.

- cluster_basis (dict).
    { "Oh/D6h": {
        "site0": [x0,y0,z0] (sympy),
        "cluster": cluster points, ndarray (N,3) (sympy),
        "basis": { top_key: ndarray (N,) (sympy) },  # orthonormal, sum_i |v_i|^2 = 1.
        "dropped": [top_key],  # harmonics linearly dependent on lower ones.
    } }
    top_key = (rank, irrep, multiplicity [1-based], component [0-based]), -1: unique.
"""

import numpy as np
import sympy as sp
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager
from multipiegen import __bin_dir__
from multipiegen.data.data_root_cluster import (
    oh_cluster_cart,
    d6h_cluster_cart,
    oh_cluster_cart_site0,
    d6h_cluster_cart_site0,
)
from multipiegen.data.data_top_harmonics_info import top_harmonics_symmetry
from multipiegen.util.util_harmonics_symbol import is_real_irrep

# ==================================================
cluster_basis_file = "cluster_basis"
_cluster_data = {
    "Oh": (oh_cluster_cart, oh_cluster_cart_site0),
    "D6h": (d6h_cluster_cart, d6h_cluster_cart_site0),
}
_xyz = sp.symbols("x y z")


# ==================================================
def _simplify(e):
    return sp.radsimp(sp.expand(e))


# ==================================================
def _dot(u, v):
    return _simplify(sum(a * b for a, b in zip(u, v)))


# ==================================================
def create_cluster(group):
    """
    Create cluster by substituting site0 into general points.

    Args:
        group (str): "Oh" or "D6h".

    Returns:
        - (ndarray) -- site0 (3,) (sympy).
        - (ndarray) -- cluster points (N,3) (sympy).
    """
    cluster, site0 = _cluster_data[group]
    site0 = str_to_sympy(site0)
    sub = dict(zip(_xyz, site0))
    cluster = np.vectorize(lambda e: _simplify(sp.sympify(e).subs(sub)), otypes=[object])(str_to_sympy(cluster))
    return site0, cluster


# ==================================================
def create_harmonics_vectors(group, cluster):
    """
    Values of real top harmonics at cluster points.

    Args:
        group (str): "Oh" or "D6h".
        cluster (ndarray): cluster points (N,3).

    Returns:
        - (dict) -- { top_key: ndarray (N,) (sympy) }.
    """
    vec = {}
    for key, expr in top_harmonics_symmetry[group].items():
        if not is_real_irrep(key[1]):
            continue
        f = str_to_sympy(expr)
        vec[key] = np.array([_simplify(f.subs(dict(zip(_xyz, p)))) for p in cluster], dtype=object)
    return vec


# ==================================================
def gram_schmidt(vectors):
    """
    Gram-Schmidt orthonormalization within each (irrep, component), from lower rank.

    Args:
        vectors (dict): { top_key: ndarray (N,) }.

    Returns:
        - (dict) -- orthonormal basis { top_key: ndarray (N,) }, in the order of vectors.
        - (list) -- dropped (linearly dependent) top_keys.
    """
    chains = {}
    for key in vectors.keys():
        chains.setdefault((key[1], key[3]), []).append(key)

    basis = {}
    dropped = []
    for keys in chains.values():
        keys = sorted(keys, key=lambda k: (k[0], k[2]))
        ortho = []  # [(u, <u|u>)], orthogonal but not normalized (keeps expressions simple).
        for key in keys:
            u = vectors[key].copy()
            for w, ww in ortho:
                c = _simplify(_dot(w, u) / ww)
                if c != 0:
                    u = np.array([_simplify(a - c * b) for a, b in zip(u, w)], dtype=object)
            uu = _dot(u, u)
            if uu == 0:
                dropped.append(key)
                continue
            ortho.append((u, uu))
            n = sp.sqrt(uu)
            basis[key] = np.array([_simplify(a / n) for a in u], dtype=object)

    order = list(vectors.keys())
    basis = {k: basis[k] for k in order if k in basis}
    dropped = [k for k in order if k in dropped]
    return basis, dropped


# ==================================================
def _check_components(basis, dropped):
    """
    All components of a multiplet must be kept or dropped together.
    """
    kept = {}
    for key in basis.keys():
        kept.setdefault((key[0], key[1], key[2]), set()).add(key[3])
    for key in dropped:
        if (key[0], key[1], key[2]) in kept:
            raise ValueError(f"components of {key[:3]} are partially dropped.")


# ==================================================
def create_cluster_basis(group):
    """
    Create orthonormal cluster basis.

    Args:
        group (str): "Oh" or "D6h".

    Returns:
        - (dict) -- { "site0", "cluster", "basis", "dropped" }.
    """
    site0, cluster = create_cluster(group)
    vectors = create_harmonics_vectors(group, cluster)
    basis, dropped = gram_schmidt(vectors)
    _check_components(basis, dropped)
    if len(basis) != len(cluster):
        raise ValueError(f"number of basis ({len(basis)}) differs from number of cluster points ({len(cluster)}).")
    return {"site0": site0, "cluster": cluster, "basis": basis, "dropped": dropped}


# ==================================================
def save_cluster_basis(groups=("Oh", "D6h"), verbose=True):
    """
    Create and save cluster basis as binary (multipiegen/binary_data/cluster_basis.pkl).

    Args:
        groups (tuple, optional): groups.
        verbose (bool, optional): verbose ?
    """
    bm = BinaryManager(topdir=__bin_dir__, verbose=verbose, comment=__doc__)
    for group in groups:
        bm[group] = create_cluster_basis(group)
    bm.save_binary(cluster_basis_file)


# ==================================================
def load_cluster_basis():
    """
    Load cluster basis.

    Returns:
        - (dict) -- { "Oh/D6h": { "site0", "cluster", "basis", "dropped" } }.
    """
    return dict(BinaryManager(cluster_basis_file, topdir=__bin_dir__))


# ==================================================
def cluster_basis_matrix(data, numeric=True):
    """
    Basis as matrix, each row is a basis vector.

    Args:
        data (dict): cluster basis of a group, load_cluster_basis()[group].
        numeric (bool, optional): float ndarray ? otherwise sympy object ndarray.

    Returns:
        - (ndarray) -- basis matrix (n_basis, N).
        - (list) -- top_key of each row.
    """
    keys = list(data["basis"].keys())
    m = np.array([data["basis"][k] for k in keys], dtype=object)
    if numeric:
        m = m.astype(float)
    return m, keys
