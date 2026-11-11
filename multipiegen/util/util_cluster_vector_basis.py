# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Orthonormal vector cluster basis (internal rank s=1) of Oh and D6h clusters (binary: cluster_vector_basis.pkl).

This is an intermediate basis to create the bond cluster basis (util_bond_cluster_basis, cluster_bond_basis.pkl); whether it is used
by itself is not decided yet.

Same procedure as the cluster basis (util_cluster_basis), for the vector harmonics. s=0 is the cluster basis itself.

- Vector harmonics. For a real top harmonics T (top key (l, irrep, mul, comp), u-matrix U) and k = -1, 0, 1,
  the vector harmonics is sum_m U[m] O_{lm}^{(s=1)}(k) (util_multipolar_harmonics, orbital rank l+k), a 3-vector
  whose components are in the order of the internal basis [px, py, pz] (data_basis.cubic_orbital_basis[1] and
  hexagonal_orbital_basis[1]).
- Type. For a polar internal basis (x="q"), it is polar (X="Q") if s+k is even and axial (X="G") if s+k is odd;
  for an axial internal basis (x="g"), Q and G are interchanged. A polar one has the irrep and component of T
  (coefficient 1); an axial one has the irrep, component and coefficient c of data_harmonics_mapping, G(key) = c T.
- Cluster vector. The values at the cluster points (N,3) (see util_cluster_basis.create_cluster).
- Orthonormalization. Within each (irrep, component), Gram-Schmidt in the order of rank l, then k = -1, 0, 1,
  then multiplicity; linearly dependent ones are dropped. The result is 3N = 144 (Oh) or 72 (D6h) vectors for each x.
  The linear dependence is first found numerically, then the kept vectors are orthonormalized exactly.

- cluster_vector_basis (dict)::

    { "Oh/D6h": {
        "site0": [x0,y0,z0] (sympy),
        "cluster": cluster points, ndarray (N,3) (sympy),
        "internal_basis": ["px","py","pz"],
        "basis": { PGMultipoleType: ndarray (dim,N,3) (sympy) },  # orthonormal, sum_{i,a} v_{ia}^2 = 1 per component.
        "dropped": [ PGMultipoleType ],  # vector harmonics linearly dependent on lower ones (whole multiplets).
    } }

    PGMultipoleType = (X, l, Gamma, n, p, s, k, x) as the vector SAMB of MultiPie Ver.2: type Q/G, rank, irrep,
        multiplicity (-1: unique), add_multiplicity p = -1 (multiplets are distinguished by k), s = 1, k = -1, 0, 1,
        internal type x = q/g. All components of a multiplet are in the value (rows, in the order of components).
"""

import numpy as np
from multipie import PGMultipoleType
import sympy as sp
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.data.data_basis import cubic_orbital_basis, hexagonal_orbital_basis
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.data.data_top_harmonics_info import top_harmonics_u_matrix
from multipiegen.util.util_cluster_basis import _dot, _simplify, create_cluster
from multipiegen.util.util_harmonics_symbol import is_real_irrep
from multipiegen.util.util_info import load_info
from multipiegen.util.util_multipolar_harmonics import InternalBasis, create_harmonics_set, get_internal_basis

# ==================================================
cluster_vector_basis_file = "cluster_vector_basis"
_internal_names = {"Oh": cubic_orbital_basis[1], "D6h": hexagonal_orbital_basis[1]}
_x, _y, _z = sp.symbols("x y z", real=True)
TOL_DEPENDENT = 1e-7  # relative residual below which a vector is regarded as linearly dependent.


# ==================================================
def internal_basis(group):
    """
    Internal basis of s=1 (the orbitals of data_basis), as InternalBasis (tesseral for other s).

    Args:
        group (str): "Oh" or "D6h".

    Returns:
        - (InternalBasis) -- internal basis.
    """
    bf = load_info()["harmonics"]["basis_function"]
    U1 = np.array([bf[n][1] for n in _internal_names[group]], dtype=object).T  # <m|alpha>.
    tesseral = get_internal_basis("tesseral")
    return InternalBasis(lambda s: U1 if s == 1 else tesseral.U(s), name=f"{group} p-orbitals")


# ==================================================
def _u_matrix(group, top_key):
    u = top_harmonics_u_matrix[group][top_key]
    return np.array(str_to_sympy(u) if isinstance(u, str) else u, dtype=object)


# ==================================================
def vector_harmonics(group):
    """
    Vector harmonics (s=1) of real top harmonics.

    Args:
        group (str): "Oh" or "D6h".

    Returns:
        - (dict) -- { (top_key, k): (3 sympy expressions in x,y,z) }, nonzero ones only.
    """
    basis = internal_basis(group)
    res = {}
    for tk in top_harmonics_u_matrix[group].keys():
        if not is_real_irrep(tk[1]):
            continue
        u = _u_matrix(group, tk)
        for k in (-1, 0, 1):
            if tk[0] + k < 0:
                continue
            v = tuple(sp.expand(sp.radsimp(e)) for e in u @ create_harmonics_set(tk[0], s=1, k=k, basis=basis))
            if any(e != 0 for e in v):
                res[(tk, k)] = v
    return res


# ==================================================
def label(group, top_key, k, x):
    """
    Key and coefficient of a vector harmonics.

    Args:
        group (str): "Oh" or "D6h".
        top_key (tuple): top key (rank, irrep, mul, comp).
        k (int): -1, 0, 1.
        x (str): internal basis type, "q" (polar) or "g" (axial).

    Returns:
        - (tuple) -- key (X, rank, irrep, mul, comp, k).
        - (sympy) -- coefficient c (the vector is c times the vector harmonics of top_key).
    """
    even = (1 + k) % 2 == 0
    X = ("Q" if even else "G") if x == "q" else ("G" if even else "Q")
    if X == "Q":
        return ("Q", *top_key, k), sp.S(1)
    for key, (c, tk) in data_harmonics_mapping[group]["G"].items():
        if tk == top_key:
            return ("G", *key, k), sp.S(c)
    raise KeyError(f"no axial harmonics for {top_key} in {group}.")


# ==================================================
def _order(key):
    """Gram-Schmidt order: rank, k (-1,0,1), multiplicity."""
    X, rank, irrep, mul, comp, k = key
    return (rank, k, mul)


# ==================================================
def _numerical_selection(vectors_num):
    """
    Numerical Gram-Schmidt within each (irrep, component).

    Args:
        vectors_num (dict): { key: float ndarray (3N,) }.

    Returns:
        - (list) -- kept keys.
        - (list) -- dropped keys.
        - (tuple) -- (min relative residual of kept, max relative residual of dropped).
    """
    groups = {}
    for key in vectors_num:
        groups.setdefault((key[2], key[4]), []).append(key)
    kept, dropped, rk, rd = [], [], [1.0], [0.0]
    for keys in groups.values():
        B = []
        for key in sorted(keys, key=_order):
            v = vectors_num[key].copy()
            v0 = np.linalg.norm(v)
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
def _gram_schmidt_exact(vectors):
    """
    Exact Gram-Schmidt within each (irrep, component) (all vectors must be independent).

    Args:
        vectors (dict): { key: object ndarray (3N,) (sympy) }.

    Returns:
        - (dict) -- orthonormal basis { key: ndarray (3N,) }.
    """
    groups = {}
    for key in vectors:
        groups.setdefault((key[2], key[4]), []).append(key)
    basis = {}
    for keys in groups.values():
        ortho = []
        for key in sorted(keys, key=_order):
            u = vectors[key].copy()
            for w, ww in ortho:
                c = _simplify(_dot(w, u) / ww)
                if c != 0:
                    u = np.array([_simplify(a - c * b) for a, b in zip(u, w)], dtype=object)
            uu = _dot(u, u)
            if uu == 0:
                raise ValueError(f"{key} is linearly dependent in exact arithmetic.")
            ortho.append((u, uu))
            n = sp.sqrt(uu)
            basis[key] = np.array([_simplify(a / n) for a in u], dtype=object)
    return basis


# ==================================================
def _multipole_keys(vectors, x):
    """
    Re-key by PGMultipoleType(X, l, Gamma, n, p=-1, s=1, k, x), with components stacked.

    Args:
        vectors (dict): { (X, rank, irrep, mul, comp, k): ndarray (N,3) or None }, complete multiplets.
        x (str): internal type, "q" or "g".

    Returns:
        - (dict) -- { PGMultipoleType: ndarray (dim,N,3) (or None) }, in the order of first appearance.
    """
    multiplets = {}
    for (X, rank, irrep, mul, comp, k), v in vectors.items():
        multiplets.setdefault((X, rank, irrep, mul, k), {})[comp] = v
    res = {}
    for (X, rank, irrep, mul, k), comps in multiplets.items():
        key = PGMultipoleType(X, rank, irrep, mul, -1, 1, k, x)
        vals = [comps[c] for c in sorted(comps)]
        res[key] = None if vals[0] is None else np.array(vals, dtype=object)
    return res


# ==================================================
def create_cluster_vector_basis(group, verbose=False):
    """
    Create orthonormal vector cluster basis (s=1) for polar (q) and axial (g) internal basis.

    Args:
        group (str): "Oh" or "D6h".
        verbose (bool, optional): print the separation of kept and dropped vectors ?

    Returns:
        - (dict) -- { "site0", "cluster", "internal_basis", "basis", "dropped" }.
    """
    site0, cluster = create_cluster(group)
    N = len(cluster)
    vh = vector_harmonics(group)
    P = np.array(cluster, dtype=float)

    # numerical values (shared by q and g, up to coefficient).
    num = {}
    for (tk, k), v in vh.items():
        f = sp.lambdify((_x, _y, _z), list(v), "numpy")
        num[(tk, k)] = np.array([np.broadcast_to(np.asarray(c, dtype=float), (N,)) for c in f(*P.T)]).T.ravel()

    exact_cache = {}

    def exact(tk, k):
        if (tk, k) not in exact_cache:
            v = vh[(tk, k)]
            exact_cache[(tk, k)] = np.array(
                [_simplify(e.xreplace({_x: p[0], _y: p[1], _z: p[2]})) for p in cluster for e in v], dtype=object
            )
        return exact_cache[(tk, k)]

    basis, dropped = {}, []
    for x in ["q", "g"]:
        lab = {(tk, k): label(group, tk, k, x) for (tk, k) in vh}
        vec_num = {lab[t][0]: float(lab[t][1]) * num[t] for t in vh}
        kept, drop, (rk, rd) = _numerical_selection(vec_num)
        if verbose:
            print(f"{group} {x}: kept {len(kept)}, dropped {len(drop)}, min residual kept {rk:.1e}, max dropped {rd:.1e}")
        if len(kept) != 3 * N:
            raise ValueError(f"number of basis ({len(kept)}) differs from 3N ({3 * N}).")
        inv = {lab[t][0]: t for t in vh}
        vec = {key: np.array([_simplify(lab[inv[key]][1] * e) for e in exact(*inv[key])], dtype=object) for key in kept}
        b = _gram_schmidt_exact(vec)
        order = [lab[t][0] for t in vh]
        basis.update(_multipole_keys({key: b[key].reshape(N, 3) for key in order if key in b}, x))
        dropped += list(_multipole_keys({key: None for key in order if key in drop}, x).keys())

    return {
        "site0": site0,
        "cluster": cluster,
        "internal_basis": list(_internal_names[group]),
        "basis": basis,
        "dropped": dropped,
    }


# ==================================================
def save_cluster_vector_basis(groups=("Oh", "D6h"), verbose=True):
    """
    Create and save vector cluster basis as binary (multipiegen/binary/cluster_vector_basis.pkl).

    Args:
        groups (tuple, optional): groups.
        verbose (bool, optional): verbose ?
    """
    bm = BinaryManager(topdir=__bin_dir__, verbose=verbose, comment=__doc__)
    for group in groups:
        bm[group] = create_cluster_vector_basis(group, verbose=verbose)
    bm.save_binary(cluster_vector_basis_file)


# ==================================================
def load_cluster_vector_basis():
    """
    Load vector cluster basis.

    Returns:
        - (dict) -- { "Oh/D6h": { "site0", "cluster", "internal_basis", "basis", "dropped" } }.
    """
    return dict(BinaryManager(cluster_vector_basis_file, topdir=__bin_dir__))


# ==================================================
def cluster_vector_basis_matrix(data, x, numeric=True):
    """
    Basis as matrix, each row is a component (N,3) flattened.

    Args:
        data (dict): vector cluster basis of a group, load_cluster_vector_basis()[group].
        x (str): "q" or "g".
        numeric (bool, optional): float ndarray ? otherwise sympy object ndarray.

    Returns:
        - (ndarray) -- basis matrix (3N, 3N).
        - (list) -- (PGMultipoleType, component index) of each row.
    """
    rows, keys = [], []
    for key, v in data["basis"].items():
        if key.x != x:
            continue
        for c, vc in enumerate(v):
            rows.append(vc.ravel())
            keys.append((key, c))
    m = np.array(rows, dtype=object)
    if numeric:
        m = m.astype(float)
    return m, keys
