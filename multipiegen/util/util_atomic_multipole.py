# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Atomic multipoles (matrices in the atomic orbital basis) and their point-group SAMB.

Based on MultiPie Ver.2 (others/tools/create_binary/create_data_atomic_multipole.py and
create_data_atomic_multipole_group.py).

- Atomic multipole X_{lm}^{(s)}(k) (X = Q, G, T, M; s = 0, 1; l <= 7) in the spinful LMs or JML basis of
  data_basis (atomic_lms_basis, atomic_jm_basis; L = 0..3). The matrix elements are those of
  multipie.util.util_atomic_multipole (Wigner-Eckart); only the elements allowed by Jz conservation are evaluated.
  The matrix of each bra-ket block (L1, L2) (data_basis.atomic_rank_block) is kept if it is nonzero. "lm" is the
  spinless (s = 0, up-spin) part.
- Point-group SAMB. Every harmonics (X, l, Gamma, n) of a point group is c T with a top harmonics T of Oh or D6h
  (data_harmonics_mapping) whose u-matrix U_T (spherical m = l..-l) gives the top SAMB sum_m U_T[m] X_{lm}^{(s)}(k)
  for X = Q, T (polar harmonics) or G, M (axial harmonics). Each component is normalized to the Frobenius norm 1
  (diagonal block) or 1/sqrt(2) (off-diagonal block), which is the Gram-Schmidt of MultiPie Ver.2, since the components are
  orthogonal. The SAMB of the point group is c times the top SAMB (key replaced), computed once per top harmonics.
  "jml" is computed from the JML matrices; "lgs" is U_bra^dagger M U_ket from the LMs matrices M, with U the orbitals
  of data_basis.atomic_lgs_basis (atomic_orbital_basis, basis_function u-matrices, spin expanded); "lg" is the
  spinless (s = 0, up-spin) part of "lgs" multiplied by sqrt(2).
- Key. PGMultipoleType(X, l, Gamma, n, p, s, k, x="q"): n is the multiplicity of the harmonics (X, l, Gamma) (its xyz
  polynomial, not renumbered unlike MultiPie Ver.2); p = -1 (no duplicates).
- Order. X (Q, G, T, M), s, k, l, Gamma (character table), n.

- atomic_multipole (dict)::

    { "lms"/"jml"/"lm": { (L1, L2): { SphericalMultipoleType(X, l, s, k, "q"): ndarray(2l+1, n1, n2) } } }

- atomic_multipole_group (dict)::

    { PG_id: { "jml"/"lgs"/"lg": { (L1, L2): { PGMultipoleType: (ndarray(dim, n1, n2), ndarray(dim)) } } } }
        matrix of each component, cartesian expression of the harmonics.
"""

from functools import lru_cache

import numpy as np
import sympy as sp
from multipie import PGMultipoleType, SphericalMultipoleType
from multipie.util.util_atomic_multipole import (
    Pl_sk_J1J2L1L2,
    Xl_l1l2,
    Xlm_lm,
    is_zero,
    pauli_matrix_n,
    qmtg_map_spinful,
)
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager
from sympy.physics.wigner import wigner_3j

from multipiegen import __bin_dir__
from multipiegen.data.data_basis import (
    atomic_jm_basis,
    atomic_lms_basis,
    atomic_orbital_basis,
    atomic_rank_block,
)
from multipiegen.data.data_character_table import character
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.data.data_top_harmonics_info import top_harmonics_u_matrix
from multipiegen.util.util_info import load_info
from multipiegen.util.util_symmetry_operation import simplify_array, simplify_value

# ==================================================
atomic_multipole_file = "atomic_multipole"
atomic_multipole_group_file = "atomic_multipole_group"
max_rank = 7
max_s_rank = 1
heads = ["Q", "G", "T", "M"]
_half = sp.Rational(1, 2)


# ==================================================
def _parse_state(state, L):
    """'(m,u)' (lms) or '(J,M)' (jml) of data_basis -> (L, m, s) or (J, M, L) in sympy."""
    a, b = state.strip("()").split(",")
    if b in ["u", "d"]:
        return (sp.S(L), sp.S(a), _half if b == "u" else -_half)
    return (sp.S(a), sp.S(b), sp.S(L))


# ==================================================
@lru_cache(maxsize=None)
def spinful_basis(b_type, L):
    """
    Spinful basis of rank L (data_basis.atomic_lms_basis or atomic_jm_basis).

    Args:
        b_type (str): "lms" [(L, m, s)] or "jml" [(J, M, L)].
        L (int): orbital rank.

    Returns:
        - (tuple) -- states as tuples of sympy.
    """
    basis = atomic_lms_basis if b_type == "lms" else atomic_jm_basis
    return tuple(_parse_state(i, L) for i in basis[L])


# ==================================================
def _element(b_type, X, l, m, s, k, b1, b2):
    """Matrix element <b1|X_{lm}^{(s)}(k)|b2> (0 when Jz is not conserved)."""
    if b_type == "lms":
        (l1, m1, s1), (l2, m2, s2) = b1, b2
        if m1 + s1 != m2 + s2 + m:
            return 0
        if s == 0:
            return Xlm_lm(X, l, m, (l1, m1), (l2, m2)) if s1 == s2 else 0
        Y = qmtg_map_spinful(X, s, k)
        n = s1 - s2
        xlm = Xlm_lm(Y, l + k, m - n, (l1, m1), (l2, m2))
        if xlm == 0:
            return 0
        v = wigner_3j(l + k, l, s, m - n, -m, n) * xlm * pauli_matrix_n(s, n, s1, s2)
        return v * sp.I ** (s + k) * (-1) ** (l + m) * sp.sqrt(2 * l + 1)
    else:
        (J1, M1, l1), (J2, M2, l2) = b1, b2
        if M1 != M2 + m:
            return 0
        Y = qmtg_map_spinful(X, s, k)
        if is_zero(Y, l + k, l1, l2) or not abs(l1 - l2) <= l + k <= l1 + l2:
            return 0
        w = wigner_3j(J1, J2, l, -M1, M2, m)
        if w == 0:
            return 0
        P = Pl_sk_J1J2L1L2(l, s, k, J1, J2, l1, l2)
        if P == 0:
            return 0
        return (-1) ** (J1 + M1) * w * P * Xl_l1l2(Y, l + k, l1, l2)


# ==================================================
def atomic_multipole_matrix(b_type, X, l, s, k):
    """
    Matrices of an atomic multipole in each bra-ket block.

    Args:
        b_type (str): "lms" or "jml".
        X (str): head, "Q", "G", "T" or "M".
        l (int): rank.
        s (int): internal rank, 0 or 1.
        k (int): -s <= k <= s.

    Returns:
        - (dict) -- {(L1, L2): ndarray(2l+1, n1, n2)} of the nonzero blocks.
    """
    Y = qmtg_map_spinful(X, s, k)
    res = {}
    for L1, L2 in atomic_rank_block:
        if is_zero(Y, l + k, L1, L2) or not abs(L1 - L2) <= l + k <= L1 + L2 or l + k < 0:
            continue
        r1, r2 = spinful_basis(b_type, L1), spinful_basis(b_type, L2)
        mat = np.full((2 * l + 1, len(r1), len(r2)), sp.S(0), dtype=object)
        for im, m in enumerate(range(l, -l - 1, -1)):
            for i, b1 in enumerate(r1):
                for j, b2 in enumerate(r2):
                    v = _element(b_type, X, l, m, s, k, b1, b2)
                    if v != 0:
                        mat[im, i, j] = simplify_value(v)
        if not (mat == 0).all():
            res[(L1, L2)] = mat
    return res


# ==================================================
def atomic_multipole_keys():
    """List of (X, l, s, k) of atomic multipoles."""
    return [(X, l, s, k) for X in heads for s in range(max_s_rank + 1) for k in range(-s, s + 1) for l in range(max_rank + 1)]


# ==================================================
def assemble_atomic_multipole(results):
    """
    Assemble atomic multipole data.

    Args:
        results (dict): {(b_type, X, l, s, k): {(L1, L2): matrix}} (atomic_multipole_matrix).

    Returns:
        - (dict) -- {"lms"/"jml"/"lm": {(L1, L2): {SphericalMultipoleType: matrix}}}.
    """
    data = {}
    for bt in ["lms", "jml"]:
        data[bt] = {}
        for block in atomic_rank_block:
            d = {}
            for X, l, s, k in atomic_multipole_keys():  # order Q/G/T/M, s, k, l.
                mat = results[(bt, X, l, s, k)].get(block)
                if mat is not None:
                    d[SphericalMultipoleType(X, l, s, k, "q")] = mat
            data[bt][block] = d
    data["lm"] = {block: {key: m[:, ::2, ::2] for key, m in d.items() if key.s == 0} for block, d in data["lms"].items()}
    return data


# ==================================================
def load_atomic_multipole():
    """
    Load atomic multipole (multipiegen/binary/atomic_multipole.pkl).

    Returns:
        - (BinaryManager) -- {"lms"/"jml"/"lm": {(L1, L2): {SphericalMultipoleType: ndarray(2l+1, n1, n2)}}}.
    """
    return BinaryManager(atomic_multipole_file, topdir=__bin_dir__)


# ==================================================
def load_atomic_multipole_group():
    """
    Load atomic multipole of point groups (multipiegen/binary/atomic_multipole_group.pkl).

    Returns:
        - (BinaryManager) -- {PG_id: {"jml"/"lgs"/"lg": {(L1, L2): {PGMultipoleType: (matrix, ex)}}}}.
    """
    return BinaryManager(atomic_multipole_group_file, topdir=__bin_dir__)


# ==================================================
@lru_cache(maxsize=None)
def lgs_u_matrix():
    """
    U-matrix from the LMs basis (data_basis.atomic_lms_basis) to the Lgs basis (data_basis.atomic_lgs_basis).

    Returns:
        - (dict) -- {L: ndarray(2(2L+1), 2(2L+1))}, ``|L gamma s> = sum |L m s> U``.
    """
    bf = load_info()["harmonics"]["basis_function"]
    U = {}
    for L, names in enumerate(atomic_orbital_basis):
        U1 = np.array([bf[n][1] for n in names], dtype=object).T
        Uf = np.full((2 * U1.shape[0], 2 * U1.shape[1]), sp.S(0), dtype=object)
        Uf[::2, ::2] = U1
        Uf[1::2, 1::2] = U1
        U[L] = Uf
    return U


# ==================================================
def _u(u):
    return np.array(str_to_sympy(u) if isinstance(u, str) else u, dtype=object)


# ==================================================
def top_samb(top, top_key, mats, diagonal):
    """
    Normalized top SAMB, sum_m U_T[m] X_{lm} for a component of the top harmonics.

    Args:
        top (str): "Oh" or "D6h".
        top_key (tuple): key of the top harmonics component.
        mats (ndarray): (2l+1, n1, n2) atomic multipole matrices.
        diagonal (bool): diagonal block ?

    Returns:
        - (ndarray) -- (n1, n2) matrix with the Frobenius norm 1 (diagonal) or 1/sqrt(2) (off-diagonal).
    """
    u = _u(top_harmonics_u_matrix[top][top_key])
    m = simplify_array(np.tensordot(u, mats, axes=(0, 0)))
    norm2 = simplify_value(sum(sp.expand(v * sp.conjugate(v)) for v in m.ravel() if v != 0))
    f = sp.sqrt(norm2) if diagonal else sp.sqrt(2 * norm2)
    return simplify_array(m / f)


# ==================================================
def _harmonics_head(X):
    return ["Q", "T"] if X == "Q" else ["G", "M"]


# ==================================================
def top_keys_used():
    """Set of (top, top_key, X) used by point group harmonics up to max_rank."""
    used = set()
    for m in data_harmonics_mapping.values():
        for X in ["Q", "G"]:
            for (l, _, _, _), (_, tk) in m[X].items():
                if l <= max_rank:
                    used.add((m["top"], tk, X))
    return used


# ==================================================
def create_top_samb(atomic_multipole, top, block, b_type):
    """
    Top SAMB of a bra-ket block.

    Args:
        atomic_multipole (dict): atomic multipole data.
        top (str): "Oh" or "D6h".
        block (tuple): (L1, L2).
        b_type (str): "jml" or "lgs".

    Returns:
        - (dict) -- {(top_key, head, s, k): ndarray(n1, n2)}.
    """
    src = atomic_multipole["jml" if b_type == "jml" else "lms"][block]
    bra, ket = block
    if b_type == "lgs":
        U = lgs_u_matrix()
        Ubra, Uket = U[bra].conjugate().T, U[ket]
    res = {}
    for t, tk, X in sorted(top_keys_used(), key=str):
        if t != top:
            continue
        l = tk[0]
        for head in _harmonics_head(X):
            for key, mats in src.items():
                if key.X != head or key.l != l:
                    continue
                m = top_samb(top, tk, mats, bra == ket)
                if b_type == "lgs":
                    m = simplify_array(Ubra @ m @ Uket)
                res[(tk, head, key.s, key.k)] = m
    return res


# ==================================================
def create_atomic_multipole_pg(pg, top_data, ex):
    """
    Atomic multipole SAMB of a point group.

    Args:
        pg (str): point group tag.
        top_data (dict): {b_type: {block: {(top_key, head, s, k): matrix}}} of the top group of pg.
        ex (dict): {(X, l, Gamma, n): ex ndarray(dim)} cartesian expressions of the harmonics (pg_harmonics).

    Returns:
        - (dict) -- {"jml"/"lgs"/"lg": {(L1, L2): {PGMultipoleType: (matrix ndarray(dim, n1, n2), ex ndarray(dim))}}}.
    """
    mapping = data_harmonics_mapping[pg]
    irreps = list(character[pg][1].keys())
    multiplets = {}
    for X in ["Q", "G"]:
        for (l, irrep, n, comp), (c, tk) in mapping[X].items():
            if l <= max_rank:
                multiplets.setdefault((X, l, irrep, n), []).append((comp, sp.S(c), tk))
    for v in multiplets.values():
        v.sort(key=lambda t: t[0])

    order = lambda key: (heads.index(key.X), key.s, key.k, key.l, irreps.index(key.Gamma), key.n, key.p)
    res = {}
    for bt in ["jml", "lgs"]:
        res[bt] = {}
        for block in atomic_rank_block:
            top = top_data[bt][block]
            d = {}
            for (X, l, irrep, n), comps in multiplets.items():
                for head in _harmonics_head(X):
                    for s in range(max_s_rank + 1):
                        for k in range(-s, s + 1):
                            if (comps[0][2], head, s, k) not in top:
                                continue
                            mat = np.array([c * top[(tk, head, s, k)] for _, c, tk in comps], dtype=object)
                            d[PGMultipoleType(head, l, irrep, n, -1, s, k, "q")] = (mat, ex[(X, l, irrep, n)])
            res[bt][block] = {key: d[key] for key in sorted(d, key=order)}
    res["lg"] = {
        block: {key: (mat[:, ::2, ::2] * sp.sqrt(2), e) for key, (mat, e) in d.items() if key.s == 0}
        for block, d in res["lgs"].items()
    }
    return res
