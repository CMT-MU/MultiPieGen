# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Multipolar harmonics (internal rank s = 1, 2, 3) of point groups with polar (q) and axial (g) internal basis.

Based on MultiPie Ver.2 (others/tools/create_binary/create_data_multipolar_harmonics.py).

- Harmonics of a point group. For each harmonics (X, l, Gamma, n) of the point group (data_harmonics_mapping,
  X(key) = c T with the top harmonics T), its cartesian expression is c T and its u-matrix is U = c U_T (m, dim).
- Internal coefficient. For s = 1, 2, 3 and k = -s..s, internal_coefficient = U^T O_l^{(s)}(k) (dim, 2s+1): the
  coefficient (polynomial in x, y, z, each element factorized by sympy.factor) of each internal basis function, where
  O_l^{(s)}(k) is the multipolar spherical harmonics (util_multipolar_harmonics, orbital rank l+k) whose internal
  components are in the order of the internal basis: the orbitals of data_basis.cubic_orbital_basis[s] (cubic
  groups, top Oh) or hexagonal_orbital_basis[s] (trigonal and hexagonal groups, top D6h).
- Internal type. x = "q" (polar internal basis) if (X = Q and s+k even) or (X = G and s+k odd), otherwise "g".
- Multipole expression. multipole_expression = internal_coefficient @ (internal symbols) (dim), where the internal
  symbols are those of info["harmonics"]["internal_basis"] ((s, x) -> [Q or G + orbital_latex_name of each orbital]).
- Key. PGMultipoleType(X, l, Gamma, n, p, s, k, x): n is the multiplicity of the harmonics (X, l, Gamma, n), which
  fixes its xyz polynomial (not renumbered, unlike MultiPie Ver.2); p numbers those with the same
  (x, l, X, Gamma, n, k) as 1, 2, ... (-1 if unique). Vanishing ones (w3j = 0) are excluded.
- Order. x (q, g), s, Gamma (character table), l, X (Q, G), k, n, p.

- multipolar_harmonics (dict)::

    { PG_id: { PGMultipoleType: (internal_coefficient, multipole_expression) } }
        internal_coefficient: ndarray (dim, 2s+1) (sympy), coefficient (polynomial in x, y, z) of each internal basis
            function, for each component.
        multipole_expression: ndarray (dim) (sympy), multipole expression with the internal symbols.
        The expression of the harmonics (X, l, Gamma, n) itself is in harmonics.pkl (util_harmonics).
"""

from functools import lru_cache

import numpy as np
import sympy as sp
from multipie import PGMultipoleType
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.data.data_basis import cubic_orbital_basis, hexagonal_orbital_basis
from multipiegen.data.data_character_table import character
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.data.data_top_harmonics_info import top_harmonics_symmetry, top_harmonics_u_matrix
from multipiegen.util.util_info import load_info
from multipiegen.util.util_multipolar_harmonics import InternalBasis, create_harmonics_set, get_internal_basis

# ==================================================
multipolar_harmonics_file = "multipolar_harmonics"
max_rank = 3
max_s_rank = 3
_x, _y, _z = sp.symbols("x y z", real=True)
_crystal = {"Oh": "cubic", "D6h": "hexagonal"}
_orbital_basis = {"cubic": cubic_orbital_basis, "hexagonal": hexagonal_orbital_basis}


# ==================================================
@lru_cache(maxsize=None)
def internal_basis(crystal):
    """
    Internal basis of s = 1, 2, 3: the orbitals of data_basis (u-matrices of info.pkl).

    Args:
        crystal (str): "cubic" or "hexagonal".

    Returns:
        - (InternalBasis) -- internal basis.
    """
    bf = load_info()["harmonics"]["basis_function"]
    U = {s: np.array([bf[n][1] for n in _orbital_basis[crystal][s]], dtype=object).T for s in range(1, max_s_rank + 1)}
    tesseral = get_internal_basis("tesseral")
    return InternalBasis(lambda s: U[s] if s in U else tesseral.U(s), name=f"{crystal} orbitals")


# ==================================================
@lru_cache(maxsize=None)
def internal_symbols(crystal):
    """
    Internal symbols, {(s, x): [sympy]} (info["harmonics"]["internal_basis"]).
    """
    return dict(load_info()["harmonics"]["internal_basis"][crystal])


# ==================================================
def _expr(e):
    e = str_to_sympy(e) if isinstance(e, str) else sp.S(e)
    return sp.expand(e.subs({s: {"x": _x, "y": _y, "z": _z}[s.name] for s in e.free_symbols}))


# ==================================================
def pg_harmonics(pg, rank=max_rank):
    """
    Harmonics of a point group up to rank, as multiplets.

    Args:
        pg (str): point group tag.
        rank (int, optional): max rank.

    Returns:
        - (dict) -- {(X, l, irrep, n): (ex ndarray (dim), U ndarray (2l+1, dim))}, components in order.
    """
    m = data_harmonics_mapping[pg]
    top = m["top"]
    blocks = {}
    for X in ["Q", "G"]:
        for (l, irrep, n, comp), (c, tk) in m[X].items():
            if l > rank:
                continue
            u = top_harmonics_u_matrix[top][tk]
            u = np.array(str_to_sympy(u) if isinstance(u, str) else u, dtype=object)
            ex = sp.expand(sp.S(c) * _expr(top_harmonics_symmetry[top][tk]))
            blocks.setdefault((X, l, irrep, n), []).append((comp, ex, sp.S(c) * u))
    res = {}
    for key, lst in blocks.items():
        lst.sort(key=lambda t: t[0])
        res[key] = (np.array([e for _, e, _ in lst], dtype=object), np.array([u for _, _, u in lst], dtype=object).T)
    return res


# ==================================================
def internal_type(X, s, k):
    """Internal type q (polar) or g (axial) of multipolar harmonics (X, s, k)."""
    even = (s + k) % 2 == 0
    return "q" if (X == "Q") == even else "g"


# ==================================================
def create_multipolar_harmonics(pg, rank=max_rank, s_rank=max_s_rank):
    """
    Multipolar harmonics of a point group.

    Args:
        pg (str): point group tag.
        rank (int, optional): max rank l.
        s_rank (int, optional): max internal rank s.

    Returns:
        - (dict) -- {PGMultipoleType: (internal_coefficient, multipole_expression)}, sorted.
    """
    crystal = _crystal[data_harmonics_mapping[pg]["top"]]
    basis = internal_basis(crystal)
    symbols = internal_symbols(crystal)
    harm = pg_harmonics(pg, rank)

    raw = {}
    for s in range(1, s_rank + 1):
        for (X, l, irrep, n), (ex, U) in harm.items():
            for k in range(-s, s + 1):
                if l + k < 0:
                    continue
                v = create_harmonics_set(l, s=s, k=k, basis=basis)  # (2l+1, 2s+1)
                coeff = np.vectorize(lambda e: sp.factor(sp.expand(e)), otypes=[object])(U.T @ v)
                if (coeff == 0).all():
                    continue
                x = internal_type(X, s, k)
                expr = np.array([sum(c * sym for c, sym in zip(row, symbols[(s, x)])) for row in coeff], dtype=object)
                raw.setdefault((x, s, l, X, irrep, n, k), []).append((coeff, expr))

    # keep the multiplicity n of the harmonics (its xyz polynomial); number duplicates by p.
    res = {}
    for (x, s, l, X, irrep, n, k), lst in raw.items():
        for no, val in enumerate(lst):
            p = -1 if len(lst) == 1 else no + 1
            res[PGMultipoleType(X, l, irrep, n, p, s, k, x)] = val

    irreps = list(character[pg][1].keys())
    order = lambda key: (
        ["q", "g"].index(key.x),
        key.s,
        irreps.index(key.Gamma),
        key.l,
        ["Q", "G"].index(key.X),
        key.k,
        key.n,
        key.p,
    )
    return {key: res[key] for key in sorted(res, key=order)}


# ==================================================
def load_multipolar_harmonics():
    """
    Load multipolar harmonics (multipiegen/binary/multipolar_harmonics.pkl).

    Returns:
        - (BinaryManager) -- {PG_id: {PGMultipoleType: (internal_coefficient, multipole_expression)}}. str() gives the
          description.
    """
    return BinaryManager(multipolar_harmonics_file, topdir=__bin_dir__)
