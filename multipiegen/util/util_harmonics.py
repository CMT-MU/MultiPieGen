# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Harmonics of all point groups: cartesian expression (symmetry) and u-matrix for each multipole key.

The binary data keyed by PGMultipoleType (multipolar_harmonics.pkl, atomic_multipole_group.pkl, cluster_samb.pkl, the
cluster bases) do not hold the expression of the harmonics; look it up here by the key.

- Harmonics. Every harmonics (X, l, Gamma, n) of a point group is c times a top harmonics T of Oh or D6h
  (data_harmonics_mapping); its expression is c T (data_top_harmonics_info.top_harmonics_symmetry) and its u-matrix is
  c U_T (top_harmonics_u_matrix, rows m = l, l-1, ..., -l, columns the components).
- Key. PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q") with X = Q or G. Lookups (harmonics_expression,
  harmonics_u_matrix) accept any PGMultipoleType (or tuple) and use its (X, l, Gamma, n) only, so the keys of the
  other binaries can be used as they are; X = T and M (atomic multipoles) use the harmonics of Q and G. n is the
  multiplicity of the harmonics (1-based, -1 if unique).
- Order. Gamma (character table), l, X (Q, G), n.
- Range. 47 point groups, rank l <= 11.

- harmonics (dict)::

    { PG_id: { PGMultipoleType: (expression, u_matrix) } }
        expression: ndarray (dim) (sympy), cartesian expression (polynomial in x, y, z) of each component.
        u_matrix: ndarray (2l+1, dim) (sympy), <m|component> for m = l, l-1, ..., -l.
"""

from functools import lru_cache

import numpy as np
from multipie import PGMultipoleType
from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.data.data_character_table import character
from multipiegen.util.util_harmonics_multipole import pg_harmonics
from multipiegen.util.util_info import load_info

# ==================================================
harmonics_file = "harmonics"
max_rank = 11


# ==================================================
def create_harmonics(pg, rank=max_rank):
    """
    Harmonics of a point group.

    Args:
        pg (str): point group tag.
        rank (int, optional): max rank.

    Returns:
        - (dict) -- {PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q"): (expression ndarray(dim), u_matrix ndarray(2l+1, dim))}.
    """
    irreps = list(character[pg][1].keys())
    harm = pg_harmonics(pg, rank)
    order = sorted(harm, key=lambda k: (irreps.index(k[2]), k[1], ["Q", "G"].index(k[0]), k[3]))
    return {PGMultipoleType(X, l, irrep, n, -1, 0, 0, "q"): harm[(X, l, irrep, n)] for X, l, irrep, n in order}


# ==================================================
def load_harmonics():
    """
    Load harmonics (multipiegen/binary/harmonics.pkl).

    Returns:
        - (BinaryManager) -- {PG_id: {PGMultipoleType: (expression, u_matrix)}}.
    """
    return BinaryManager(harmonics_file, topdir=__bin_dir__)


# ==================================================
@lru_cache(maxsize=1)
def _data():
    return load_harmonics()


def _pg_id(pg):
    return pg if pg.startswith("PG:") else load_info()["id"][pg]


# harmonics type of each multipole type: T (magnetic toroidal) has the polar harmonics Q, M (magnetic) the axial G.
harmonics_type = {"Q": "Q", "G": "G", "T": "Q", "M": "G"}


def _harmonics(pg, key):
    X, l, Gamma, n = tuple(key)[:4]
    X = harmonics_type[X]
    return _data()[_pg_id(pg)][PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q")]


# ==================================================
def harmonics_expression(pg, key):
    """
    Cartesian expression of the harmonics of a multipole key.

    Args:
        pg (str): point group tag or id, e.g. "Oh" or "PG:32".
        key (PGMultipoleType or tuple): multipole key; only (X, l, Gamma, n) is used.

    Returns:
        - (ndarray) -- (dim) expression (sympy) of each component.
    """
    return _harmonics(pg, key)[0]


# ==================================================
def harmonics_u_matrix(pg, key):
    """
    u-matrix of the harmonics of a multipole key.

    Args:
        pg (str): point group tag or id, e.g. "Oh" or "PG:32".
        key (PGMultipoleType or tuple): multipole key; only (X, l, Gamma, n) is used.

    Returns:
        - (ndarray) -- (2l+1, dim) u-matrix (sympy), <m|component> for m = l, l-1, ..., -l.
    """
    return np.asarray(_harmonics(pg, key)[1])
