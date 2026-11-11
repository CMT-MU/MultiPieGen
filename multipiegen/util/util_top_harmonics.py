# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Top harmonics information (polar), derived from data_top_harmonics.

- top_harmonics_symmetry: { "tesseral"/"Oh"/"D6h": { tess_key/top_key: symmetry (str, cartesian expression) } }.
- top_harmonics_u_matrix: { "tesseral"/"Oh"/"D6h": { tess_key/top_key: u-matrix (str), <m|gamma>, m = l, ..., -l } }.

tess_key = (l, m), negative m is for S, otherwise for C. top_key = (rank, irrep, multiplicity [1-based], component
[0-based]); -1 in multiplicity and component represents a unique value. The axial irrep is obtained by interchanging
"g" and "u".

"Oh" and "D6h" are the expressions and u-matrices of data_top_harmonics in a compact form. "tesseral" is computed from
the multipolar spherical harmonics (rank <= 11) when it is first used (a few seconds).
"""

import re
from collections.abc import Mapping

import sympy as sp

from multipiegen.data.data_top_harmonics import top_harmonics
from multipiegen.util.util_multipolar_harmonics import Olm, _create_unitary_matrix_tesseral, max_rank


# ==================================================
def _remove_star(s):
    s = str(s).replace("**", "\0")
    s = re.sub(r"(\d)\*([a-zA-Z])", r"\1\2", s)
    s = s.replace("*", " ").replace("\0", "**")
    return re.sub(r"\s*([()+-])\s*", r"\1", s).strip()


# ==================================================
def _top():
    sym, u_mat = {"Oh": {}, "D6h": {}}, {"Oh": {}, "D6h": {}}
    for (top, rank, irrep, mul, comp), (_, e, u) in top_harmonics.items():
        sym[top][(rank, irrep, mul, comp)] = _remove_star(e)
        u_mat[top][(rank, irrep, mul, comp)] = _remove_star(u)
    return sym, u_mat


_tesseral_cache = {}


def _tesseral():
    if not _tesseral_cache:
        sym, u_mat = {}, {}
        for l in range(max_rank + 1):
            Umg = _create_unitary_matrix_tesseral(l)
            for m in range(l, 0, -1):
                c = sp.simplify((sp.S(-1) ** m * Olm(l, m) + Olm(l, -m)) / sp.sqrt(2))
                s = sp.simplify((sp.S(-1) ** m * Olm(l, m) - Olm(l, -m)) / sp.sqrt(2) / sp.I)
                sym[(l, m)] = _remove_star(str(c))
                sym[(l, -m)] = _remove_star(str(s))
                u_mat[(l, m)] = _remove_star(str(Umg[:, 2 * (l - m)].tolist()))
                u_mat[(l, -m)] = _remove_star(str(Umg[:, 2 * (l - m) + 1].tolist()))
            sym[(l, 0)] = _remove_star(str(sp.simplify(Olm(l, 0))))
            u_mat[(l, 0)] = _remove_star(str(Umg[:, 2 * l].tolist()))
        _tesseral_cache["symmetry"], _tesseral_cache["u_matrix"] = sym, u_mat
    return _tesseral_cache


# ==================================================
class _TopHarmonics(Mapping):
    """{"tesseral"/"Oh"/"D6h": {...}}, with "tesseral" computed on first access."""

    def __init__(self, top, kind):
        self._top, self._kind = top, kind

    def __getitem__(self, key):
        if key == "tesseral":
            return _tesseral()[self._kind]
        return self._top[key]

    def __iter__(self):
        return iter(["tesseral", "Oh", "D6h"])

    def __len__(self):
        return 3


_sym, _u = _top()
top_harmonics_symmetry = _TopHarmonics(_sym, "symmetry")
top_harmonics_u_matrix = _TopHarmonics(_u, "u_matrix")
