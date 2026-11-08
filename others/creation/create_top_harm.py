# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create top harmonics information.
    - top_harmonics_symmetry: symmetry (expression) of harmonics.
    - top_harmonics_u_matrix: u-matrix to spherical harmonics.

note:
    - u-matrix <m|gamma> is in descending order, i.e., m=l,l-1,...,-l.
    - axial irrep is obtained by interchanging "g" and "u".
    - representation matrices of all symmetry operations are common for the same irrep. irrespective of polar or axial (checked).
"""

import os
import sympy as sp
import re
from multipiegen.data.data_top_harmonics import top_harmonics
from multipiegen.util.util_multipolar_harmonics import Olm, max_rank, _create_unitary_matrix_tesseral

# ==================================================
header = """
Top harmonics infomation (polar) [axial can be obtained by interchanging g and u in irrep].

- top_harmonics_symmetry.
    { "tesseral/Oh/D6h": { tess_key/top_key: symmetry } }.
- top_harmonics_u_matrix, <m|gamma>.
    { "tesseral/Oh/D6h": { tess_key/top_key: u-matrix } }.

tess_key = (l, m), negative m is for S, otherwise for C.
top_key = (rank, irrep, multiplicity [1-based], component [0-based]).
-1 in multiplicity and component represent unique value.
"""


# ==================================================
def _remove_star(s):
    s = str(s).replace("**", "\0")
    s = re.sub(r"(\d)\*([a-zA-Z])", r"\1\2", s)
    s = s.replace("*", " ").replace("\0", "**")
    return re.sub(r"\s*([()+-])\s*", r"\1", s).strip()


# ==================================================
def create_symmetry_umat():
    sym = {"tesseral": {}, "Oh": {}, "D6h": {}}
    u_mat = {"tesseral": {}, "Oh": {}, "D6h": {}}
    for k, v in top_harmonics.items():
        top, rank, irrep, mul, comp = k
        t, e, u = v
        key = (rank, irrep, mul, comp)
        sym[top][key] = _remove_star(e)
        u_mat[top][key] = _remove_star(u)

    for l in range(max_rank + 1):
        Umg = _create_unitary_matrix_tesseral(l)
        for m in range(l, 0, -1):
            c = sp.simplify((sp.S(-1) ** m * Olm(l, m) + Olm(l, -m)) / sp.sqrt(2))
            s = sp.simplify((sp.S(-1) ** m * Olm(l, m) - Olm(l, -m)) / sp.sqrt(2) / sp.I)
            sym["tesseral"][(l, m)] = _remove_star(str(c))
            sym["tesseral"][(l, -m)] = _remove_star(str(s))
            u_mat["tesseral"][(l, m)] = _remove_star(str(Umg[:, 2 * (l - m)].tolist()))
            u_mat["tesseral"][(l, -m)] = _remove_star(str(Umg[:, 2 * (l - m) + 1].tolist()))
        sym["tesseral"][(l, 0)] = _remove_star(str(sp.simplify(Olm(l, 0))))
        u_mat["tesseral"][(l, 0)] = _remove_star(str(Umg[:, 2 * l].tolist()))

    return sym, u_mat


# ==================================================
def create_file():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "data_top_harmonics_info.py")

    sym, u_mat = create_symmetry_umat()

    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("top_harmonics_symmetry =", sym, file=f)
        print("", file=f)
        print("top_harmonics_u_matrix =", u_mat, file=f)


# ==================================================
if __name__ == "__main__":
    create_file()
