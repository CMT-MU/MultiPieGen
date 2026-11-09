# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Response tensors of magnetic point groups in terms of active multipoles (up to rank 4).

Ported from MultiPie (Group.response_tensor_all and Group.response_tensor), without Group: the active multipoles
are taken from group.pkl and the cartesian-multipole relation and the LaTeX names from info.pkl. The cartesian
decomposition of tensor components (P0-P4), mp_string and simplify_tensor are those of multipie.util.util_response_tensor.

- response tensor of a group (dict) = {(X, rank, opt): (matrix, expression)}
    - X: tensor type, "Q" (polar), "G" (axial), "T" (magnetic polar), "M" (magnetic axial).
    - rank: 0-4, opt: "" (rank 0, 1), "s/a" (rank 2), "ss/sa/a" (rank 3), "sss/ssa/aas/aaa/sa/as/s/a/t" (rank 4).
    - matrix: ndarray (sympy) of tensor components (symbols such as Q_{xyz}, 0 for vanishing components).
    - expression: {component symbol: expression in terms of active multipoles}.
"""

import re

import numpy as np
import sympy as sp
from multipie.util.util_binary import BinaryManager
from multipie.util.util_response_tensor import P0, P1, P2, P3, P4, mp_string, simplify_tensor, split_symbol

from multipiegen import __bin_dir__

response_tensor_file = "response_tensor"
tensor_type = ["Q", "G", "T", "M"]

_c = [1, 2, 3]
_v = [(1, 1), (2, 2), (3, 3), (2, 3), (3, 1), (1, 2)]
_u = [(1, 1), (2, 2), (3, 3), (2, 3), (3, 1), (1, 2), (3, 2), (1, 3), (2, 1)]
_a = [(2, 3), (3, 1), (1, 2)]
_t = [(1, 1, 1), (2, 2, 2), (3, 3, 3), (2, 2, 3), (3, 3, 1), (1, 1, 2), (2, 3, 3), (3, 1, 1), (1, 2, 2), (1, 2, 3)]

# (rank, opt): (row indices, column indices), as MultiPie Group.response_tensor_all.
tensor_shape = {
    (0, ""): ([()], [()]),
    (1, ""): ([()], [(i,) for i in _c]),
    (2, "s"): ([(i,) for i in _c], [(j,) for j in _c]),
    (2, "a"): ([(i,) for i in _c], [(j,) for j in _c]),
    (3, "ss"): (_v, [(j,) for j in _c]),
    (3, "sa"): (_v, [(j,) for j in _c]),
    (3, "a"): (_a, [(j,) for j in _c]),
    (4, "sss"): (_v, _v),
    (4, "ssa"): (_v, _v),
    (4, "aas"): (_a, _a),
    (4, "aaa"): (_a, _a),
    (4, "sa"): (_v, _a),
    (4, "as"): (_a, _v),
    (4, "s"): (_v, _u),
    (4, "a"): (_a, _u),
    (4, "t"): (_t, [(j,) for j in _c]),
}

_P = {0: P0, 1: P1, 2: P2, 3: P3, 4: P4}
_to_axial = {"Q": "G", "G": "Q", "T": "M", "M": "T"}
_to_magnetic = {"Q": "T", "G": "M", "T": "Q", "M": "G"}


# ==================================================
def active_dict(active, cartesian_mp):
    """
    Active multipoles for each cartesian multipole.

    Args:
        active (list): active multipoles of MPG, e.g. ["Qs", "Gdu", ...].
        cartesian_mp (dict): relation between cartesian component and multipole, {comp: [(coeff, name)]}.

    Returns:
        - (dict) -- {X+comp: [(coeff, (X, name))]}.
    """
    d = {}
    for cc, lst in cartesian_mp.items():
        for X in tensor_type:
            ex = [(c, (X, v)) for c, v in lst if X + v in active]
            if ex:
                d[X + cc] = ex
    return d


# ==================================================
def multipole_latex(X, name, ss, latex):
    """
    LaTeX symbol of a multipole, e.g. \\mathbb{Q}_{x}^{\\alpha(1)}.

    Args:
        X (str): Q/G/T/M.
        name (str): orbital name (info["harmonics"]["latex"]["orbital"]).
        ss (str): superscript of the cartesian decomposition, e.g. "(1)".
        latex (dict): LaTeX suffix of orbital names, e.g. "_{x}^{\\alpha}".

    Returns:
        - (str) -- LaTeX.
    """
    m = re.fullmatch(r"_\{(.*?)\}(?:\^\{(.*)\})?", latex[name])
    sub, sup = m.group(1), m.group(2) or ""
    return r"\mathbb{" + X + "}_{" + sub + "}^{" + sup + ss + "}"


# ==================================================
def response_tensor_component(X, rank, comp, opt, act, latex):
    """
    Tensor component in terms of active multipoles.

    Args:
        X (str): tensor type, Q/G/T/M.
        rank (int): rank, 0-4.
        comp (tuple): component indices (1:x, 2:y, 3:z).
        opt (str): tensor type option.
        act (dict): active dict (active_dict).
        latex (dict): LaTeX suffix of orbital names.

    Returns:
        - (sympy) -- expression.
    """
    if rank == 0:
        rt = P0()
    elif rank == 1:
        rt = P1(*comp)
    else:
        rt = _P[rank](*comp, opt)

    s = sp.S(0)
    for c, v in rt:
        Xc, ss, sub = split_symbol(v)
        if X in ["G", "M"]:
            Xc = _to_axial[Xc]
        if X in ["T", "M"]:
            Xc = _to_magnetic[Xc]
        for mc, (Xm, name) in act.get(Xc + sub, []):
            s += c * mc * sp.Symbol(multipole_latex(Xm, name, ss, latex))
    return s


# ==================================================
def response_tensor_all(active, cartesian_mp, latex, X):
    """
    All response tensors of type X up to rank 4.

    Args:
        active (list): active multipoles of MPG.
        cartesian_mp (dict): relation between cartesian component and multipole (cubic or hexagonal).
        latex (dict): LaTeX suffix of orbital names.
        X (str): tensor type, Q/G/T/M.

    Returns:
        - (dict) -- {(X, rank, opt): (matrix, expression)}.
    """
    act = active_dict(active, cartesian_mp)
    d = {}
    for (rank, opt), (rows, cols) in tensor_shape.items():
        m = np.zeros((len(rows), len(cols)), dtype=object)
        for ii, i in enumerate(rows):
            for jj, j in enumerate(cols):
                comp = (*i, *j)
                m[ii, jj] = (
                    mp_string(list(comp), X, False) if comp else X,
                    response_tensor_component(X, rank, comp, opt, act, latex),
                )
        d[(X, rank, opt)] = simplify_tensor(m)
    return d


# ==================================================
def load_response_tensor():
    """
    Load response tensors (multipiegen/binary/response_tensor.pkl).

    Returns:
        - (BinaryManager) -- {MPG_id: {(X, rank, opt): (matrix, expression)}}. str() gives the description.
    """
    return BinaryManager(response_tensor_file, topdir=__bin_dir__)
