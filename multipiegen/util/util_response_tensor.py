# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Response tensors of magnetic point groups in terms of active multipoles (up to rank 4).

Ported from MultiPie Ver.2 (Group.response_tensor_all and Group.response_tensor), without Group: the active multipoles
are taken from group.pkl and the cartesian-multipole relation and the LaTeX names from info.pkl. The cartesian
decomposition of tensor components (P0-P4), mp_string and simplify_tensor are those of multipie.util.util_response_tensor,
except for the rank-4 t part (T1234), which is defined here (fully symmetric; that of MultiPie Ver.2 and 3.0-3.2.0 is
not), so that the tensors do not depend on the installed MultiPie.

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

# (rank, opt): (row indices, column indices), as MultiPie Ver.2 Group.response_tensor_all.
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


# ==================================================
def T1234(i1, i2, i3, i4):
    """
    T part of rank 4 tensor component (fully symmetric in i1, i2, i3).

    The same as multipie.util.util_response_tensor.T1234 of MultiPie 3.2.1 or later (CMT-MU/MultiPie#255), so that the
    tensors do not depend on the installed MultiPie (that of MultiPie Ver.2 and 3.0-3.2.0 is not fully symmetric).

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """

    def delta(i, j):
        return int(i == j)

    def epsilon(i, j, k):
        return int(sp.LeviCivita(i, j, k))

    # (d12 d34 + d13 d24 + d14 d23) Q^{(1)}.
    d = delta(i1, i2) * delta(i3, i4) + delta(i1, i3) * delta(i2, i4) + delta(i1, i4) * delta(i2, i3)
    l1 = [(d, "Q^{(1)}")] if d != 0 else []

    # d12 Q^{(1)}_{34} + d13 Q^{(1)}_{24} + d23 Q^{(1)}_{14}.
    l2 = [(1, mp_string(k, "Q^{(1)}")) for a, b, k in [(i1, i2, [i3, i4]), (i1, i3, [i2, i4]), (i2, i3, [i1, i4])] if a == b]

    # d14 Q^{(2)}_{23} + d24 Q^{(2)}_{13} + d34 Q^{(2)}_{12}.
    l3 = [(1, mp_string(k, "Q^{(2)}")) for a, b, k in [(i1, i4, [i2, i3]), (i2, i4, [i1, i3]), (i3, i4, [i1, i2])] if a == b]

    # (d12 e345 + d13 e245 + d23 e145) G^{(1)}_{5}.
    l4 = []
    for i5 in range(1, 4):
        e = delta(i1, i2) * epsilon(i3, i4, i5) + delta(i1, i3) * epsilon(i2, i4, i5) + delta(i2, i3) * epsilon(i1, i4, i5)
        if e != 0:
            l4.append((e, mp_string([i5], "G^{(1)}")))

    # e345 G^{(1)}_{125} + e145 G^{(1)}_{235} + e245 G^{(1)}_{135}.
    l5 = []
    for i5 in range(1, 4):
        for a, b, c in [(i1, i2, i3), (i2, i3, i1), (i1, i3, i2)]:
            e = epsilon(c, i4, i5)
            if e != 0:
                l5.append((e, mp_string([a, b, i5], "G^{(1)}")))

    # Q^{(1)}_{1234}.
    l6 = [(1, mp_string([i1, i2, i3, i4], "Q^{(1)}"))]

    return l1 + l2 + l3 + l4 + l5 + l6


# ==================================================
def _P4(i1, i2, i3, i4, opt=None):
    """
    Rank 4 tensor component, P4 of MultiPie with the t part of T1234.
    """
    if opt == "t":
        return T1234(i1, i2, i3, i4)
    return P4(i1, i2, i3, i4, opt)


_P = {0: P0, 1: P1, 2: P2, 3: P3, 4: _P4}
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
