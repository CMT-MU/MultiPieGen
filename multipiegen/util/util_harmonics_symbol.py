# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Polar and axial symbols of harmonics for all point groups.

Each harmonics of a point group comes from a top (Oh/D6h) harmonics T (data_harmonics_mapping).
    polar : Q(key_p) = cp T,
    axial : G(key_a) = ca T,
so that the polar harmonics Q(key_p) and the axial harmonics (cp ca) G(key_a) have the same
expression cp T. The axial symbol differs from the polar one only by irrep, component, and sign.

- harmonics key = (rank, irrep, multiplicity [1-based], component [0-based]), -1: unique.
"""

import sympy as sp
from multipie.util.util import str_to_sympy
from multipiegen.data.data_character_table import character
from multipiegen.util.util_harmonics_mapping import data_harmonics_mapping
from multipiegen.util.util_top_harmonics import top_harmonics_symmetry


# ==================================================
def is_real_irrep(irrep):
    """
    Is real irrep of top group (Oh/D6h) ? (complex combinations have irrep names ending with "a" or "b").

    Args:
        irrep (str): irrep name of top group.

    Returns:
        - (bool) -- real irrep ?
    """
    return irrep[-1] in "gu"


# ==================================================
def irrep_latex(irrep):
    """
    Irrep in LaTeX (same format as MultiPie Ver.2).
        e.g., "T1u" -> "T_{1u}", "A2''" -> "A_{2}^{\\prime\\prime}", "E1ua" -> "E_{1u}^{(a)}", "Ea'" -> "E^{\\prime(a)}".

    Args:
        irrep (str): irrep name.

    Returns:
        - (str) -- irrep in LaTeX.
    """
    body = irrep.rstrip("'")
    prime = "\\prime" * (len(irrep) - len(body))
    cmplx = ""
    if len(body) > 1 and body[-1] in "ab":  # complex irrep, e.g., Ea, E1ga.
        cmplx = f"({body[-1]})"
        body = body[:-1]
    s = body[0] + (f"_{{{body[1:]}}}" if len(body) > 1 else "")
    if prime or cmplx:
        s += "^{" + prime + cmplx + "}"
    return s


# ==================================================
def multipole_symbol(X, key, sign=1):
    """
    Multipole symbol in LaTeX (same format as MultiPie Ver.2), e.g., "\\mathbb{Q}_{1,1}(T_{1u})".

    Args:
        X (str): "Q" (polar) or "G" (axial).
        key (tuple): harmonics key (rank, irrep, multiplicity, component).
        sign (int, optional): sign of symbol (+1/-1).

    Returns:
        - (str) -- multipole symbol in LaTeX.
    """
    rank, irrep, mul, comp = key
    sub = f"{rank}" if comp == -1 else f"{rank},{comp + 1}"
    arg = irrep_latex(irrep) if mul == -1 else f"{irrep_latex(irrep)},{mul}"
    s = f"\\mathbb{{{X}}}_{{{sub}}}({arg})"
    return ("-" if sign < 0 else "") + s


# ==================================================
def polar_to_axial(group, key):
    """
    Axial key and sign of the axial harmonics having the same expression as a polar harmonics.

    Args:
        group (str): point group tag.
        key (tuple): polar harmonics key (rank, irrep, multiplicity, component).

    Returns:
        - (tuple) -- axial harmonics key (rank, irrep, multiplicity, component).
        - (int) -- sign (+1/-1), polar Q(key) and sign*G(axial key) have the same expression.
    """
    mapping = data_harmonics_mapping[group]
    cp, top_key = mapping["Q"][key]
    for a_key, (ca, t_key) in mapping["G"].items():
        if t_key == top_key:
            return a_key, int(cp) * int(ca)
    raise KeyError(f"no axial harmonics of {group} comes from {top_key}.")


# ==================================================
def _sort_key(group):
    order = {irrep: i for i, irrep in enumerate(character[group][1].keys())}
    return lambda k: (k[0], order[k[1]], k[2], k[3])


# ==================================================
def harmonics_symbol_table(group):
    """
    Table of harmonics with polar and axial symbols.

    Args:
        group (str): point group tag.

    Returns:
        - (list) -- [(key, polar symbol, axial symbol, expression in LaTeX)],
          sorted by (rank, irrep in character table, multiplicity, component) as in MultiPie Ver.2.
    """
    mapping = data_harmonics_mapping[group]
    top = mapping["top"]
    axial = {t_key: (a_key, ca) for a_key, (ca, t_key) in mapping["G"].items()}

    tbl = []
    for key in sorted(mapping["Q"].keys(), key=_sort_key(group)):
        cp, top_key = mapping["Q"][key]
        a_key, ca = axial[top_key]
        expr = sp.latex(cp * str_to_sympy(top_harmonics_symmetry[top][top_key]))
        tbl.append((key, multipole_symbol("Q", key), multipole_symbol("G", a_key, int(cp) * int(ca)), expr))
    return tbl
