# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Mapping of harmonics keys of each point group to the top harmonics, derived from data_compatibility_relation and
data_top_harmonics.top_harmonics_tag.

- data_harmonics_mapping: { PG: { "top": "Oh/D6h", "Q/G": { pg_harm_key: (coeff, top_harm_key) } } }.
    pg_harm_key = (rank, irrep, multiplicity (1-based), component (0-based)) [-1: no multiplicity/single component].
    top_harm_key = (rank, irrep, multiplicity, component).

For each top harmonics (in the order of top_harmonics_tag) and its component, the compatibility relation gives the
polar (Q) and axial (G) irrep and component of the subgroup with a coefficient. For a complex irrep of the subgroup,
the complex top harmonics (irrep + "a"/"b", single component) is referred to instead. The multiplicity of the subgroup
is n = -1 if only one top harmonics belongs to (rank, irrep, comp), otherwise n = 1, 2, ... in the order of
top_harmonics_tag.
"""

from multipiegen.data.data_compatibility_relation import compatibility_relation
from multipiegen.data.data_top_harmonics import top_harmonics_tag

HARMONICS_TYPES = ("Q", "G")


# ==================================================
def _top_irrep_comp(irrep, comp, sub_irrep):
    """(irrep, comp) of the top harmonics to be referred (complex top harmonics for a complex subgroup irrep)."""
    if "a" in sub_irrep or "b" in sub_irrep:
        if comp == 0:
            return irrep + "a", -1
        if comp == 1:
            return irrep + "b", -1
    return irrep, comp


# ==================================================
def create_mapping(group_tag):
    """
    Mapping of the harmonics of a point group to the top harmonics.

    Args:
        group_tag (str): point group tag (with "-c" for complex expression).

    Returns:
        - (dict) -- {"top": "Oh/D6h", "Q/G": {(l, irrep, n, comp): (coeff, top_key)}}.
    """
    top, cr = compatibility_relation[group_tag]
    mapping = {X: {} for X in HARMONICS_TYPES}
    for rank, irrep, mul, comp in top_harmonics_tag[top]:
        polar, axial = cr[(irrep, comp)]
        for X, (sub_irrep, sub_comp, coeff) in zip(HARMONICS_TYPES, (polar, axial)):
            top_irrep, top_comp = _top_irrep_comp(irrep, comp, sub_irrep)
            mapping[X].setdefault((rank, sub_irrep, sub_comp), []).append((coeff, (rank, top_irrep, mul, top_comp)))

    res = {X: {} for X in HARMONICS_TYPES}
    res["top"] = top
    for X, mapping_X in mapping.items():
        for (rank, irrep, comp), data in mapping_X.items():
            if len(data) == 1:
                res[X][(rank, irrep, -1, comp)] = data[0]
            else:
                for n, v in enumerate(data, start=1):
                    res[X][(rank, irrep, n, comp)] = v
    return res


# ==================================================
data_harmonics_mapping = {tag: create_mapping(tag) for tag in compatibility_relation}
