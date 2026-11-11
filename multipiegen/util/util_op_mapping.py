# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Mapping of symmetry operations to the top parent group, derived from data_parent_pg and data_parent_sg.

- op_mapping_d6h_oh_cluster (dict): mapping of point groups to the top Oh/D6h cluster (operation no. 1-based).
    - PG_tag : (top group, mapping), mapping[coset][parent op no. of child op. from 1 to g].
    - the top group is the child itself when the child is Oh or D6h.
    - the mapping of a space group to the Oh/D6h cluster is that of its point group (same point parts of operations).
- op_mapping_parent_sg (dict): mapping of space groups to the top parent space group (operation no. 1-based).
    - SG_no : (top SG_no, mapping), mapping[coset][parent op no. of child op. from 1 to g].
    - the top SG is 221-230 (Oh) or 191-194 (D6h), through the chain of maximal subgroups (data_parent_sg).
    - the axes may be changed along the chain, so the point part of a child operation may differ from that of the
      parent operation (e.g. 2[010] -> 2[100]); not a mapping to the Oh/D6h cluster.

The chain of representative parents (data_representative_parent_pg/sg) is followed from the top group down to the
child, composing the coset mappings of each step (data_parent_pg/sg): mapping[coset] lists, for each operation of the
child, the operation number of the top group.
"""

from multipiegen.data.data_parent_pg import data_parent_pg, data_representative_parent_pg
from multipiegen.data.data_parent_sg import data_parent_sg, data_representative_parent_sg

_pg_top = {"oh": ["Oh"], "d6h": ["D6h"]}
_sg_top = {"oh": [221, 222, 223, 224, 225, 226, 227, 228, 229, 230], "d6h": [191, 192, 193, 194]}
_mapping_top = {"oh": [list(range(1, 49))], "d6h": [list(range(1, 25))]}


# ==================================================
def _chain_mapping(c_group, rep, data, oh, d6h):
    lst = []
    while True:
        p_group = rep[c_group]
        lst.append(c_group)
        if p_group in oh + d6h:
            lst.append(p_group)
            break
        c_group = p_group
    lst = lst[::-1]

    mapping = _mapping_top["oh"] if lst[0] in oh else _mapping_top["d6h"]
    for i in range(1, len(lst)):
        mp = data[lst[i]][lst[i - 1]][0]
        mapping = [[m[no - 1] for no in mpr] for m in mapping for mpr in mp]
    return lst[0], mapping


# ==================================================
def create_op_mapping(rep, data, top):
    """
    Operation mapping of all groups to their top group.

    Args:
        rep (dict): representative parent of each group.
        data (dict): {child: {parent: ([[parent SO]], R)}}.
        top (dict): {"oh": [top groups], "d6h": [top groups]}.

    Returns:
        - (dict) -- {group: (top group, mapping)}.
    """
    dic = {c: _chain_mapping(c, rep, data, top["oh"], top["d6h"]) for c in rep.keys()}
    for kind in ["oh", "d6h"]:
        for i in top[kind]:
            dic[i] = (i, _mapping_top[kind])
    return dic


# ==================================================
op_mapping_d6h_oh_cluster = create_op_mapping(data_representative_parent_pg, data_parent_pg, _pg_top)
op_mapping_parent_sg = create_op_mapping(data_representative_parent_sg, data_parent_sg, _sg_top)
