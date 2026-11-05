# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from multipiegen import __top_dir__
import os

os.chdir(__top_dir__)

import numpy as np
import sympy as sp
from multipie import Group
from multipie.util.util import str_to_sympy
from multipiegen.data.data_root_cluster import root_cluster_cart, oh_cluster_frac, d6h_cluster_frac, root_cluster_cart_site0
from multipiegen.data.data_op_mapping import op_mapping_root_cluster, op_mapping_d6h_oh_cluster
from multipiegen.data.data_parent_pg import data_parent_pg, data_representative_parent_pg
from multipiegen.data.data_parent_sg import data_parent_sg, data_representative_parent_sg

# ==================================================
chiral_pg = ["C1", "C2", "D2", "C4", "D4", "C3", "D3", "D3-1", "C6", "D6", "T", "O"]
sohncke_sg = [  # 65 - 22 = 43
    1,
    3,
    4,
    5,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    75,
    77,
    79,
    80,
    89,
    90,
    93,
    94,
    97,
    98,
    143,
    146,
    149,
    150,
    155,
    168,
    173,
    177,
    182,
    195,
    196,
    197,
    198,
    199,
    207,
    208,
    209,
    210,
    211,
    214,
]
chiral_sg = [  # 11x2 = 22, (right, left).
    (76, 78),
    (91, 95),
    (92, 96),
    (144, 145),
    (151, 153),
    (152, 154),
    (169, 170),
    (171, 172),
    (178, 179),
    (180, 181),
    (212, 213),
]


# ==================================================
def _shift(a, var):
    """
    shift for space group.

    Args:
        a (ndarray): wyckoff positions.
        var (tuple): (x,y,z) variables.

    Returns:
        - (ndarray) -- shifted wyckoff positions.
    """
    site = a.shape[1] == 3
    if site:
        a = np.vectorize(
            lambda e: sp.expand(e - e.subs({var[0]: 0, var[1]: 0, var[2]: 0})) + e.subs({var[0]: 0, var[1]: 0, var[2]: 0}) % 1,
            otypes=[object],
        )(a)
    else:
        pass

    return a


# ==================================================
def _inversion(a, var):
    """
    inversion + minus sign with shift.

    Args:
        a (ndarray): wyckoff positions.
        var (tuple): (x,y,z) variables.

    Returns:
        - (ndarray) -- inversion + minus sign.
    """
    site = a.shape[1] == 3
    if site:
        a = np.vectorize(lambda i: -i.subs({"x": -var[0], "y": -var[1], "z": -var[2]}))(a)
    else:
        pass

    return _shift(a, var)


# ==================================================
def create_coset_so():
    """
    create symmetry operations of coset.
    """
    dic = {}
    for c_group in data_parent_pg.keys():
        p_group = list(data_parent_pg[c_group].keys())[0]
        mapping, _ = data_parent_pg[c_group][p_group]
        pg = Group(p_group)
        so = [pg.symmetry_operation["tag"][mp[0] - 1] for mp in mapping[1:]]
        dic[c_group] = (p_group, so)
    for c_group in data_parent_sg.keys():
        p_group = list(data_parent_sg[c_group].keys())[0]
        mapping, _ = data_parent_sg[c_group][p_group]
        pg = Group(p_group)
        so = [pg.symmetry_operation["tag"][mp[0] - 1] for mp in mapping[1:]]
        dic[c_group] = (p_group, so)

    outfile = "data_coset_so.py"
    f = open(outfile, "w")
    print('"""\nSymmetry operations of coset.\n\n- data_coset_so = { child: (parent, [SO]) }\n"""\n', file=f)
    print("data_coset_so =", dic, file=f)


# ==================================================
def check_wyckoff_inversion():
    """
    check whether -wyckofff(-x,-y,-z) coincides with wyckoff (x,y,z) ? [for chiral_pg and sohncke_sg]
    """

    var = Group.global_info()["harmonics"]["variable"]

    for group in [chiral_pg, sohncke_sg]:
        for c_group in group:
            print("===", c_group, "===")
            cg = Group(c_group)
            print("--- site ---")
            for wp in cg.wyckoff["site"].keys():
                wyckoff = cg.wyckoff["site"][wp]["conventional"]
                wyckoff_bar = _inversion(wyckoff, var)
                print(" ", wp, str(wyckoff.tolist()).replace(" ", "") == str(wyckoff_bar.tolist()).replace(" ", ""))
            """
            print("--- bond ---")
            for wp in cg.wyckoff["bond"].keys():
                wyckoff = cg.wyckoff["bond"][wp]["conventional"]
                wyckoff_bar = _inversion(wyckoff, var)
                print(" ", wp, str(wyckoff.tolist()).replace(" ", "") == str(wyckoff_bar.tolist()).replace(" ", ""))
            """


# ==================================================
# create_coset_so()
check_wyckoff_inversion()
