# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create operation mapping to D6h and Oh, and root cluster.
    - op_mapping_d6h_oh_cluster: operation mapping of point groups to D6h and Oh cluster.
    - op_mapping_parent_sg: operation mapping of space groups to the top parent space group (221-230, 191-194).
    - op_mapping_root_cluster: operation mapping to root cluster.
"""

import os
from multipiegen.data.data_parent_pg import data_representative_parent_pg as p_rep
from multipiegen.data.data_parent_pg import data_parent_pg as p_data
from multipiegen.data.data_parent_sg import data_representative_parent_sg as s_rep
from multipiegen.data.data_parent_sg import data_parent_sg as s_data

# ==================================================
mapping_rc_oh = [
    [
        1,
        92,
        67,
        45,
        2,
        43,
        93,
        68,
        3,
        69,
        44,
        91,
        85,
        8,
        62,
        39,
        86,
        37,
        9,
        63,
        87,
        61,
        38,
        7,
        13,
        80,
        55,
        33,
        14,
        31,
        81,
        56,
        15,
        57,
        32,
        79,
        73,
        20,
        50,
        27,
        74,
        25,
        21,
        51,
        75,
        49,
        26,
        19,
    ],
    [
        23,
        76,
        53,
        30,
        22,
        29,
        78,
        52,
        24,
        54,
        28,
        77,
        83,
        16,
        58,
        36,
        82,
        35,
        18,
        60,
        84,
        59,
        34,
        17,
        11,
        88,
        65,
        42,
        10,
        41,
        90,
        64,
        12,
        66,
        40,
        89,
        95,
        4,
        70,
        48,
        94,
        47,
        6,
        72,
        96,
        71,
        46,
        5,
    ],
]
mapping_rc_d6h = [
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
    [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48],
    [49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72],
    [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96],
]

mapping_oh = [list(range(1, 49))]
mapping_d6h = [list(range(1, 25))]


p_oh = ["Oh"]
p_d6h = ["D6h"]
s_oh = [221, 222, 223, 224, 225, 226, 227, 228, 229, 230]
s_d6h = [191, 192, 193, 194]


# ==================================================
def create_mapping_root_cluster(c_group, rep, data, oh, d6h, mapping_oh, mapping_d6h):
    lst = []
    while 1:
        p_group = rep[c_group]
        lst.append(c_group)
        if p_group in oh + d6h:
            lst.append(p_group)
            break
        c_group = p_group
    lst = lst[::-1]
    # print(lst)

    mapping = mapping_oh if lst[0] in oh else mapping_d6h
    # print(lst[0], "=>", lst[1], ":", mapping)
    for i in range(1, len(lst)):
        p_group = lst[i - 1]
        c_group = lst[i]
        mp = data[c_group][p_group][0]
        mapping1 = []
        for mapping2 in mapping:
            for mpr in mp:
                mapping1.append([mapping2[no - 1] for no in mpr])
        # print(p_group, "=>", c_group, ":", mapping1)
        mapping = mapping1

    return lst[0], mapping


# ==================================================
def create_mapping(rep, data, oh, d6h, mapping_oh, mapping_d6h, rc):
    dic = {}
    for c_group in rep.keys():
        top, mapping = create_mapping_root_cluster(c_group, rep, data, oh, d6h, mapping_oh, mapping_d6h)
        if rc:
            dic[c_group] = ("RC", mapping)
        else:
            dic[c_group] = (top, mapping)
    if rc:
        for i in oh:
            dic[i] = ("RC", mapping_oh)
        for i in d6h:
            dic[i] = ("RC", mapping_d6h)
    else:
        for i in oh:
            dic[i] = (i, mapping_oh)
        for i in d6h:
            dic[i] = (i, mapping_d6h)

    return dic


# ==================================================
def create_mapping_all():
    # mapping to root cluster.
    rc_p_mapping = create_mapping(p_rep, p_data, p_oh, p_d6h, mapping_rc_oh, mapping_rc_d6h, rc=True)
    rc_s_mapping = create_mapping(s_rep, s_data, s_oh, s_d6h, mapping_rc_oh, mapping_rc_d6h, rc=True)
    root_mapping = rc_p_mapping | rc_s_mapping

    # mapping to D6h and Oh (PG), and to the top parent SG (SG).
    p_mapping = create_mapping(p_rep, p_data, p_oh, p_d6h, mapping_oh, mapping_d6h, rc=False)
    s_mapping = create_mapping(s_rep, s_data, s_oh, s_d6h, mapping_oh, mapping_d6h, rc=False)

    return root_mapping, p_mapping, s_mapping


# ==================================================
def create_file():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "data_op_mapping.py")

    root_mapping, p_mapping, s_mapping = create_mapping_all()

    with open(out_file, "w", encoding="utf-8") as f:
        print("op_mapping_d6h_oh_cluster =", p_mapping, file=f)
        print("", file=f)
        print("op_mapping_parent_sg =", s_mapping, file=f)
        print("", file=f)
        print("op_mapping_root_cluster =", root_mapping, file=f)


# ==================================================
if __name__ == "__main__":
    create_file()
