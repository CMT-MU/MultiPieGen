"""
Create group mapping to root cluster.

note:
- copy data_parent_pg.py in this directory.
- output is moved to data_root_cluster.py
"""

from data_parent_pg import data_representative_parent_pg as rep
from data_parent_pg import data_parent_pg as data

# ==================================================
mapping_oh = [
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
mapping_d6h = [
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
    [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48],
    [49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72],
    [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96],
]


# ==================================================
def create_mapping_root_cluster(c_group):
    lst = []
    while 1:
        p_group = rep[c_group]
        lst.append(c_group)
        if p_group in ["Oh", "D6h"]:
            lst.append(p_group)
            break
        c_group = p_group
    lst = lst[::-1]
    # print(lst)

    mapping = mapping_oh if lst[0] == "Oh" else mapping_d6h
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

    return mapping


# ==================================================
def create_mapping():
    dic = {}
    for c_group in rep.keys():
        mapping = create_mapping_root_cluster(c_group)
        dic[c_group] = mapping
    dic["Oh"] = mapping_oh
    dic["D6h"] = mapping_d6h

    return dic


# ==================================================
print(create_mapping())
