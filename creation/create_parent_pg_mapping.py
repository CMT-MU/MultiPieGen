"""
Create parent point group mapping.
"""

import os
import sympy as sp
import numpy as np
from multipie import Group
from multipie.util.util import str_to_sympy

# ==================================================
header = """
Parent group information for all point groups.

- data_parent_pg.
    { child PG_id: { parent PG_id: ([[parent SO]], R) } }.
    SO: 1-based, [parent SO] in order of child-PG SO.
    R: rotaion matrix for converting principal axis (fractional).

- data_representative_parent_pg.
    { child PG_id: parent_PG_id }.
"""

# ==================================================
# rotaion matrices for converting principal axes.
Rc = "[[-1/sqrt(2), -1/sqrt(6), 1/sqrt(3)], [1/sqrt(2), -1/sqrt(6), 1/sqrt(3)], [0, 2/sqrt(6), 1/sqrt(3)]]"
Rc1 = "[[-1/sqrt(6), 1/sqrt(2), 1/sqrt(3)], [-1/sqrt(6), -1/sqrt(2), 1/sqrt(3)], [2/sqrt(6), 0, 1/sqrt(3)]]"
Rt = "[[sqrt(2)/2, sqrt(2)/2, 0], [-sqrt(2)/2, sqrt(2)/2, 0], [0, 0, 1]]"
Rx = "[[1, 0, 0], [0, 0, 1], [0, -1, 0]]"
Ry = "[[0, 0, -1], [0, 1, 0], [1, 0, 0]]"
Rz = "[[0, 1, 0], [-1, 0, 0], [0, 0, 1]]"
T = "[[1,-1/2,0],[0,sqrt(3)/2,0],[0,0,1]]"

Rc = sp.Matrix(str_to_sympy(Rc))
Rc1 = sp.Matrix(str_to_sympy(Rc1))
Rt = sp.Matrix(str_to_sympy(Rt))
Rx = sp.Matrix(str_to_sympy(Rx))
Ry = sp.Matrix(str_to_sympy(Ry))
Rz = sp.Matrix(str_to_sympy(Rz))
T = sp.Matrix(str_to_sympy(T))

# ==================================================
# child -> (parent, rotation matrix if necessary).
lst = [
    # triclinic.
    ("C1", ["Ci", "C2", "Cs", "C3"]),
    ("Ci", ["C2h", "C3i"]),
    # monoclinic.
    ("C2", ["C2h", "D2", ("C2v", Rx), ("C4", Rx), ("S4", Rx), ("D3", T.inv() * Rz), ("D3-1", T.inv()), ("C6", T.inv() * Rx)]),
    ("Cs", ["C2h", "C2v", ("C3v", T.inv() * Rz), ("C3v-1", T.inv()), ("C3h", T.inv() * Rx)]),
    ("C2h", ["D2h", ("C4h", Rx), ("C6h", T.inv() * Rx), ("D3d", T.inv() * Rz), ("D3d-1", T.inv())]),
    # orthorhombic.
    ("D2", ["D2h", "D4", "D2d", ("D2d-1", Rt), ("D6", T.inv()), "T"]),
    ("C2v", ["D2h", "C4v", ("D2d", Rt), "D2d-1", ("C6v", T.inv()), ("D3h", T.inv() * Rx), ("D3h-1", T.inv() * Ry)]),
    ("D2h", ["D4h", ("D6h", T.inv()), "Th"]),
    # tetragonal.
    ("C4", ["C4h", "D4", "C4v"]),
    ("S4", ["C4h", "D2d", "D2d-1"]),
    ("C4h", ["D4h"]),
    ("D4", ["D4h", "O"]),
    ("C4v", ["D4h"]),
    ("D2d", ["D4h", "Td"]),
    ("D2d-1", ["D4h", ("Td", Rt)]),
    ("D4h", ["Oh"]),
    # trigonal.
    ("C3", ["C3i", "D3", "D3-1", "C3v", "C3v-1", "C6", "C3h", ("T", Rc * T)]),
    ("C3i", ["D3d", "D3d-1", "C6h", ("Th", Rc * T)]),
    ("D3", ["D3d", ("D3d-1", T.inv() * Rz * T), "D6", ("D3h", T.inv() * Rz * T), "D3h-1", ("O", Rc * T)]),
    ("D3-1", ["D3d-1", ("D3d", T.inv() * Rz * T), "D6", "D3h", ("D3h-1", T.inv() * Rz * T), ("O", Rc1 * T)]),
    ("C3v", ["D3d", ("D3d-1", T.inv() * Rz * T), "C6v", "D3h", ("D3h-1", T.inv() * Rz * T), ("Td", Rc * T)]),
    ("C3v-1", ["D3d-1", ("D3d", T.inv() * Rz * T), "C6v", ("D3h", T.inv() * Rz * T), "D3h-1", ("Td", Rc1 * T)]),
    ("D3d", ["D6h", ("Oh", Rc * T)]),
    ("D3d-1", ["D6h", ("Oh", Rc1 * T)]),
    # hexagonal.
    ("C6", ["C6h", "D6", "C6v"]),
    ("C3h", ["C6h", "D3h", "D3h-1"]),
    ("C6h", ["D6h"]),
    ("D6", ["D6h"]),
    ("C6v", ["D6h"]),
    ("D3h", ["D6h"]),
    ("D3h-1", ["D6h"]),
    # ("D6h", []),
    # cubic.
    ("T", ["Th", "O", "Td"]),
    ("Th", ["Oh"]),
    ("O", ["Oh"]),
    ("Td", ["Oh"]),
    # ("Oh", []),
]


# ==================================================
def create_parent_group_part(c_group, p_group, add_lst=[], R=None):
    g_p = Group(p_group)
    so_p = g_p.symmetry_operation["fractional"]
    g = Group(c_group)
    so = g.symmetry_operation["fractional"]
    if R is not None:
        so = np.array([sp.simplify(R * sp.Matrix(o) * R.inv()).tolist() for o in so])

    lst = so_p.tolist()
    result = {}
    for no, o in enumerate(so.tolist()):
        o1 = lst.index(o)
        result[no + 1] = [o1 + 1]
        for add in add_lst:
            add_so = so_p[add]
            so_add = add_so @ so[no]
            result[no + 1].append(lst.index(so_add.tolist()) + 1)
    result = np.array([r for r in result.values()]).T.tolist()
    return result


# ==================================================
def create_parent_group(c_group, p_group, R):
    add_lst = []
    top = 0
    while top != -1:
        result = create_parent_group_part(c_group, p_group, add_lst, R)
        num = sorted(sum(result, []))
        top = next((i + 1 for i in range(len(num) - 1) if num[i + 1] != num[i] + 1), -1)
        n = len(set(sum(result, [])))
        if n < len(Group(p_group).symmetry_operation["tag"]) and top == -1:
            top = n
        add_lst.append(top)

    assert len(set(sum(result, []))) == len(Group(p_group).symmetry_operation["tag"])
    return result


# ==================================================
def create():
    dic = {}  # 1-based.
    rep = {}
    for c_group, p_groups in lst:
        top = p_groups[0]
        rep[c_group] = top[0] if type(top) is tuple else top
        for p_group in p_groups:
            pp = p_group[0] if type(p_group) == tuple else p_group
            print(c_group, "<=", pp)
            R = None
            if type(p_group) is tuple:
                p_group, R = p_group
            mapping = create_parent_group(c_group, p_group, R)
            if R is not None:
                R = str(R.tolist()).replace(" ", "")
            else:
                R = "[[1,0,0],[0,1,0],[0,0,1]]"
            if c_group in dic.keys():
                dic[c_group][p_group] = (mapping, R)
            else:
                dic[c_group] = {p_group: (mapping, R)}

    return dic, rep


# ==================================================
def create_file():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "data_parent_pg.py")
    dic, rep = create()
    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("data_parent_pg =", dic, file=f)
        print("", file=f)
        print("data_representative_parent_pg =", rep, file=f)


create_file()
