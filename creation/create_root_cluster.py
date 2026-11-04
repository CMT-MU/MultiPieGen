"""
Create root cluster.
    - oh_cluster_frac: points of D6h cluster in fractional coordinate (=d6h0 when principal axis is taken to be [111]).
    - d6h_cluster_frac: points of Oh cluster in fractional coordinate (=oh0).
    - root_cluster_cart: points of root cluster in cartesian coordinate.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import os
import numpy as np
import sympy as sp
from multipie import Group
from multipie.util.util import str_to_sympy

# ==================================================
Rc = "[[-1/sqrt(2), -1/sqrt(6), 1/sqrt(3)], [1/sqrt(2), -1/sqrt(6), 1/sqrt(3)], [0, 2/sqrt(6), 1/sqrt(3)]]"
Rx = "[[1,0,0],[0,0,-1],[0,-1,0]]"
Ry = "[[0,0,-1],[0,1,0],[-1,0,0]]"
Rz = "[[0,-1,0],[-1,0,0],[0,0,1]]"
T = "[[1, -1/2, 0], [0, sqrt(3)/2, 0], [0, 0, 1]]"

Rc = sp.Matrix(str_to_sympy(Rc))
Rci = Rc.inv()
Rx = str_to_sympy(Rx)
Ry = str_to_sympy(Ry)
Rz = str_to_sympy(Rz)
T = str_to_sympy(T)
Ti = np.array(sp.Matrix(T).inv())


# ==================================================
def create_root_cluster():
    d6h = Group("D6h")
    so = d6h.symmetry_operation["cartesian"]
    so111 = np.array([sp.simplify(Rc * i * Rci) for i in so])

    x, y, z = sp.symbols("x y z", real=True)
    d6h0 = so111 @ [x, y, z]
    d6h1 = np.array([Rx @ i for i in d6h0])
    d6h2 = np.array([Ry @ i for i in d6h0])
    d6h3 = np.array([Rz @ i for i in d6h0])

    cluster = str(sum([[i.tolist() for i in d] for d in [d6h0, d6h1, d6h2, d6h3]], [])).replace(" ", "").replace("*", "")

    return cluster


# ==================================================
def create_d6h_001_cluster():
    x, y, z = sp.symbols("x y z", real=True)
    x0 = Ti @ np.array(Rci) @ np.array([x, y, z])
    gr = Group("D6h")
    xi = gr.symmetry_operation["fractional"] @ x0
    xi = str(xi.tolist()).replace(" ", "").replace("*", "")
    return xi


# ==================================================
def create_oh_cluster():
    x, y, z = sp.symbols("x y z", real=True)
    x0 = np.array([x, y, z])
    gr = Group("Oh")
    xi = gr.symmetry_operation["fractional"] @ x0
    xi = str(xi.tolist()).replace(" ", "").replace("*", "")
    return xi


# ==================================================
def check_oh(cluster):
    M = "[[2/3, -1/3, 2/3], [-1/3, 2/3, 2/3], [2/3, 2/3, -1/3]]"
    M = str_to_sympy(M)

    oh = Group("Oh")
    so = oh.symmetry_operation["cartesian"]

    x, y, z = sp.symbols("x y z", real=True)
    oh0 = so @ [x, y, z]
    oh1 = so @ (M @ [x, y, z])

    oh_cluster = str(sum([[i.tolist() for i in d] for d in [oh0, oh1]], [])).replace(" ", "")

    cluster = sorted([str(i.tolist()) for i in str_to_sympy(cluster)])
    oh_cluster = sorted([str(i.tolist()) for i in str_to_sympy(oh_cluster)])

    return cluster == oh_cluster


# ==================================================
def create_file():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "data_root_cluster.py")

    root_cluster = create_root_cluster()
    d6h_cluster = create_d6h_001_cluster()
    oh_cluster = create_oh_cluster()

    with open(out_file, "w", encoding="utf-8") as f:
        print("oh_cluster_frac =", oh_cluster, file=f)
        print("", file=f)
        print("d6h_cluster_frac =", d6h_cluster, file=f)
        print("", file=f)
        print("root_cluster_cart =", root_cluster, file=f)

    if check_oh(root_cluster):
        print("no error")
    else:
        print("d6h-based cluster n.e. oh-based one.")


# ==================================================
create_file()
