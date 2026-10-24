"""
Create root cluster.
"""

import numpy as np
import sympy as sp
from multipie import Group
from multipie.util.util import str_to_sympy


# ==================================================
def create_cluster():
    Rc = "[[-1/sqrt(2), -1/sqrt(6), 1/sqrt(3)], [1/sqrt(2), -1/sqrt(6), 1/sqrt(3)], [0, 2/sqrt(6), 1/sqrt(3)]]"
    Rx = "[[1,0,0],[0,0,-1],[0,-1,0]]"
    Ry = "[[0,0,-1],[0,1,0],[-1,0,0]]"
    Rz = "[[0,-1,0],[-1,0,0],[0,0,1]]"
    Rc = sp.Matrix(str_to_sympy(Rc))
    Rci = Rc.inv()
    Rx = str_to_sympy(Rx)
    Ry = str_to_sympy(Ry)
    Rz = str_to_sympy(Rz)

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
cluster = create_cluster()
print(cluster)
if check_oh(cluster):
    print("no error")
else:
    print("d6h-based cluster n.e. oh-based one.")
