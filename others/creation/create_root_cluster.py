# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create root cluster.
    - d6h_cluster_cart: points of D6h cluster in cartesian coordinate (=d6h0 when principal axis is taken to be [111]). [same as general Wyckoff points]
    - oh_cluster_cart: points of Oh cluster in cartesian coordinate (=oh0). [same as general Wyckoff points]
    - root_cluster_cart: points of root cluster in cartesian coordinate.
"""

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
    # x0 = np.array(Rci) @ np.array([x, y, z])
    x0 = np.array([x, y, z])
    gr = Group("D6h")
    xi = gr.symmetry_operation["cartesian"] @ x0
    xi = np.array([sp.simplify(i) for i in xi])
    xi = str(xi.tolist()).replace(" ", "").replace("*", "")
    return xi


# ==================================================
def create_oh_cluster():
    x, y, z = sp.symbols("x y z", real=True)
    x0 = np.array([x, y, z])
    gr = Group("Oh")
    xi = gr.symmetry_operation["cartesian"] @ x0
    xi = np.array([sp.simplify(i) for i in xi])
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

    site0 = "[1,5,6]"
    vector0 = "[5,1,6]"
    site0h = str((np.array(Rci) @ str_to_sympy(site0)).tolist()).replace(" ", "")
    vector0h = str((np.array(Rci) @ str_to_sympy(vector0)).tolist()).replace(" ", "")

    with open(out_file, "w", encoding="utf-8") as f:

        print(f'root_cluster_cart_site0 = "{site0}"', file=f)
        print(f'oh_cluster_cart_site0 = "{site0}"', file=f)
        print(f'd6h_cluster_cart_site0 = "{site0h}"', file=f)
        print("", file=f)

        print(f'root_cluster_cart_vector0 = "{vector0}"', file=f)
        print(f'oh_cluster_cart_vector0 = "{vector0}"', file=f)
        print(f'd6h_cluster_cart_vector0 = "{vector0h}"', file=f)
        print("", file=f)

        print(f"cluster_scale = 5", file=f)
        print("", file=f)

        print(f'oh_cluster_cart = "{oh_cluster}"', file=f)
        print("", file=f)
        print(f'd6h_cluster_cart = "{d6h_cluster}"', file=f)
        print("", file=f)
        print(f'root_cluster_cart = "{root_cluster}"', file=f)

    if check_oh(root_cluster):
        print("no error")
    else:
        print("d6h-based cluster n.e. oh-based one.")


# ==================================================
if __name__ == "__main__":
    create_file()
