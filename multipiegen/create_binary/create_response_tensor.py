# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create response_tensor.pkl (response tensors of all magnetic point groups in terms of active multipoles).

    python -m multipiegen.create_binary.create_response_tensor [--serial]

The data are taken from info.pkl (cartesian-multipole relation, LaTeX names) and group.pkl (active multipoles).
See multipiegen/util/util_response_tensor.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info
from multipiegen.util.util_response_tensor import response_tensor_all, response_tensor_file, tensor_type

h_response_tensor = """
* Response tensors of all MPG (up to rank 4).
- MPG_id (str): (dict) response tensors, "MPG:{PG}.{no}.{ID}".
  - (X, rank, opt) (str,int,str): ((ndarray(sympy), dict)) tensor matrix and expression.
    - X: tensor type, "Q" (polar), "G" (axial), "T" (magnetic polar), "M" (magnetic axial).
    - rank: tensor rank, 0-4.
    - opt: "" (rank 0,1), "s/a" (rank 2), "ss/sa/a" (rank 3), "sss/ssa/aas/aaa/sa/as/s/a/t" (rank 4).
    - matrix: tensor components, e.g. Q_{xyz} (sympy symbol), 0 for vanishing components.
      - rank 0: (1,1), rank 1: (1,3) [x,y,z].
      - rank 2 (s/a): (3,3) [x,y,z] x [x,y,z].
      - rank 3: (6,3) [xx,yy,zz,yz,zx,xy] x [x,y,z] (ss/sa), (3,3) [yz,zx,xy] x [x,y,z] (a).
      - rank 4: (6,6) (sss/ssa), (3,3) (aas/aaa), (6,3) (sa), (3,6) (as), (6,9) (s), (3,9) (a), (10,3) (t),
                with pairs v=[xx,yy,zz,yz,zx,xy], a=[yz,zx,xy], u=[xx,yy,zz,yz,zx,xy,zy,xz,yx],
                and triples t=[xxx,yyy,zzz,yyz,zzx,xxy,yzz,zxx,xyy,xyz].
    - expression: {component symbol: expression in terms of active multipoles (LaTeX symbols)}.
NOTE:
  - multipole symbol: \\mathbb{X}_{sub}^{sup(n)} with the LaTeX names of info["harmonics"]["latex"]["orbital"],
    (n) labels the independent cartesian decomposition.
  - same as MultiPie Group.response_tensor_all except for the multipole names.
"""


# ==================================================
_shared = {}


def _data():
    if not _shared:
        _shared["info"] = load_info()
        _shared["group"] = load_group()
    return _shared["info"], _shared["group"]


# ==================================================
def create_response_tensor_group(mpg_id):
    """
    Response tensors of a magnetic point group.

    Args:
        mpg_id (str): MPG id.

    Returns:
        - (dict) -- {(X, rank, opt): (matrix, expression)} for X = Q, G, T, M.
    """
    info, group = _data()
    g = group[mpg_id]
    crystal = "hexagonal" if g["info"].hexagonal_g else "cubic"
    cartesian_mp = info["response_tensor"]["cartesian_multipole"][crystal]
    latex = info["harmonics"]["latex"]["orbital"]
    d = {}
    for X in tensor_type:
        d.update(response_tensor_all(g["active_multipole"], cartesian_mp, latex, X))
    return d


def _worker(mpg_id):
    return mpg_id, create_response_tensor_group(mpg_id)


# ==================================================
def create_response_tensor(ids=None, serial=False, verbose=True):
    """
    Create response tensors (BinaryManager).

    Args:
        ids (list, optional): MPG ids (all if None).
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- response tensors.
    """
    info, _ = _data()
    if ids is None:
        ids = info["id_set"]["MPG"]["all"]

    rt = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    if serial:
        results = map(_worker, ids)
        for k, v in results:
            rt[k] = v
    else:
        with Pool(os.cpu_count()) as pool:
            for k, v in pool.imap(_worker, ids):
                rt[k] = v

    rt.add_comment(h_response_tensor)
    return rt


# ==================================================
def save_response_tensor(serial=False):
    """
    Create and save response_tensor.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_response_tensor(serial=serial, verbose=False).save_binary(response_tensor_file)
    print(f"response_tensor.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create response_tensor.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_response_tensor(parser.parse_args().serial)
