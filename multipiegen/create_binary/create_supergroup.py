# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create supergroup.pkl (direct supergroups of each real PG and SG, with the embeddings). Needs group.pkl and info.pkl.

    python -m multipiegen.create_binary.create_supergroup [--serial]

See multipiegen/util/util_supergroup.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util.util_info import load_info
from multipiegen.util.util_supergroup import direct_supergroups, exact_matrix, supergroup_file

h_supergroup = """
* Direct supergroups of all real PG and SG (the supergroups P in which H is a maximal subgroup).
- PG_id/SG_id (str): ([dict]) direct supergroups of a group H, in order of the index and the group number ([] for none).
  - "parent" (str): (str) supergroup P.
  - "A" (str): (Matrix(sympy)) embedding x_P = A x_H (4x4, conventional coordinates; R = A[:3,:3] is right-handed).
  - "index" (str): (int) [P:H] modulo the lattice of H.
  - "index_parent" (str): (int) [P:H] modulo the lattice of P (= order of the point group of P / that of H).
NOTE:
  - maximal: a t-subgroup ([T_P:T_H] = 1, point-group index 2, 3 or 4) or a k-subgroup (the same point group,
    [T_P:T_H] = 2, 3 or 4) with no intermediate group; isomorphic supergroups (P of the type of H) are not included.
  - the embeddings are given up to the normalizer of P (the same set of operations of P in the coordinates of H) and up
    to the affine normalizer of H (A and A n give the same supergroup, n H n^-1 = H): the other embeddings are A n.
  - not always the "direct parent" of coset.pkl (the representative parent of MultiPieGen, not always maximal).
"""


# ==================================================
def _worker(gid):
    return gid, [
        {"parent": d["parent"], "A": exact_matrix(d["A"]), "index": d["index"], "index_parent": d["index_parent"]}
        for d in direct_supergroups(gid)
    ]


# ==================================================
def create_supergroup_all(serial=False, verbose=True):
    """
    Create the direct supergroups of all real PG and SG (BinaryManager).

    Args:
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- {PG_id/SG_id: [{"parent", "A", "index", "index_parent"}]}.
    """
    info = load_info()
    ids = info["id_set"]["PG"]["all"] + info["id_set"]["SG"]["all"]
    if serial:
        results = dict(map(_worker, ids))
    else:
        with Pool(os.cpu_count()) as pool:
            results = dict(pool.imap_unordered(_worker, ids))
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    for gid in ids:
        bm[gid] = results[gid]
    bm.add_comment(h_supergroup)
    return bm


# ==================================================
def save_supergroup(serial=False):
    """
    Create and save supergroup.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_supergroup_all(serial=serial, verbose=False).save_binary(supergroup_file)
    print(f"supergroup.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create supergroup.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_supergroup(parser.parse_args().serial)
