# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create cluster_samb.pkl (cluster SAMB of each Wyckoff site and bond of the real point groups and the space groups).
Needs group.pkl, cluster_scalar_basis.pkl, and cluster_bond_basis.pkl.

    python -m multipiegen.create_binary.create_cluster_samb [--serial]

See multipiegen/util/util_cluster_samb.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util.util_cluster_samb import cluster_samb_file, create_cluster_samb
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

h_cluster_samb = """
* Cluster SAMB for all real PG and SG.
- PG_id/SG_id (str): (dict) SAMB data.
  - "site" (str): (dict) site-cluster SAMB data.
    - sw_tag (str): (dict) SAMB data at site Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy)) SAMB for each component.
  - "bond_s" (str): (dict) symmetric-bond-cluster SAMB data.
    - bw_tag (str): (dict) SAMB data at bond Wyckoff, as "site".
  - "bond_a" (str): (dict) anti-symmetric-bond-cluster SAMB data.
    - bw_tag (str): (dict) SAMB data at bond Wyckoff, as "site".
NOTE:
  - the expression of the harmonics (X,l,Gamma,n) of a key is in harmonics.pkl.
  - created from the general-point SAMB of each group, i.e., cluster_scalar_basis ("site", "bond_s") and
    cluster_bond_basis (x = "q", "bond_a") of Oh or D6h restricted by the operation mapping (coset 0), and averaged
    by the mapping of each Wyckoff site and bond (with the sign of the bond reversal for "bond_a").
  - keys of the top group with the irrep and component of the group; n is the multiplicity of the harmonics (X,l,Gamma).
  - "-1" in multiplicity represents no multiplicity.
  - each component is normalized on the ns sites (bonds).
  - multipoles are sorted in order of [Gamma, l, n, p] ("site", "bond_s") or [Gamma, l, k, Q/G, n, p] ("bond_a").
"""


# ==================================================
def _worker(arg):
    gid, pg, pg_id = arg
    if "group" not in _group:
        _group["group"] = load_group()
    g = _group["group"]
    return gid, create_cluster_samb(
        pg, g[pg_id]["symmetry_operation"]["tag"], g[gid]["symmetry_operation"]["tag"], g[gid]["wyckoff"]
    )


_group = {}


# ==================================================
def group_args():
    """[(group id, point group tag, point group id)] of all real PG and SG."""
    info = load_info()
    args = [(gid, info["tag"][gid], gid) for gid in info["id_set"]["PG"]["all"]]
    group = load_group()
    for gid in info["id_set"]["SG"]["all"]:
        pg_id = group[gid]["info"].PG
        args.append((gid, info["tag"][pg_id], pg_id))
    return args


# ==================================================
def create_cluster_samb_all(serial=False, verbose=True):
    """
    Create cluster SAMB of all real PG and SG (BinaryManager).

    Args:
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- {PG_id/SG_id: {"site"/"bond_s"/"bond_a": {wyckoff: {PGMultipoleType: samb}}}}.
    """
    args = group_args()
    if serial:
        results = dict(map(_worker, args))
    else:
        with Pool(os.cpu_count()) as pool:
            results = dict(pool.imap_unordered(_worker, args))
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    for gid, _, _ in args:
        bm[gid] = results[gid]
    bm.add_comment(h_cluster_samb)
    return bm


# ==================================================
def save_cluster_samb(serial=False):
    """
    Create and save cluster_samb.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_cluster_samb_all(serial=serial, verbose=False).save_binary(cluster_samb_file)
    print(f"cluster_samb.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create cluster_samb.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_cluster_samb(parser.parse_args().serial)
