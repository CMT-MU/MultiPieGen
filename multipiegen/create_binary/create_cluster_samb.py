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

import numpy as np

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util.util_cluster_samb import cluster_samb_file, create_cluster_samb
from multipiegen.util.util_cluster_samb_image import cluster_samb_reference, inversion_pairs, orbit_correspondence
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

h_cluster_samb = """
* Cluster SAMB for all real PG and SG.
- PG_id/SG_id (str): (dict) SAMB data.
  - "site" (str): (dict) site-cluster SAMB data.
    - sw_tag (str): (dict) SAMB data at site Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy)) SAMB for each component.
  - "bond" (str): (dict) bond-cluster SAMB data.
    - bw_tag (str): (dict) SAMB data at bond Wyckoff, as "site".
NOTE:
  - the expression of the harmonics (X,l,Gamma,n) of a key is in harmonics.pkl (T and M use those of Q and G).
  - created from the general-point SAMB of each group, i.e., cluster_scalar_basis ("site", symmetric bond) and
    cluster_bond_basis (x = "q", anti-symmetric bond) of Oh or D6h restricted by the operation mapping (coset 0), and
    averaged by the mapping of each Wyckoff site and bond (with the sign of the bond reversal for anti-symmetric bond).
  - "bond" = symmetric bond SAMB (X = Q, real) + i anti-symmetric bond SAMB (X = Q -> T, G -> M, pure imaginary).
    Q and T/M are orthonormalized separately (different time-reversal parity).
  - keys of the top group with the irrep and component of the group; n is the multiplicity of the harmonics (X,l,Gamma).
  - "-1" in multiplicity represents no multiplicity.
  - each component is normalized on the ns sites (bonds).
  - multipoles are sorted in order of [Gamma, l, n, p] ("site"), and of Gamma, then [l, n, p] for Q and
    [l, k, T/M, n, p] for T/M ("bond", Q before T/M in each Gamma).
"""


# ==================================================
def _worker(arg):
    gid, pg, pg_id, reference = arg
    if "group" not in _group:
        _group["group"] = load_group()
    g = _group["group"]
    if reference is None:
        return gid, create_cluster_samb(
            pg, g[pg_id]["symmetry_operation"]["tag"], g[gid]["symmetry_operation"]["tag"], g[gid]["wyckoff"]
        )
    point_ops = [np.array(R, dtype=float)[:3, :3] for R in g[gid]["symmetry_operation"]["fractional"]]
    only = {kind: list(refs) for kind, refs in reference.items()}
    return gid, create_cluster_samb(
        pg,
        g[pg_id]["symmetry_operation"]["tag"],
        g[gid]["symmetry_operation"]["tag"],
        g[gid]["wyckoff"],
        reference=reference,
        point_ops=point_ops,
        only=only,
    )


_group = {}


# ==================================================
class _Serial:
    """stand-in for Pool in a single process."""

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def imap_unordered(self, f, args):
        return map(f, args)


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
        - (BinaryManager) -- {PG_id/SG_id: {"site"/"bond": {wyckoff: {PGMultipoleType: samb}}}}.
    """
    args = [a + (None,) for a in group_args()]
    with Pool(os.cpu_count()) if not serial else _Serial() as pool:
        results = dict(pool.imap_unordered(_worker, args))
        # (A): reference sites (bonds) for the images under the inversion coset.
        group = load_group()
        corr = orbit_correspondence(group, inversion_pairs(group, load_info()), None if serial else pool)
        reference = cluster_samb_reference(group, results, corr)
        args2 = [(gid, pg, pg_id, reference[gid]) for gid, pg, pg_id, _ in args if gid in reference]
        for gid, res in pool.imap_unordered(_worker, args2):
            for kind, d in res.items():
                results[gid][kind].update(d)
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    for gid, _, _, _ in args:
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
