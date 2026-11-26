# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create coset.pkl (cosets of each real PG and SG in its parents). Needs group.pkl and cluster_samb.pkl.

    python -m multipiegen.create_binary.create_coset [--serial]

See multipiegen/util/util_coset.py and multipiegen/util/util_cluster_samb_image.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

import numpy as np
import sympy as sp
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util import util_coset
from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_cluster_samb_image import cluster_samb_transport, inversion_pairs, orbit_correspondence
from multipiegen.util.util_coset import coset_file
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

h_coset = """
* Cosets of all real PG and SG in their parents.
- PG_id/SG_id (str): (dict) coset data of a group H (only for a group with a parent).
  - "parent" (str): ([str]) parents along the chain of representative parents, the direct parent first.
  - "coset" (str): (dict) left cosets g H in each parent.
    - parent_id (str): ([[str]]) operations of each coset, its representative first, without H itself and the cosets of
      translations only (PG: SO tag of the parent; SG: "point part of the SO tag of the parent:[translation in the
      coordinates of H modulo its lattice]").
  - "first" (str): (dict) first coset of H in the direct parent (the first one other than H).
    - "parent" (str): (str) direct parent.
    - "R" (str): (Matrix(sympy)) x_parent = R x_H (data_parent_pg/sg).
    - "index" (str): (int) number of cosets [P:H].
    - "ops" (str): ([str]) SO tags of the parent in the coset.
    - "rep" (str): (str) SO tag of the representative g (first operation).
    - "rep_matrix" (str): (Matrix(sympy)) g in the coordinates of H (4x4 for SG, 3x3 for PG).
    - "order" (str): (int) smallest n with g^n in H (modulo the lattice of H).
    - "kind" (str): (str) "rotation" or "inversion".
    - "frame" (str): (Matrix(sympy)) frame of the image system (g in the coordinates of H) for "rotation", None otherwise.
    - "inversion" (str): (dict) for "inversion", None otherwise.
      - "image" (str): (str) image group H' = g H g^-1.
      - "map" (str): (Matrix(sympy)) M: coordinates of H -> those of H'.
      - "shift" (str): (Matrix(sympy)) origin shift k (coordinates of the parent) with k H' k^-1 = g H g^-1 (SG only).
      - "wyckoff" (str): (dict) W: (W', {parameter of W': expression of the parameters of W}).
      - "site"/"bond" (str): (dict) image of each site/bond Wyckoff position.
        - W (str): (dict) "image": (H', W'), "perm": row i of W -> row perm[i] of W', "sign": -1 for a reversed bond,
          "diag": the cluster SAMB are mapped label by label, "eps": {PGMultipoleType: [+-1 for each component]} (diag)
          or "U": {block: (rows of W, rows of W', U)} (otherwise), rows = [(PGMultipoleType, component)], T z_a = sum_b
          U_ab z'_b within each block (irrep, Q and T/M separately for bond).
    - "parent_limit" (str): (dict) parent structures, W: [{"param": {parameter: value}, "position": representative
      position, "parent_wyckoff": Wyckoff position of the parent, "preserved_by": indices of the cosets (in the order of
      the coset images) that keep it}].
NOTE:
  - "-c" PG are not included (they share the cosets of the real PG).
  - the cluster SAMB of the image of the inversion are those of cluster_samb.pkl, after the change of the reference
    sites (bonds) (util_cluster_samb_image).
"""


# ==================================================
def _exact(M):
    """float matrix -> exact (rational, or with sqrt(2), sqrt(3))."""
    return sp.Matrix(np.asarray(M, dtype=float)).applyfunc(
        lambda v: sp.nsimplify(round(float(v), 12), [sp.sqrt(2), sp.sqrt(3)], tolerance=1e-9)
    )


# ==================================================
def coset_data(gid, transport=None):
    """
    Coset data of a group.

    Args:
        gid (str): group id ("PG:n" of a real PG, or "SG:n").
        transport (dict, optional): cluster_samb_transport of the group ({"site", "bond"}), for "inversion".

    Returns:
        - (dict) -- {"parent", "coset", "first"}, None for a group without a parent.
    """
    chain = util_coset.parent_chain(gid)
    if not chain:
        return None
    sg = gid.startswith("SG:")
    if sg:
        f = util_coset.first_coset_sg(int(gid[3:]))
        limit = util_coset.parent_limit_sg(int(gid[3:]))
        to_id = lambda g: f"SG:{g}"
    else:
        tag = util_coset._pg_tag(gid)
        f = util_coset.first_coset_pg(tag)
        limit = util_coset.parent_limit_pg(tag)
        to_id = lambda g: util_coset.pg_id[g]
    rep = _exact(f["rep_child"])
    inversion = None
    if f["kind"] == "inversion":
        inversion = {"image": to_id(f["image"]), "map": sp.Matrix(str_to_sympy(f["map"])), "wyckoff": f["wyckoff"]}
        if sg:
            inversion["shift"] = sp.Matrix(str_to_sympy(f["shift"]))
        inversion.update(transport)
    first = {
        "parent": to_id(f["parent"]),
        "R": sp.Matrix(str_to_sympy(f["R"])),
        "index": f["index"],
        "ops": f["ops"],
        "rep": f["rep"],
        "rep_matrix": rep,
        "order": f["order"],
        "kind": f["kind"],
        "frame": rep if f["kind"] == "rotation" else None,
        "inversion": inversion,
        "parent_limit": limit,
    }
    return {
        "parent": [p for p, _ in chain],
        "coset": {p: util_coset.coset_list(gid, R, p) for p, R in chain},
        "first": first,
    }


# ==================================================
def _worker(arg):
    return arg[0], coset_data(*arg)


# ==================================================
def create_coset_all(serial=False, verbose=True):
    """
    Create coset data of all real PG and SG (BinaryManager).

    Args:
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- {PG_id/SG_id: {"parent", "coset", "first"}}.
    """
    info, group = load_info(), load_group()
    ids = info["id_set"]["PG"]["all"] + info["id_set"]["SG"]["all"]
    pool = None if serial else Pool(os.cpu_count())
    try:
        corr = orbit_correspondence(group, inversion_pairs(group, info), pool)
        transport = cluster_samb_transport(load_cluster_samb(), corr)
        args = [(gid, transport.get(gid)) for gid in ids]
        results = dict(map(_worker, args) if pool is None else pool.imap_unordered(_worker, args))
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    for gid in ids:
        if results[gid] is not None:
            bm[gid] = results[gid]
    bm.add_comment(h_coset)
    return bm


# ==================================================
def save_coset(serial=False):
    """
    Create and save coset.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_coset_all(serial=serial, verbose=False).save_binary(coset_file)
    print(f"coset.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create coset.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_coset(parser.parse_args().serial)
