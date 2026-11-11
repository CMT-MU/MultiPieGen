# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create harmonics.pkl (expression and u-matrix of the harmonics of all 47 point groups, rank <= 11).

    python -m multipiegen.create_binary.create_harmonics [--serial]

See multipiegen/util/util_harmonics.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util.util_harmonics import create_harmonics, harmonics_file, max_rank
from multipiegen.util.util_info import load_info

h_harmonics = f"""
* Harmonics of all PG (rank l <= {max_rank}): cartesian expression and u-matrix for each multipole key.
- PG_id (str): (dict) harmonics.
  - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,sympy), ndarray(2l+1,dim,sympy)) [expression] and [u-matrix].
NOTE:
  - X = Q/G, p = -1, s = k = 0, x = "q" in the key. Look up the harmonics of a key of the other binaries by its
    (X,l,Gamma,n), with T -> Q and M -> G (multipiegen.util.util_harmonics.harmonics_expression, harmonics_u_matrix).
  - n is the multiplicity of the harmonics (X,l,Gamma), 1-based; "-1" represents no multiplicity.
  - expression: polynomial in x, y, z of each component; u-matrix: <m|component> for m = l, l-1, ..., -l.
  - the harmonics is c times a top harmonics of Oh or D6h (data_harmonics_mapping).
  - sorted in order of [Gamma, l, Q/G, n].
"""


# ==================================================
def _worker(arg):
    gid, tag = arg
    return gid, create_harmonics(tag)


# ==================================================
def create_harmonics_all(serial=False, verbose=True):
    """
    Create harmonics of all point groups (BinaryManager).

    Args:
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- {PG_id: {PGMultipoleType: (expression, u_matrix)}}.
    """
    info = load_info()
    args = [(gid, info["tag"][gid]) for gid in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]]
    if serial:
        results = dict(map(_worker, args))
    else:
        with Pool(min(os.cpu_count(), 8)) as pool:
            results = dict(pool.imap_unordered(_worker, args))
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    for gid, _ in args:
        bm[gid] = results[gid]
    bm.add_comment(h_harmonics)
    return bm


# ==================================================
def save_harmonics(serial=False):
    """
    Create and save harmonics.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_harmonics_all(serial=serial, verbose=False).save_binary(harmonics_file)
    print(f"harmonics.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create harmonics.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_harmonics(parser.parse_args().serial)
