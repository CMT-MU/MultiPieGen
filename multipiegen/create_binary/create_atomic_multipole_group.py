# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create atomic_multipole_group.pkl (atomic multipole SAMB of all 47 point groups, needs atomic_multipole.pkl).

    python -m multipiegen.create_binary.create_atomic_multipole_group [--serial]

See multipiegen/util/util_atomic_multipole.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.data.data_basis import atomic_rank_block
from multipiegen.util.util_harmonics_mapping import data_harmonics_mapping
from multipiegen.util.util_atomic_multipole import (
    atomic_multipole_group_file,
    create_atomic_multipole_pg,
    create_top_samb,
    load_atomic_multipole,
    max_rank,
)
from multipiegen.util.util_info import load_info

h_atomic_multipole_group = f"""
* Matrix of atomic multipole SAMB for all PG (rank l <= {max_rank}, s = 0, 1, orbitals of L = 0..3).
- PG_id (str): (dict) matrix data.
  - "jml/lgs/lg" (str): (dict) matrix of atomic multipole in (J,M;L)/(L,gamma,s)/(L,gamma) basis.
    - (L1,L2) (int,int): (dict) atomic multipole data for each bra-ket block.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,2(2L1+1),2(2L2+1),sympy)) matrix for each component.
NOTE:
  - the expression of the harmonics (X,l,Gamma,n) of a key is in harmonics.pkl.
  - (L,gamma) basis: orbitals of info["harmonics"]["atomic_basis"]["spinless"]["lg"].
  - half size of matrix for spinless "lg" (s = 0, up-spin part of "lgs" times sqrt(2)).
  - each component is normalized to the Frobenius norm 1 (diagonal block) or 1/sqrt(2) (off-diagonal block).
  - n is the multiplicity of the harmonics (X,l,Gamma) (its xyz polynomial, not renumbered); p = -1; x = "q".
  - multipoles are sorted in order of [Q/G/T/M, s, k, l, Gamma, n, p].
"""

_am = {}


# ==================================================
def _top_worker(arg):
    if "am" not in _am:
        _am["am"] = load_atomic_multipole()
    top, block, bt = arg
    return arg, create_top_samb(_am["am"], top, block, bt)


# ==================================================
def create_atomic_multipole_group(serial=False, verbose=True):
    """
    Create atomic multipole SAMB of all point groups (BinaryManager).

    Args:
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- {PG_id: {"jml"/"lgs"/"lg": {(L1, L2): {PGMultipoleType: matrix}}}}.
    """
    args = [(top, block, bt) for top in ["Oh", "D6h"] for block in atomic_rank_block for bt in ["jml", "lgs"]]
    if serial:
        results = dict(map(_top_worker, args))
    else:
        with Pool(os.cpu_count()) as pool:
            results = dict(pool.imap_unordered(_top_worker, args))
    top_data = {}
    for (top, block, bt), d in results.items():
        top_data.setdefault(top, {}).setdefault(bt, {})[block] = d

    info = load_info()
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    for gid in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        pg = info["tag"][gid]
        bm[gid] = create_atomic_multipole_pg(pg, top_data[data_harmonics_mapping[pg]["top"]])
    bm.add_comment(h_atomic_multipole_group)
    return bm


# ==================================================
def save_atomic_multipole_group(serial=False):
    """
    Create and save atomic_multipole_group.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_atomic_multipole_group(serial=serial, verbose=False).save_binary(atomic_multipole_group_file)
    print(f"atomic_multipole_group.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create atomic_multipole_group.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_atomic_multipole_group(parser.parse_args().serial)
