# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create atomic_multipole.pkl (matrices of atomic multipoles in the LMs/JML/LM basis, L <= 3).

    python -m multipiegen.create_binary.create_atomic_multipole [--serial]

See multipiegen/util/util_atomic_multipole.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util.util_atomic_multipole import (
    assemble_atomic_multipole,
    atomic_multipole_file,
    atomic_multipole_keys,
    atomic_multipole_matrix,
    max_rank,
)

h_atomic_multipole = f"""
* Matrix of atomic multipole (rank l <= {max_rank}, s = 0, 1) in the orbitals of L = 0..3.
- "lms/jml/lm" (str): (dict) matrix in LMs/JML/LM basis (info["harmonics"]["atomic_basis"]).
  - (L1,L2) (int,int): (dict) atomic multipole data for each bra-ket block.
    - (X,l,s,k,x) (SphericalMultipoleType): (ndarray(2l+1,2(2L1+1),2(2L2+1),sympy)) matrix for m = l, l-1, ..., -l.
NOTE:
  - half size of matrix for spinless "lm" (s = 0, up-spin part of "lms").
  - x = "q" for all.
  - multipoles are sorted in order of [Q/G/T/M, s, k, l].
"""


# ==================================================
def _worker(arg):
    return arg, atomic_multipole_matrix(*arg)


# ==================================================
def create_atomic_multipole(serial=False, verbose=True):
    """
    Create atomic multipole (BinaryManager).

    Args:
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- {"lms"/"jml"/"lm": {(L1, L2): {SphericalMultipoleType: matrix}}}.
    """
    args = [(bt, *key) for bt in ["lms", "jml"] for key in atomic_multipole_keys()]
    if serial:
        results = dict(map(_worker, args))
    else:
        with Pool(os.cpu_count()) as pool:
            results = dict(pool.imap_unordered(_worker, args))
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    for k, v in assemble_atomic_multipole(results).items():
        bm[k] = v
    bm.add_comment(h_atomic_multipole)
    return bm


# ==================================================
def save_atomic_multipole(serial=False):
    """
    Create and save atomic_multipole.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_atomic_multipole(serial=serial, verbose=False).save_binary(atomic_multipole_file)
    print(f"atomic_multipole.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create atomic_multipole.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_atomic_multipole(parser.parse_args().serial)
