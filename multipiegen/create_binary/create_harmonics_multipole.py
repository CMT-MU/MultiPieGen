# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create harmonics_multipole.pkl (multipolar harmonics, s = 1, 2, 3, up to rank 3, of all 47 point groups).

    python -m multipiegen.create_binary.create_harmonics_multipole [--serial]

See multipiegen/util/util_harmonics_multipole.py.
"""

import argparse
import os
import time
from multiprocessing import Pool

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.util.util_harmonics_multipole import (
    create_harmonics_multipole,
    harmonics_multipole_file,
    max_rank,
    max_s_rank,
)
from multipiegen.util.util_info import load_info

h_harmonics_multipole = f"""
* Multipolar harmonics of all PG (s = 1..{max_s_rank}, rank l <= {max_rank}) with polar (q) and axial (g) internal basis.
- PG_id (str): (dict) multipolar harmonics.
  - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,sympy), ndarray(dim,2s+1,sympy), ndarray(dim,sympy))
    [cartesian ex.] of the harmonics (X,l,Gamma), [coefficient of each internal basis], and [multipole ex.].
NOTE:
  - internal basis: orbitals of data_basis.cubic_orbital_basis[s] (cubic) or hexagonal_orbital_basis[s] (trigonal,
    hexagonal), with the symbols Q/G + orbital_latex_name (info["harmonics"]["internal_basis"]).
  - x = q if (X = Q and s+k even) or (X = G and s+k odd), otherwise g.
  - multiplicity n is renumbered among the same (x, l, X, Gamma, k); "-1" represents no multiplicity; p = -1.
  - sorted in order of [q/g, s, Gamma, l, Q/G, k, n].
"""


# ==================================================
def _worker(arg):
    gid, tag = arg
    return gid, create_harmonics_multipole(tag)


# ==================================================
def create_harmonics_multipole_all(serial=False, verbose=True):
    """
    Create multipolar harmonics of all point groups (BinaryManager).

    Args:
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- {PG_id: {PGMultipoleType: (ex, mharm, m_ex)}}.
    """
    info = load_info()
    ids = info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]
    args = [(i, info["tag"][i]) for i in ids]
    bm = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    if serial:
        results = map(_worker, args)
        for gid, d in results:
            bm[gid] = d
    else:
        with Pool(os.cpu_count()) as pool:
            for gid, d in pool.imap(_worker, args):
                bm[gid] = d
    bm.add_comment(h_harmonics_multipole)
    return bm


# ==================================================
def save_harmonics_multipole(serial=False):
    """
    Create and save harmonics_multipole.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_harmonics_multipole_all(serial=serial, verbose=False).save_binary(harmonics_multipole_file)
    print(f"harmonics_multipole.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create harmonics_multipole.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_harmonics_multipole(parser.parse_args().serial)
