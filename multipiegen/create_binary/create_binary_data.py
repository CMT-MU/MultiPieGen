# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create binary_data for MultiPie (multipiegen/binary_data), in the layout of MultiPie's multipie/binary_data.
    - info.pkl
    - {PG,SG,MPG,MSG}/{id}.pkl (one file for each group, ":" of the id replaced by "_", e.g. PG/PG_32.pkl)

Needs the binaries in multipiegen/binary (info, group, harmonics, multipolar_harmonics, atomic_multipole_group,
cluster_samb, response_tensor); run their create_binary modules first.

    python -m multipiegen.create_binary.create_binary_data [--serial]

See multipiegen/util/util_binary_data.py.
"""

import argparse
import os
import shutil
import time
from multiprocessing import Pool

from multipiegen.util.util_binary_data import (
    binary_data_dir,
    create_group_data,
    create_info_data,
    group_file,
    group_ids,
    group_types,
)


# ==================================================
def _worker(gid):
    create_group_data(gid).save_binary(group_file(gid))
    return gid


# ==================================================
def save_binary_data(serial=False):
    """
    Create and save binary_data (info.pkl and the files of all groups) in multipiegen/binary_data.

    Args:
        serial (bool, optional): run in a single process ?

    Note:
        - the subdirectories PG, SG, MPG, MSG are removed first, so that no file of an old group remains.
    """
    t0 = time.time()
    for tp in group_types:
        shutil.rmtree(os.path.join(binary_data_dir, tp), ignore_errors=True)
    create_info_data().save_binary("info")
    ids = group_ids()
    if serial:
        done = list(map(_worker, ids))
    else:
        with Pool(os.cpu_count()) as pool:
            done = list(pool.imap_unordered(_worker, ids, chunksize=16))
    print(f"binary_data: info + {len(done)} groups, {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create binary_data for MultiPie")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_binary_data(parser.parse_args().serial)
