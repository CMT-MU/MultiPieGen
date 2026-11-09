# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Load group data (multipiegen/binary_data/group.pkl).

The binary is created by multipiegen/create_binary/create_group.py.
"""

from multipie.util.util_binary import BinaryManager
from multipiegen import __bin_dir__

group_file = "group"


# ==================================================
def load_group():
    """
    Load group.

    Returns:
        - (BinaryManager) -- group (dict-like), {group_id: data}. str(group) gives the description of the contents.
    """
    return BinaryManager(group_file, topdir=__bin_dir__)
