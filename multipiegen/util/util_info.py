# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Load overall information of groups (multipiegen/binary/info.pkl).

The binary is created by multipiegen/create_binary/create_info.py.
"""

from multipie.util.util_binary import BinaryManager
from multipiegen import __bin_dir__

info_file = "info"


# ==================================================
def load_info():
    """
    Load info.

    Returns:
        - (BinaryManager) -- info (dict-like). str(info) gives the description of the contents.
    """
    return BinaryManager(info_file, topdir=__bin_dir__)
