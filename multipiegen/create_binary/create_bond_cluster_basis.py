# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create orthonormal bond cluster basis of Oh and D6h from the vector cluster basis, and save it as binary.
    - multipiegen/binary/bond_cluster_basis.pkl
      (load by multipiegen.util.util_bond_cluster_basis.load_bond_cluster_basis).

Run:  python -m multipiegen.create_binary.create_bond_cluster_basis   (needs cluster_vector_basis.pkl, a few s)
"""

from multipiegen.util.util_bond_cluster_basis import save_bond_cluster_basis

# ==================================================
if __name__ == "__main__":
    save_bond_cluster_basis()
