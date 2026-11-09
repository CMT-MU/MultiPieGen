# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: representation matrices of a point group from the harmonics.

    python docs/examples/representation_matrix.py

The PDFs (multipiegen/pdf/PG/{no}-{tag}/representation_matrix.pdf) are created by
multipiegen/create_pdf/create_pdf_group.py (python -m multipiegen.create_pdf.create_pdf_group --group PG).
"""

from multipiegen.util.util_representation_matrix import (
    check_representation_matrix,
    harmonics_blocks,
    representation_matrix,
)

pg = "D3h"

# ==================================================
# 1. harmonics blocks (X, rank, irrep, multiplicity) -> [(coeff, top key)] of each component.
blocks = harmonics_blocks(pg)
print(f"{pg}: {len(blocks)} harmonics blocks up to rank 11, e.g.")
for b in list(blocks)[:4]:
    print("  ", b, blocks[b])

# ==================================================
# 2. all harmonics of the same irrep have the same representation matrices.
print("\nirrep: number of harmonics blocks, max deviation from the first one, max closure residual")
for irrep, d in check_representation_matrix(pg).items():
    print(f"  {irrep:4s}: {len(d['blocks']):3d}, {d['deviation']:.1e}, {d['closure']:.1e}")

# ==================================================
# 3. exact representation matrices (in the order of symmetry operations).
rm = representation_matrix(pg)
D = rm["matrix"]["E'"]
for tag, d in list(zip(rm["tag"], D))[:3]:
    print(f"\nE' {tag}:\n{d}")
