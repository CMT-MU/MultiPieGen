# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: response tensors of a magnetic point group in terms of active multipoles.

    python docs/examples/response_tensor.py

The binary data (multipiegen/binary/response_tensor.pkl) is created by
multipiegen/create_binary/create_response_tensor.py (python -m multipiegen.create_binary.create_response_tensor),
and the PDFs (multipiegen/pdf/MPG/{id}/response_tensor_X.pdf) by multipiegen/create_pdf/create_pdf_group.py.
"""

from multipiegen.util.util_group import load_group
from multipiegen.util.util_response_tensor import load_response_tensor

rt = load_response_tensor()
group = load_group()

mpg = "MPG:24.5.149"
print(mpg, group[mpg]["info"].international, ", active multipoles:", group[mpg]["active_multipole"])

# ==================================================
# 1. keys: (X, rank, opt), X = Q (polar), G (axial), T (magnetic polar), M (magnetic axial).
d = rt[mpg]
print("number of tensors:", len(d), ", e.g.", list(d)[:4])

# ==================================================
# 2. non-vanishing tensors and their expressions in terms of active multipoles.
for (X, rank, opt), (m, ex) in d.items():
    if ex and rank <= 3:
        print(f"\n{X} rank {rank} {opt}:\n{m}")
        for c, e in ex.items():
            print(f"  {c} = {e}")
