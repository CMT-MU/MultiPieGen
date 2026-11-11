# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Example: atomic multipole matrices and the atomic multipole SAMB of a point group.

    python docs/examples/atomic_multipole.py

The binary data (multipiegen/binary/atomic_multipole.pkl, atomic_multipole_group.pkl) are created by
python -m multipiegen.create_binary.create_atomic_multipole and create_atomic_multipole_group.
"""

import numpy as np

from multipiegen.data.data_basis import atomic_lgs_basis
from multipiegen.util.util_atomic_multipole import load_atomic_multipole, load_atomic_multipole_group
from multipiegen.util.util_harmonics import harmonics_expression
from multipiegen.util.util_info import load_info

# ==================================================
# atomic multipole: {"lms"/"jml"/"lm": {(L1, L2): {SphericalMultipoleType(X, l, s, k, "q"): ndarray(2l+1, n1, n2)}}}.
am = load_atomic_multipole()
for key, mat in am["lms"][(1, 1)].items():
    print(key, mat.shape)

# ==================================================
# atomic multipole SAMB of Oh: {"jml"/"lgs"/"lg": {(L1, L2): {PGMultipoleType: matrix (dim, n1, n2)}}}.
# the expression of the harmonics (X, l, Gamma, n) is looked up by the key in harmonics.pkl.
info = load_info()
pg = load_atomic_multipole_group()[info["id"]["Oh"]]
print("lgs basis of p:", atomic_lgs_basis[1])
for key, mat in list(pg["lgs"][(1, 1)].items())[:3]:
    print(key, "symmetry:", list(harmonics_expression("Oh", key)))
    print(mat[0])

# orthonormality: Re tr(A^dagger B) = delta (diagonal block) or delta/2 (off-diagonal block).
for block in [(1, 1), (1, 2)]:
    v = np.array([np.asarray(c, dtype=complex).ravel() for m in pg["lgs"][block].values() for c in m])
    g = (v.conj() @ v.T).real
    print(block, len(v), "matrices, orthonormal:", np.allclose(g, np.eye(len(v)) * (1 if block[0] == block[1] else 0.5)))
