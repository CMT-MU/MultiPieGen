# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# ==================================================
# regular atomic orbital basis order (same as tesseral functions up to rank 3).
atomic_orbital_basis = [
    ["s"],  # C0,0
    ["px", "py", "pz"],  # C1,1, S1,1, C1,0
    ["dv", "dxy", "dxz", "dyz", "du"],  # C2,2, S2,2, C2,1, S2,1, C2,0
    ["f2", "f1", "fbz", "f3", "f3x", "f3y", "faz"],  # C3,3, S3,3, C3,2, S3,2, C3,1, S3,1, C3,0
]

# ==================================================
# regular atomic spin basis order.
atomic_spin_basis = ["u", "d"]  # up, down

# ==================================================
# regular atomic LM basis order (M: descending).
atomic_lm_basis = [["0"], ["1", "0", "-1"], ["2", "1", "0", "-1", "-2"], ["3", "2", "1", "0", "-1", "-2", "-3"]]

# ==================================================
# regular atomic JM basis order (J: ascending, M: descending).
atomic_jm_basis = [
    ["(1/2,1/2)", "(1/2,-1/2)"],
    ["(1/2,1/2)", "(1/2,-1/2)", "(3/2,3/2)", "(3/2,1/2)", "(3/2,-1/2)", "(3/2,-3/2)"],
    [
        "(3/2,3/2)",
        "(3/2,1/2)",
        "(3/2,-1/2)",
        "(3/2,-3/2)",
        "(5/2,5/2)",
        "(5/2,3/2)",
        "(5/2,1/2)",
        "(5/2,-1/2)",
        "(5/2,-3/2)",
        "(5/2,-5/2)",
    ],
    [
        "(5/2,5/2)",
        "(5/2,3/2)",
        "(5/2,1/2)",
        "(5/2,-1/2)",
        "(5/2,-3/2)",
        "(5/2,-5/2)",
        "(7/2,7/2)",
        "(7/2,5/2)",
        "(7/2,3/2)",
        "(7/2,1/2)",
        "(7/2,-1/2)",
        "(7/2,-3/2)",
        "(7/2,-5/2)",
        "(7/2,-7/2)",
    ],
]

# ==================================================
# regular atomic spinful basis order.
atomic_lms_basis = [sum([[f"({m},{s})" for s in atomic_spin_basis] for m in rank], []) for rank in atomic_lm_basis]
atomic_lgs_basis = [sum([[f"({m},{s})" for s in atomic_spin_basis] for m in rank], []) for rank in atomic_orbital_basis]

# ==================================================
# bra-ket rank block, (L1,L2) for matrix element of atomic SAMB (uppper triangular block only).
atomic_rank_block = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]

# ==================================================
# convert to regular atomic latex name.
atomic_latex_name = {
    # spin
    "u": r"\uparrow",
    "d": r"\downarrow",
    # orbital rank 0
    "s": "s",
    # orbital rank 1
    "px": "p_x",
    "py": "p_y",
    "pz": "p_z",
    # orbital rank 2
    "dv": "d_v",
    "dxy": "d_{xy}",
    "dxz": "d_{xz}",
    "dyz": "d_{yz}",
    "du": "d_u",
    # orbital rank 3
    "f2": "f_2",
    "f1": "f_1",
    "fbz": r"f_z^\beta",
    "f3": "f_3",
    "f3x": "f_{3x}",
    "f3y": "f_{3y}",
    "faz": r"f_z^\alpha",
}

# ==================================================
# alias for Wannier90: wannier90 -> regular atomic orbital.
wannier90_orbital = {
    "dx2-y2": "dv",
    "dz2": "du",
    "fx(x2-3y2)": "f2",
    "fy(3x2-y2)": "f1",
    "fz(x2-y2)": "fbz",
    "fxyz": "f3",
    "fxz2": "f3x",
    "fyz2": "f3y",
    "fz3": "faz",
}

# ==================================================
# cubic basis order (multipolar PG harmonics internal basis).
cubic_orbital_basis = [
    ["s"],
    ["px", "py", "pz"],
    ["du", "dv", "dyz", "dxz", "dxy"],
    ["f3", "fax", "fay", "faz", "fbx", "fby", "fbz"],
]

# ==================================================
# hexagonal basis order (multipolar PG harmonics internal basis).
hexagonal_orbital_basis = [
    ["s"],
    ["px", "py", "pz"],
    ["du", "dxz", "dyz", "dxy", "dv"],
    ["faz", "f1", "f2", "f3x", "f3y", "f3", "fbz"],
]

# ==================================================
# cubic and hexagonal spinful basis order (optional).
cubic_basis = [sum([[f"({m},{s})" for s in atomic_spin_basis] for m in rank], []) for rank in cubic_orbital_basis]
hexagonal_basis = [sum([[f"({m},{s})" for s in atomic_spin_basis] for m in rank], []) for rank in hexagonal_orbital_basis]

# ==================================================
# convert to latex name up to for internal basis of multipolar PG harmonics and response tensor (add prefix Q/G/T/M).
orbital_latex_name = {
    # common for orbital(tesseral), cubic, hexagonal.
    "s": "_{0}",
    "px": "_{x}",
    "py": "_{y}",
    "pz": "_{z}",
    "dv": "_{v}",
    "dxy": "_{xy}",
    "dxz": "_{xz}",
    "dyz": "_{yz}",
    "du": "_{u}",
    "f2": "_{3a}",
    "f1": "_{3b}",
    "fbz": r"_{z}^{\beta}",
    "f3": "_{xyz}",
    "f3x": "_{3u}",
    "f3y": "_{3v}",
    "faz": r"_{z}^{\alpha}",
    # for cubic.
    "fax": r"_{x}^{\alpha}",
    "fay": r"_{y}^{\alpha}",
    "fbx": r"_{x}^{\beta}",
    "fby": r"_{y}^{\beta}",
    "g": "_{4}",
    "gu": "_{4u}",
    "gv": "_{4v}",
    "gax": r"_{4x}^{\alpha}",
    "gay": r"_{4y}^{\alpha}",
    "gaz": r"_{4z}^{\alpha}",
    "gbx": r"_{4x}^{\beta}",
    "gby": r"_{4y}^{\beta}",
    "gbz": r"_{4z}^{\beta}",
    "hu": "_{5u}",
    "hv": "_{5v}",
    "i": "_{6a}",
    "l": "_{9}",
    # for hexagonal.
    "g0": "_{40}",
    "ga": "_{4b}",
    "gb": "_{4a}",
    "gc": r"_{4u}^{\beta1}",
    "gau": r"_{4v}^{\alpha}",
    "gav": r"_{4u}^{\alpha}",
    "i0": "_{6b}",
    "j": "_{7}",
}

# ==================================================
# convert from name to key, ("tesseral", (l,m)), ("Oh", (rank, irrep, mul, comp)), ("D6h", (rank, irrep, mul, comp)).
# negative m represents S, otherwise C.
harmonics_key = {
    # common for orbital(tesseral), cubic, hexagonal.
    "s": ("tesseral", (0, 0)),
    "px": ("tesseral", (1, 1)),
    "py": ("tesseral", (1, -1)),
    "pz": ("tesseral", (1, 0)),
    "dv": ("tesseral", (2, 2)),
    "dxy": ("tesseral", (2, -2)),
    "dxz": ("tesseral", (2, 1)),
    "dyz": ("tesseral", (2, -1)),
    "du": ("tesseral", (2, 0)),
    "f2": ("tesseral", (3, 3)),
    "f1": ("tesseral", (3, -3)),
    "fbz": ("tesseral", (3, 2)),
    "f3": ("tesseral", (3, -2)),
    "f3x": ("tesseral", (3, 1)),
    "f3y": ("tesseral", (3, -1)),
    "faz": ("tesseral", (3, 0)),
    # for cubic.
    "fax": ("Oh", (3, "T1u", -1, 0)),
    "fay": ("Oh", (3, "T1u", -1, 1)),
    "fbx": ("Oh", (3, "T2u", -1, 0)),
    "fby": ("Oh", (3, "T2u", -1, 1)),
    "g": ("Oh", (4, "A1g", -1, -1)),
    "gu": ("Oh", (4, "Eg", -1, 0)),
    "gv": ("Oh", (4, "Eg", -1, 1)),  # -gv
    "gax": ("Oh", (4, "T1g", -1, 0)),
    "gay": ("Oh", (4, "T1g", -1, 1)),
    "gaz": ("Oh", (4, "T1g", -1, 2)),
    "gbx": ("Oh", (4, "T2g", -1, 0)),
    "gby": ("Oh", (4, "T2g", -1, 1)),
    "gbz": ("Oh", (4, "T2g", -1, 2)),
    "hu": ("Oh", (5, "Eu", -1, 0)),
    "hv": ("Oh", (5, "Eu", -1, 1)),
    "i": ("Oh", (6, "A2g", -1, -1)),
    "l": ("Oh", (9, "A1u", -1, -1)),
    # for hexagonal.
    "g0": ("D6h", (4, "A1g", -1, -1)),
    "ga": ("D6h", (4, "B1g", -1, -1)),
    "gb": ("D6h", (4, "B2g", -1, -1)),
    "gc": ("D6h", (4, "E2g", 1, 0)),
    "gau": ("D6h", (4, "E1g", -1, 0)),
    "gav": ("D6h", (4, "E1g", -1, 1)),
    "i0": ("D6h", (6, "A2g", -1, -1)),
    "j": ("D6h", (7, "A1u", -1, -1)),
}
