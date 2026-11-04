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
# convert to cubic latex name up to rank 4 (+alpha) for multipolar PG harmonics and response tensor (add prefix Q/G/T/M).
cubic_latex_name = {
    # rank 0
    "s": "_0",
    # rank 1
    "px": "_x",
    "py": "_y",
    "pz": "_z",
    # rank 2
    "du": "_u",
    "dv": "_v",
    "dyz": "_{yz}",
    "dxz": "_{xz}",
    "dxy": "_{xy}",
    # rank 3
    "f3": "_{xyz}",
    "fax": r"_x^\alpha",
    "fay": r"_y^\alpha",
    "faz": r"_z^\alpha",
    "fbx": r"_x^\beta",
    "fby": r"_y^\beta",
    "fbz": r"_z^\beta",
    # rank 4
    "g": "_4",
    "gu": "_{4u}",
    "gv": "_{4v}",
    "gax": r"_{4x}^\alpha",
    "gay": r"_{4y}^\alpha",
    "gaz": r"_{4z}^\alpha",
    "gbx": r"_{4x}^\beta",
    "gby": r"_{4y}^\beta",
    "gbz": r"_{4z}^\beta",
    # rank 5
    "hu": "_{5u}",
    "hv": "_{5v}",
    # rank 6
    "i": "_{6a}",
    # rank 9
    "l": "_9",
}

# ==================================================
# convert to hexagonal latex name up to rank 4 (+alpha) for multipolar PG harmonics and response tensor.
hexagonal_latex_name = {
    # rank 0
    "s": "_0",
    # rank 1
    "px": "_x",
    "py": "_y",
    "pz": "_z",
    # rank 2
    "du": "_u",
    "dv": "_v",
    "dyz": "_{yz}",
    "dxz": "_{xz}",
    "dxy": "_{xy}",
    # rank 3
    "f1": "_{3b}",
    "f2": "_{3a}",
    "f3": "_{xyz}",
    "f3x": "_{3u}",
    "f3y": "_{3v}",
    "faz": r"_z^\alpha",
    "fbz": r"_z^\beta",
    # rank 4
    "g0": "_{40}",
    "ga": "_{4b}",
    "gb": "_{4a}",
    "gc": r"_{4u}^{\beta1}",
    "gv": r"_{4u}^{\beta2}",
    "gau": r"_{4v}^\alpha",
    "gav": r"_{4u}^\alpha",
    "gaz": r"_{4v}^{\beta1}",
    "gbz": r"_{4v}^{\beta2}",
    # rank 6
    "i0": "_{6b}",
    # rank 7
    "j": "_7",
}
