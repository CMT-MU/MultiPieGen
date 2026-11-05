"""
Create top harmonics information.
    - top_harmonics_symmetry: symmetry (expression) of harmonics.
    - top_harmonics_u_matrix: u-matrix to spherical harmonics.

note:
    - u-matrix <m|gamma> is in descending order, i.e., m=l,l-1,...,-l.
    - axial irrep is obtained by interchanging "g" and "u".
    - representation matrices of all symmetry operations are common for the same irrep. irrespective of polar or axial (checked).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import os
from multipiegen.data.data_top_harmonics import top_harmonics

# ==================================================
header = """
Top harmonics infomation (polar) [axial can be obtained by interchanging g and u in irrep].

- top_harmonics_symmetry.
    { "Oh/D6h": { (rank, irrep, multiplicity, component): symmetry } }.
- top_harmonics_u_matrix, <m|gamma>.
    { "Oh/D6h": { (rank, irrep, multiplicity, component): u-matrix } }.
"""


# ==================================================
def create_symmetry_umat():
    sym = {"Oh": {}, "D6h": {}}
    u_mat = {"Oh": {}, "D6h": {}}
    for k, v in top_harmonics.items():
        top, rank, irrep, mul, comp = k
        t, e, u = v
        key = (rank, irrep, mul, comp)
        sym[top][key] = e.replace("**", "^").replace("*", "").replace(" ", "").replace("^", "**")
        u_mat[top][key] = u.replace("*", "").replace(" ", "")

    return sym, u_mat


# ==================================================
def create_file():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "data_top_harmonics_info.py")

    sym, u_mat = create_symmetry_umat()

    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("top_harmonics_symmetry =", sym, file=f)
        print("", file=f)
        print("top_harmonics_u_matrix =", u_mat, file=f)


# ==================================================
create_file()
