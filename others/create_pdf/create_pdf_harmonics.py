# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create LaTeX and PDF of harmonics with polar and axial symbols side by side, for all 47 point groups.
    - others/pdf/PG/{no}-{tag}/harmonics.{tex,pdf}

The expression is that of the polar harmonics. The axial symbol is obtained from the polar one
by changing irrep, interchanging components, and adding a sign (compatibility relation).
"""

import os
from multipie import Group
from multipie.util.util_pdf_latex import PDFviaLaTeX
from multipiegen import __top_dir__
from multipiegen.data.data_group_pg import pg_id_set, pg_tag
from multipiegen.util.util_harmonics_symbol import harmonics_symbol_table

pdfdir = os.path.join(__top_dir__, "others", "pdf")
max_rank = 11


# ==================================================
def group_dir(tag):
    """
    Output directory of a point group, e.g., others/pdf/PG/032-Oh.
    """
    group = Group(tag)
    no = int(group.info.PG.split(":")[1])
    return os.path.join(pdfdir, "PG", f"{no:03d}-{tag}")


# ==================================================
def create_pdf_harmonics(tag):
    """
    Create harmonics.tex and harmonics.pdf of a point group.

    Args:
        tag (str): point group tag.
    """
    group = Group(tag)
    h_dir = group_dir(tag)
    os.makedirs(h_dir, exist_ok=True)

    pdf = PDFviaLaTeX("harmonics", dir=h_dir, pt=8, style="narrowest", landscape=True, english=True)
    pdf.title(group.latex(True))

    rows = harmonics_symbol_table(tag)
    col = ["polar", "axial", "expression"]
    no = 1
    for l in range(max_rank + 1):
        tbl = []
        row = []
        for key, p_sym, a_sym, expr in rows:
            if key[0] != l:
                continue
            tbl.append([p_sym, a_sym, expr])
            row.append(str(no))
            no += 1
        cap = f"Harmonics for rank {l}."
        pdf.table(tbl, row, col, "No.", cpos="cccl", caption=cap, stretch=1.3, rmath=True, tmath=True)

    pdf.build()
    print(f"created {os.path.join(h_dir, 'harmonics.pdf')}")


# ==================================================
def create_pdf_harmonics_all():
    """
    Create harmonics PDF of all 47 point groups (32 crystallographic + 5 settings + 10 complex).
    """
    for id_s in pg_id_set["all"] + pg_id_set["complex"]:
        create_pdf_harmonics(pg_tag[id_s])


# ==================================================
if __name__ == "__main__":
    create_pdf_harmonics_all()
