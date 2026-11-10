# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create generated pages of the documentation.
    - docs/src/{pg,sg,mpg,msg}_pdf.md: lists of the PDFs of each group (multipiegen/pdf/{PG,SG,MPG,MSG}), sub-pages of
      "Group Information" (docs/src/group_info.md). These pages are meant to be linked from outside (e.g. MultiPie),
      so they describe the content only, not how it is generated.
    - docs/src/examples.md: code of the examples listed in docs/examples/README.md.

    python docs/create_pages.py    (run by "make build" in docs)
"""

import os
import re

from multipiegen import __top_dir__
from multipiegen.create_pdf.create_pdf_group import group_dir, has_pdf, kinds
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

doc_dir = os.path.dirname(os.path.abspath(__file__))
crystal_name = ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"]
label = {
    "harmonics": "harmonics",
    "character_table": "character",
    "symmetry_operation": "symmetry operation",
    "wyckoff_site": "Wyckoff site",
    "wyckoff_bond": "Wyckoff bond",
    "representation_matrix": "representation matrix",
    "response_tensor_Q": "Q tensor",
    "multipolar_harmonics_s1_q": "s=1 q",
    "multipolar_harmonics_s1_g": "s=1 g",
    "multipolar_harmonics_s2_q": "s=2 q",
    "multipolar_harmonics_s2_g": "s=2 g",
    "multipolar_harmonics_s3_q": "s=3 q",
    "multipolar_harmonics_s3_g": "s=3 g",
    "atomic_multipole_lg": "atomic lg",
    "atomic_multipole_lgs": "atomic lgs",
    "atomic_multipole_jml": "atomic jml",
    "response_tensor_G": "G tensor",
    "response_tensor_T": "T tensor",
    "response_tensor_M": "M tensor",
}

head = {
    "PG": """# Point groups

All 47 point groups: 32 crystallographic point groups, 5 other settings ("-1"), and 10 complex groups ("-c").

- harmonics: harmonics up to rank 11, with the polar (Q) and axial (G) symbols side by side (the expression is that
  of the polar harmonics).
- character: character table, polar-axial conversion, and symmetric and anti-symmetric products.
- symmetry operation, Wyckoff site and bond, representation matrix: as the names say.
- s=1..3 q/g: multipolar harmonics (rank <= 3) with internal rank s and polar (q) or axial (g) internal basis.
- atomic lg/lgs/jml: atomic multipole SAMB (matrices in the bra-ket blocks of s, p, d, f orbitals) in the
  (L,gamma), (L,gamma,s), and (J,M;L) basis.
""",
    "SG": """# Space groups

All 230 space groups: symmetry operations, Wyckoff sites, and Wyckoff bonds.

The numbers, Hermann-Mauguin symbols, symmetry operations, and Wyckoff positions are consistent with the
International Tables for Crystallography, Vol. A (ITA) and the Bilbao Crystallographic Server. Only the standard
setting of each space group is given (the default setting of the Bilbao Crystallographic Server):

- monoclinic groups: unique axis b, cell choice 1.
- groups with two origin choices: origin choice 2 (inversion center at the origin).
- rhombohedral groups: hexagonal axes.
""",
    "MPG": """# Magnetic point groups

All 157 magnetic point groups: symmetry operations, Wyckoff sites, and response tensors (Q, G, T, M up to rank 4).
The magnetic tensors (T, M) are not given for type II (gray) groups.

ID "a.b.c":

- a: number of the point group (PG column), 1-32 for the crystallographic point groups and 33-37 for the other
  settings.
- b: number in the point group, 1 for type I (the point group itself), 2 for type II (gray group, with the time
  reversal 1'), and 3 or more for type III (black-white groups).
- c: serial number of the magnetic point groups, 1-157.

Type: I (ordinary), II (gray), III (black-white, some operations combined with the time reversal).
""",
    "MSG": """# Magnetic space groups

All 1651 magnetic space groups: symmetry operations and Wyckoff sites. Only the standard setting of the BNS
description is given for each magnetic space group.

ID "a.b" is the BNS number:

- a: number of the space group.
- b: serial number in each crystal family, starting from 1 in triclinic, monoclinic, orthorhombic, tetragonal,
  trigonal and hexagonal, and cubic groups.

The BNS and OG columns are the Belov-Neronova-Smirnova and Opechowski-Guccione symbols. Type: I (ordinary),
II (gray), III (black-white, some operations combined with the time reversal), IV (black-white, the time reversal
combined with a translation).
""",
}


# ==================================================
def _links(gid, tag, names):
    cells = []
    for k in names:
        path = os.path.relpath(os.path.join(group_dir(gid, tag), k + ".pdf"), os.path.join(__top_dir__, "docs", "src"))
        if k in kinds["MPG"] and not has_pdf(k, gid):  # no magnetic tensor for type II.
            cells.append("-")
            continue
        assert os.path.isfile(os.path.join(__top_dir__, "docs", "src", path)), path
        cells.append(f"[PDF]({path})")
    return cells


# ==================================================
def create_group_pages():
    info, group = load_info(), load_group()
    for tp in ["PG", "SG", "MPG", "MSG"]:
        names = (["harmonics"] if tp == "PG" else []) + kinds[tp]
        ids = info["id_set"][tp]["all"] + (info["id_set"]["PG"]["complex"] if tp == "PG" else [])
        if tp in ["PG", "SG"]:
            cols = ["No.", "Schoenflies", "international", "tag"]
        elif tp == "MPG":
            cols = ["ID", "international", "PG", "type"]
        else:
            cols = ["ID", "BNS", "OG", "type"]
        cols += [label[k] for k in names]

        lines = [head[tp]]
        for crystal in crystal_name:
            lst = [i for i in ids if group[i]["info"].crystal == crystal]
            lines += [
                f"\n## {crystal.capitalize()}\n",
                "| " + " | ".join(cols) + " |",
                "|" + "|".join([":---"] * len(cols)) + "|",
            ]
            for gid in lst:
                g = group[gid]["info"]
                tag = info["tag"][gid]
                no = gid.split(":")[1]
                if tp in ["PG", "SG"]:
                    cells = [no, f"${g.schoenflies}$", f"${g.international}$", f"`{tag}`"]
                elif tp == "MPG":
                    cells = [no, f"${g.international}$", "$" + group[g.PG]["info"].schoenflies + "$", g.type]
                else:
                    cells = [no, f"${g.BNS}$", f"${g.OG}$", g.type]
                lines.append("| " + " | ".join(cells + _links(gid, tag, names)) + " |")
        write(f"{tp.lower()}_pdf.md", lines)


# ==================================================
def create_examples():
    lines = [
        "# Examples",
        "",
        "Sample code for each feature, in `docs/examples`. Run each file from anywhere after `pip install -e .`, e.g.",
        "`python docs/examples/group.py`.",
    ]
    readme = open(os.path.join(doc_dir, "examples", "README.md")).read()
    for name, desc in re.findall(r"\| \[(\w+\.py)\]\(\S+\) \| (.+) \|", readme):
        first = next(l for l in open(os.path.join(doc_dir, "examples", name)) if l.startswith(("from ", "import ")))
        lines += ["", f"## {name}", "", f"{desc[0].upper() + desc[1:]}.", ""]
        lines += [f"```{{literalinclude}} ../examples/{name}", ":language: python", f":start-at: {first.strip()}", "```"]
    write("examples.md", lines)


# ==================================================
def write(name, lines):
    out = os.path.join(doc_dir, "src", name)
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"created {out}")


# ==================================================
if __name__ == "__main__":
    create_group_pages()
    create_examples()
