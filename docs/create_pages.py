# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create generated pages of the documentation.
    - docs/src/harmonics_pdf.md: list of harmonics PDFs of all 47 point groups (multipiegen/pdf/PG).
    - docs/src/examples.md: code of the examples listed in docs/examples/README.md.

    python docs/create_pages.py    (run by "make build" in docs)
"""

import os
import re

from multipiegen import __top_dir__
from multipiegen.data.data_group_pg import pg_id_set, pg_info

head = """# Harmonics of point groups

Harmonics of all 47 point groups (32 crystallographic point groups, 5 other settings "-1", and 10 complex groups "-c").
Each PDF lists the harmonics up to rank 11, with the polar (Q) and axial (G) symbols side by side.
The expression is that of the polar harmonics.

The PDFs are created by `python -m multipiegen.create_pdf.create_pdf_harmonics` in `multipiegen/pdf/PG/{no}-{tag}/`.
"""

doc_dir = os.path.dirname(os.path.abspath(__file__))
crystal_name = ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"]


# ==================================================
def entry(pg_id):
    info = pg_info[pg_id]
    no = int(pg_id.split(":")[1])
    path = f"../../multipiegen/pdf/PG/{no:03d}-{info.tag}/harmonics.pdf"
    assert os.path.isfile(os.path.join(__top_dir__, path[6:])), path
    return f"| {no} | ${info.schoenflies}$ | ${info.international}$ | `{info.tag}` | [PDF]({path}) |"


# ==================================================
def create_pdf_list():
    lines = [head]
    for crystal in crystal_name:
        lines.append(f"\n## {crystal.capitalize()}\n")
        lines.append("| No. | Schoenflies | international | tag | harmonics |")
        lines.append("|---:|:---|:---|:---|:---|")
        lines += [entry(i) for i in pg_id_set["crystal"][crystal]]

    write("harmonics_pdf.md", lines)


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
    create_pdf_list()
    create_examples()
