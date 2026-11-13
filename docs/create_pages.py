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
    - docs/src/data.md: list of the data tables (multipiegen/data) with their variables, the scripts in others/creation
      that create them, and the modules that use them.
    - docs/src/binary.md: list of the binary data (multipiegen/binary) with their top-level keys, and the description
      of the dictionary keys (the header of each binary).

    python docs/create_pages.py    (run by "make build" in docs)
"""

import ast
import importlib
import os
import re

from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__, __top_dir__
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
    """
    The page of the examples, as a notebook (MyST-NB) executed when the docs are built: the code of each example after
    its docstring, split into cells at the "# ===..." lines, with the output below each cell.
    """
    lines = [
        "---",
        "jupytext:",
        "  text_representation:",
        "    format_name: myst",
        "kernelspec:",
        "  display_name: Python 3",
        "  language: python",
        "  name: python3",
        "---",
        "",
        "# Examples",
        "",
        "Sample code for each feature, in `docs/examples`. Run each file from anywhere after `pip install -e .`, e.g.",
        "`python docs/examples/group.py`. The outputs below are made when the documentation is built.",
    ]
    readme = open(os.path.join(doc_dir, "examples", "README.md")).read()
    for name, desc in re.findall(r"\| \[(\w+\.py)\]\(\S+\) \| (.+) \|", readme):
        src = open(os.path.join(doc_dir, "examples", name)).read()
        body = ast.parse(src).body
        start = body[1].lineno if isinstance(body[0], ast.Expr) else body[0].lineno  # after the module docstring.
        code = src.split("\n")
        first = start - 1
        cells = [[]]
        for l in code[first:]:
            if l.startswith("# ===="):
                cells.append([])
            else:
                cells[-1].append(l)
        lines += ["", f"## {name}", "", f"{desc[0].upper() + desc[1:]}."]
        for cell in cells:
            cell = "\n".join(cell).strip("\n")
            if cell:
                lines += ["", "```{code-cell} ipython3", cell, "```"]
    write("examples.md", lines)


# ==================================================
data_content = {
    "data_basis.py": "orbital and basis definitions: atomic orbitals, spin, LM/JM/LMs/Lgs basis, bra-ket blocks, "
    "LaTeX names, Wannier90 aliases, cubic/hexagonal orbitals (internal basis of multipolar harmonics), harmonics keys, "
    "max ranks and rank names",
    "data_character_table.py": "character tables of all 47 point groups",
    "data_compatibility_relation.py": "compatibility relations of the Oh and D6h subgroups",
    "data_group_pg.py": "group information of the point groups",
    "data_group_sg.py": "group information of the space groups",
    "data_group_mpg.py": "group information of the magnetic point groups",
    "data_group_msg.py": "group information of the magnetic space groups",
    "data_parent_pg.py": "parent point groups of each point group",
    "data_parent_sg.py": "parent space groups of each space group",
    "data_response_tensor.py": "multipoles in the response tensors and their cartesian components",
    "data_root_cluster.py": "Oh and D6h clusters with the representative site and vector",
    "data_top_harmonics.py": "harmonics of the top groups (Oh, D6h)",
    "data_wyckoff_site.py": "Wyckoff sites of the point and space groups",
    "data_wyckoff_bond_sym.py": "Wyckoff bonds of the point and space groups with the bond symmetry",
    "data_wyckoff_mpg.py": "Wyckoff sites of the magnetic point groups",
    "data_wyckoff_msg.py": "Wyckoff sites of the magnetic space groups",
}
data_from_multipie = [
    "data_character_table.py",
    "data_compatibility_relation.py",
    "data_group_mpg.py",
    "data_group_msg.py",
    "data_group_pg.py",
    "data_group_sg.py",
    "data_top_harmonics.py",
    "data_wyckoff_site.py",
    "data_wyckoff_bond_sym.py",
    "data_wyckoff_mpg.py",
    "data_wyckoff_msg.py",
    "data_response_tensor.py",
]
data_var_note = {
    ("data_basis.py", "atomic_lgs_basis"): "regular atomic spinful basis order of the orbitals (Lgs)",
    ("data_basis.py", "cubic_basis"): "spinful basis order of cubic_orbital_basis",
    ("data_basis.py", "hexagonal_basis"): "spinful basis order of hexagonal_orbital_basis",
    ("data_character_table.py", "character"): "{ PG: ( [[SO in conjugacy class]], { irrep: [character] } ) }",
    ("data_compatibility_relation.py", "compatibility_relation"): "{ PG: (top PG, { (top irrep, comp): ((polar irrep, "
    "comp, coeff), (axial irrep, comp, coeff)) }) }",
    ("data_top_harmonics.py", "top_harmonics_tag"): "tags of the top harmonics",
    ("data_top_harmonics.py", "top_harmonics"): "top harmonics (cartesian expression) of Oh and D6h",
    ("data_wyckoff_bond_sym.py", "wyckoff_bond_data"): "{ group: { Wyckoff: (vector, center, mapping, bond_symmetry) } }",
    ("data_wyckoff_mpg.py", "wyckoff_mpg_data"): "{ MPG_id: { Wyckoff: (position, site-symmetry) } }",
    ("data_wyckoff_msg.py", "wyckoff_msg_data"): "{ MSG_id: { Wyckoff: (position, site-symmetry) } }",
    ("data_wyckoff_site.py", "wyckoff_site_data"): "{ PG/SG: { Wyckoff: (position, site-symmetry) } } (point groups first)",
}

binary_content = [
    (
        "harmonics.pkl",
        "expression and u-matrix of the harmonics of all PG (rank <= 11), looked up by a multipole key",
        "create_harmonics",
        "util_harmonics.load_harmonics",
    ),
    (
        "info.pkl",
        "overall information of all groups, harmonics, basis, root cluster, response tensor",
        "create_info",
        "util_info.load_info",
    ),
    (
        "group.pkl",
        "data of each PG, SG, MPG, MSG: symmetry operations, character, Wyckoff sites and bonds",
        "create_group",
        "util_group.load_group",
    ),
    (
        "response_tensor.pkl",
        "response tensors (rank 0-4) of all MPG",
        "create_response_tensor",
        "util_response_tensor.load_response_tensor",
    ),
    (
        "multipolar_harmonics.pkl",
        "multipolar harmonics (s=1-3, l<=3) of all PG",
        "create_multipolar_harmonics",
        "util_harmonics_multipole.load_multipolar_harmonics",
    ),
    (
        "atomic_multipole.pkl",
        "matrices of atomic multipoles in the LMs/JML/LM basis",
        "create_atomic_multipole",
        "util_atomic_multipole.load_atomic_multipole",
    ),
    (
        "atomic_multipole_group.pkl",
        "atomic multipole SAMB of all PG in the JML/Lgs/Lg basis",
        "create_atomic_multipole_group",
        "util_atomic_multipole.load_atomic_multipole_group",
    ),
    (
        "cluster_scalar_basis.pkl",
        "orthonormal (scalar) cluster basis of the Oh and D6h clusters",
        "create_cluster_scalar_basis",
        "util_cluster_basis.load_cluster_basis",
    ),
    (
        "cluster_vector_basis.pkl",
        "orthonormal vector cluster basis (s=1) of the Oh and D6h clusters",
        "create_cluster_vector_basis",
        "util_cluster_vector_basis.load_cluster_vector_basis",
    ),
    (
        "cluster_bond_basis.pkl",
        "orthonormal bond cluster basis of the Oh and D6h clusters",
        "create_cluster_bond_basis",
        "util_bond_cluster_basis.load_bond_cluster_basis",
    ),
    (
        "cluster_samb.pkl",
        "cluster SAMB of each Wyckoff site and bond of the real PG and SG",
        "create_cluster_samb",
        "util_cluster_samb.load_cluster_samb",
    ),
]


# ==================================================
def _var_comment(src_lines, lineno):
    """Comment lines just above a top-level assignment (without separators)."""
    lst = []
    i = lineno - 2
    while i >= 0 and src_lines[i].startswith("#"):
        line = src_lines[i].lstrip("#").strip()
        if line and not set(line) <= set("="):
            lst.insert(0, line)
        i -= 1
    return " ".join(lst)


def _doc_note(doc, name):
    """Description of a variable in the module docstring ("- name ...: description")."""
    lines = (doc or "").splitlines()
    for i, line in enumerate(lines):
        m = re.match(r"^\s*-\s*" + re.escape(name) + r"\b[^:]*:?\s*(.*)$", line)
        if m:
            note = m.group(1).strip().rstrip(".")
            if not note and i + 1 < len(lines):  # "- name." followed by the structure.
                note = lines[i + 1].strip().rstrip(".")
            return note
    return ""


def _size(v):
    t = type(v).__name__
    return f"{t} ({len(v)})" if hasattr(v, "__len__") and not isinstance(v, str) else t


def _creators():
    """{data file: [scripts in others/creation that write it]}."""
    d = {}
    cdir = os.path.join(__top_dir__, "others", "creation")
    for f in sorted(os.listdir(cdir)):
        if f.startswith("create_") and f.endswith(".py"):
            for name in set(re.findall(r"os\.path\.join\(cur_dir, \"(data_[a-z_]+\.py)\"\)", open(os.path.join(cdir, f)).read())):
                d.setdefault(name, []).append(f"others/creation/{f}")
    return d


def _users():
    """{data module name: [modules of multipiegen that import it]}."""
    d = {}
    pdir = os.path.join(__top_dir__, "multipiegen")
    for sub in ["util", "create_binary", "create_pdf"]:
        for f in sorted(os.listdir(os.path.join(pdir, sub))):
            if f.endswith(".py"):
                for name in set(re.findall(r"multipiegen\.data\.(data_[a-z_]+)", open(os.path.join(pdir, sub, f)).read())):
                    d.setdefault(name, []).append(f"{sub}/{f[:-3]}")
    return d


# ==================================================
def create_data_page():
    creators, users = _creators(), _users()
    ddir = os.path.join(__top_dir__, "multipiegen", "data")
    files = sorted(f for f in os.listdir(ddir) if f.startswith("data_") and f.endswith(".py"))
    lines = [
        "# Data tables",
        "",
        "The data tables in `multipiegen/data` are the source of all generated data. They are curated by hand; some of "
        "them are first created by a script in `others/creation` (it writes next to itself, and the result is copied "
        'into `multipiegen/data`). "origin" is "MultiPie Ver.2.4.0" for the tables taken from MultiPie Ver.2.4.0 '
        '(some of them modified), and "new" otherwise. The tables taken from MultiPie Ver.2.4.0 were copied once and '
        "have been maintained here since then. The generation never reads the data of MultiPie Ver.2; it uses only its "
        "general functions and utilities.",
        "",
        "Data that can be derived from these tables are computed in `multipiegen` instead of being stored:",
        "",
        "| data | computed in | from |",
        "|---|---|---|",
        "| `*_tag`, `*_id`, `*_id_set`, `group_name_list` | `util/util_group_data` | `*_info` of `data_group_pg/sg/mpg/msg` |",
        "| `top_harmonics_symmetry`, `top_harmonics_u_matrix` | `util/util_top_harmonics` | `data_top_harmonics` "
        "(tesseral part from the multipolar harmonics) |",
        "| `data_harmonics_mapping` | `util/util_harmonics_mapping` | `data_compatibility_relation`, "
        "`data_top_harmonics.top_harmonics_tag` |",
        "| `op_mapping_d6h_oh_cluster`, `op_mapping_parent_sg` | `util/util_op_mapping` | `data_parent_pg`, `data_parent_sg` |",
        "| polar-axial conversion, symmetric and anti-symmetric products | `create_binary/create_group.character_plus` | "
        "`data_character_table` and the symmetry operations |",
        "",
        "| file | content | origin | created by | used by |",
        "|---|---|---|---|---|",
    ]
    for f in files:
        name = f[:-3]
        assert f in data_content, f"no description of {f} in docs/create_pages.py"
        created = ", ".join(f"`{c}`" for c in creators.get(f, [])) or "-"
        used = ", ".join(f"`{u}`" for u in users.get(name, [])) or "-"
        origin = "MultiPie Ver.2.4.0" if f in data_from_multipie else "new"
        lines.append(f"| [`{f}`](#{name.replace('_', '-')}-py) | {data_content[f]} | {origin} | {created} | {used} |")

    for f in files:
        name = f[:-3]
        src = open(os.path.join(ddir, f)).read()
        tree = ast.parse(src)
        doc = ast.get_docstring(tree)
        mod = importlib.import_module(f"multipiegen.data.{name}")
        lines += ["", f"## {f}", "", data_content[f][0].upper() + data_content[f][1:] + "."]
        created = creators.get(f, [])
        if created:
            lines += ["", "Created by " + ", ".join(f"`{c}`" for c in created) + "."]
        lines += ["", "| variable | type (size) | description |", "|---|---|---|"]
        seen = set()
        for node in tree.body:
            if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                v = node.targets[0].id
                if v in seen:
                    continue
                seen.add(v)
                note = (
                    data_var_note.get((f, v))
                    or _var_comment(src.splitlines(), node.lineno)
                    or _doc_note(doc, v)
                    or "see the description below"
                )
                lines.append(f"| `{v}` | {_size(getattr(mod, v))} | {note} |")
        if doc:
            lines += ["", "```text", doc, "```"]
    write("data.md", lines)


# ==================================================
def create_binary_page():
    lines = [
        "# Binary data",
        "",
        "The binary data in `multipiegen/binary` are created from the data tables by the modules in "
        "`multipiegen/create_binary`. They use the format of MultiPie Ver.2's `BinaryManager` (gzip-compressed pickle "
        'with a header), so they can be read without MultiPieGen, e.g. `BinaryManager("info", topdir=...)`. '
        "`str(data)` shows the header, the description of the dictionary keys given in each section below.",
        "",
        "| file | content | top-level keys | create | load |",
        "|---|---|---|---|---|",
    ]
    sections = []
    for f, content, create, load in binary_content:
        name = f[:-4]
        bm = BinaryManager(name, topdir=__bin_dir__)
        keys = list(bm.keys())
        if len(keys) <= 10:
            top = ", ".join(f"`{k}`" for k in keys)
        else:
            top = f"{len(keys)} keys: " + ", ".join(f"`{k}`" for k in keys[:2]) + f", ..., `{keys[-1]}`"
        anchor = name.replace("_", "-") + "-pkl"
        cmd = f"`python -m multipiegen.create_binary.{create}`"
        ld = f"`multipiegen.util.{load}()`"
        lines.append(f"| [`{f}`](#{anchor}) | {content} | {top} | {cmd} | {ld} |")
        size = os.path.getsize(os.path.join(__bin_dir__, f))
        sections += [
            "",
            f"## {f}",
            "",
            f"{content[0].upper() + content[1:]}.",
            "",
            f"- create: {cmd}",
            f"- load: {ld}",
            f"- top-level keys: {top}",
            f"- size: {size / 1024:,.0f} KB",
            "",
            "Description of the dictionary keys (header of the binary):",
            "",
            "```text",
            str(bm).strip().strip('"').strip(),
            "```",
        ]
    write("binary.md", lines + sections + binary_data_section())


# ==================================================
def binary_data_section():
    """Section of binary.md on the binary data for MultiPie (multipiegen/binary_data)."""
    from multipiegen.util.util_binary_data import binary_data_dir, group_file, group_types, load_binary_data

    lines = [
        "",
        "## Binary data for MultiPie (binary_data)",
        "",
        "`multipiegen/binary_data` holds the binary data to be built into MultiPie Ver.3 and later, in the layout of "
        "MultiPie's `multipie/binary_data`: `info.pkl` (the same as `binary/info.pkl`) and one file for each group, "
        "`{PG,SG,MPG,MSG}/{id}.pkl`. They are made from the binaries above, regrouped by group "
        "(the representation matrices of PG are added), by `python -m multipiegen.create_binary.create_binary_data`. "
        'Load a file with `BinaryManager("PG/PG_32", topdir=...)` (`:` of the id is replaced by `_` in the file names) or `multipiegen.util.util_binary_data.load_binary_data("PG:32")`. '
        "The SAMB and the multipolar harmonics hold no expression of the harmonics: look it up by `(X, l, Gamma, n)` of "
        'the key in `"harmonics"` of the PG file (T -> Q, M -> G); for SG and MPG, the PG file is that of `info.PG`.',
        "",
        "| group | files | keys | size |",
        "|---|---|---|---|",
    ]
    examples = {}
    for tp in group_types:
        d = os.path.join(binary_data_dir, tp)
        files = sorted(os.listdir(d))
        size = sum(os.path.getsize(os.path.join(d, f)) for f in files)
        gid = {"PG": "PG:32", "SG": "SG:221"}.get(tp, files[0][:-4].replace("_", ":", 1))  # file name -> id.
        bm = load_binary_data(gid)
        examples[tp] = (gid, bm)
        keys = ", ".join(f"`{k}`" for k in bm.keys())
        lines.append(f"| {tp} | {len(files)} (`{tp}/{files[0]}`, ...) | {keys} | {size / 1024:,.0f} KB |")
    lines += ["", "The keys of the complex PG (`-c`) lack `cluster_samb`.", ""]
    for tp, (gid, bm) in examples.items():
        lines += [
            f"### {tp} files",
            "",
            f"Description of the dictionary keys (header of `{group_file(gid)}.pkl`):",
            "",
            "```text",
            str(bm).strip().strip('"').strip(),
            "```",
            "",
        ]
    return lines


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
    create_data_page()
    create_binary_page()
