# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create LaTeX and PDF of each group (based on MultiPie others/tools/create_PDF/create_pdf_character.py,
create_pdf_symmetry_operation.py, create_pdf_wyckoff.py), from group.pkl.

    - multipiegen/pdf/PG/{no}-{tag}/{character_table,symmetry_operation,wyckoff_site,wyckoff_bond,
      representation_matrix}.{tex,pdf}
    - multipiegen/pdf/SG/{no}-{tag}/{symmetry_operation,wyckoff_site,wyckoff_bond}.{tex,pdf} ("^" in tag -> "_")
    - multipiegen/pdf/MPG/{id}/{symmetry_operation,wyckoff_site,response_tensor_X}.{tex,pdf}
      (X = Q, G, and T, M except for type II, from response_tensor.pkl)
    - multipiegen/pdf/MSG/{id}/{symmetry_operation,wyckoff_site}.{tex,pdf}

Differences from MultiPie:

    - complex PG ("-c") also has wyckoff_site and wyckoff_bond (same as the real PG).
    - the caption of a Wyckoff bond shows the bond symmetry.

Run:  python -m multipiegen.create_pdf.create_pdf_group [--group PG,SG,MPG,MSG] [--kind ...] [--serial]
"""

import argparse
import os
import time
from multiprocessing import Pool

import numpy as np
import sympy as sp
from multipie.util.util import to_latex
from multipie.util.util_pdf_latex import PDFviaLaTeX
from multipie.util.util_tag import TagIrrep, TagSymmetryOperation

from multipiegen import __top_dir__
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info
from multipiegen.util.util_representation_matrix import representation_matrix
from multipiegen.util.util_response_tensor import load_response_tensor

pdfdir = os.path.join(__top_dir__, "multipiegen", "pdf")
kinds = {
    "PG": ["character_table", "symmetry_operation", "wyckoff_site", "wyckoff_bond", "representation_matrix"],
    "SG": ["symmetry_operation", "wyckoff_site", "wyckoff_bond"],
    "MPG": [
        "symmetry_operation",
        "wyckoff_site",
        "response_tensor_Q",
        "response_tensor_G",
        "response_tensor_T",
        "response_tensor_M",
    ],
    "MSG": ["symmetry_operation", "wyckoff_site"],
}

_shared = {}


# ==================================================
def _data():
    """info and group (loaded once per process)."""
    if "info" not in _shared:
        _shared["info"] = load_info()
        _shared["group"] = load_group()
    return _shared["info"], _shared["group"]


def _response_tensor():
    """response tensors (loaded once per process)."""
    if "response_tensor" not in _shared:
        _shared["response_tensor"] = load_response_tensor()
    return _shared["response_tensor"]


# ==================================================
def has_pdf(kind, gid):
    """Is the PDF of kind created for the group ? (magnetic tensors T, M are not for type II MPG, as MultiPie)."""
    if kind in ["response_tensor_T", "response_tensor_M"]:
        return _data()[1][gid]["info"].type != "II"
    return True


# ==================================================
def group_dir(gid, tag):
    """
    Output directory of a group.

    Args:
        gid (str): group id, e.g. "PG:32".
        tag (str): group tag.

    Returns:
        - (str) -- directory, e.g. multipiegen/pdf/PG/032-Oh, multipiegen/pdf/SG/221-Oh_1, multipiegen/pdf/MSG/1.1.
    """
    tp, no = gid.split(":")
    if tp in ["PG", "SG"]:
        return os.path.join(pdfdir, tp, f"{int(no):03d}-{tag.replace('^', '_')}")
    return os.path.join(pdfdir, tp, no)


# ==================================================
def group_latex(gid, g_info):
    """Title of a group in LaTeX (as MultiPie Group.latex(True))."""
    tp, ID = gid.split(":")
    sp_ = r"\quad"
    if tp in ["PG", "SG"]:
        r = f"{tp} No. {ID}{sp_}${g_info.schoenflies}${sp_}${g_info.international}$"
    elif tp == "MPG":
        r = f"{tp} No. {ID}{sp_}${g_info.international}$"
    else:
        r = f"{tp} No. {ID}{sp_}${g_info.BNS}$"
    if g_info.setting != "":
        r += f"{sp_}({g_info.setting} setting)"
    if tp in ["PG", "SG"]:
        r += f"{sp_}[ {g_info.crystal} ]"
    else:
        r += f"{sp_}[ Type {g_info.type}, {g_info.crystal} ]"
    return r


# ==================================================
def _pdf(name, gid, tag, landscape=False, suffix=""):
    _, group = _data()
    h_dir = group_dir(gid, tag)
    os.makedirs(h_dir, exist_ok=True)
    pdf = PDFviaLaTeX(name, dir=h_dir, pt=8, style="narrowest", landscape=landscape, english=True)
    pdf.title(group_latex(gid, group[gid]["info"]) + suffix)
    return pdf


# ==================================================
def _plus_set(so):
    return r"* plus set:\,\,\," + r",\quad ".join(["$+" + to_latex(i, "vector") + "$" for i in so["plus_set"]])


# ==================================================
def create_character_table(gid, tag):
    _, group = _data()
    character = group[gid]["character"]
    irreps = list(character["table"].keys())
    ir_latex = [TagIrrep.latex(i) for i in irreps]
    pdf = _pdf("character_table", gid, tag, landscape=True)

    col = [TagSymmetryOperation.latex(i[0]) + f"({len(i)})" for i in character["conjugacy"]]
    wp, wm = sp.symbols(r"\omega \omega^*")
    sub = {sp.Symbol("wp"): wp, sp.Symbol("wm"): wm}
    has_w = any(c.has(sp.Symbol("wp"), sp.Symbol("wm")) for v in character["table"].values() for c in v)
    tbl = [[sp.latex(sp.S(c).subs(sub)) for c in v] for v in character["table"].values()]

    tt = "\n* character table"
    if has_w:
        tt += r" ($\omega=e^{2\pi i/3}$)"
    pdf.text(tt)
    pdf.text(r"\begin{quote}")
    rc = group[gid]["info"].schoenflies
    pdf.table(tbl, ir_latex, col, rc, rcmath=True, rmath=True, cmath=True, tmath=True, cpos="c" + "r" * len(col))
    pdf.text(r"\end{quote}")

    pdf.text("\n" + r"* polar $\leftrightarrow$ axial conversion")
    pdf.text(r"\begin{quote}")
    tbl = [[TagIrrep.latex(i) + r"\,\,(" + TagIrrep.latex(j) + ")" for i, j in character["polar_axial_conversion"].items()]]
    pdf.simple_table(tbl, tmath=True)
    pdf.text(r"\end{quote}")

    def decomp(lst):
        return sp.latex(sum([c * sp.Symbol(TagIrrep.latex(e)) for c, e in lst]))

    tbl = [
        [decomp(character["symmetric_product"][(ir1, ir2)]) if j >= i else "" for j, ir2 in enumerate(irreps)]
        for i, ir1 in enumerate(irreps)
    ]
    pdf.text("\n* symmetric product")
    pdf.text(r"\begin{quote}")
    pdf.table(tbl, ir_latex, ir_latex, rmath=True, cmath=True, tmath=True, cpos="c|" + "c" * len(irreps))
    pdf.text(r"\end{quote}")

    tbl = [decomp(character["anti_symmetric_product"][ir]) for ir in irreps]
    tbl = ["-" if v == "0" else v for v in tbl]
    pdf.text("\n* anti-symmetric product")
    pdf.text(r"\begin{quote}")
    pdf.table([tbl], [""], ir_latex, rmath=True, cmath=True, tmath=True)
    pdf.text(r"\end{quote}")

    pdf.build()


# ==================================================
def create_symmetry_operation(gid, tag):
    _, group = _data()
    tp = gid.split(":")[0]
    so = group[gid]["symmetry_operation"]
    pdf = _pdf("symmetry_operation", gid, tag)

    if tp in ["PG", "SG"]:
        gen = [TagSymmetryOperation.latex(i) for i in so["generator"]]
        pdf.text("\n* generator : $" + r",\,\,".join(gen) + "$\n")

    if tp == "PG":
        cc = [[TagSymmetryOperation.latex(i) for i in j] for j in group[gid]["character"]["conjugacy"]]
        pdf.text("* conjugacy class")
        pdf.text(r"\begin{quote}")
        pdf.simple_table([["[" + j[0] + r"]:", r",\,\,\,".join(j)] for j in cc], tmath=True, cpos="rl")
        pdf.text(r"\end{quote}" + "\n")

    ss = r"* symmetry operation"
    if tp == "SG":
        ss += r"\quad" + r",\quad ".join(["$+" + to_latex(i, "vector") + "$" for i in so["plus_set"]])
    pdf.text(ss)

    ops = [TagSymmetryOperation.latex(i) for i in so["tag"]]
    if tp in ["PG", "MPG"]:
        mat = [to_latex(i, "matrix") for i in so["fractional"]]
    else:
        mat = [to_latex(i[0:3, :], "matrix") for i in so["fractional"]]
    det = [str(i) for i in so["det"]]

    row = [str(i + 1) for i in range(len(ops))]
    if tp in ["PG", "SG"]:
        tbl = [list(i) for i in zip(ops, mat, det)]
        col, cpos = ["tag", "matrix (polar)", "det"], "cccc"
    else:
        tbl = [list(i) for i in zip(ops, mat, det, [str(i) for i in so["tr_sign"]])]
        col, cpos = ["tag", "matrix (polar)", "det", "TR"], "ccccc"

    cap = "Symmetry operations for 3d polar vector."
    pdf.table(tbl, row, col, "No.", rmath=True, tmath=True, stretch=1.3, caption=cap, cpos=cpos, long=True, hl=True)
    pdf.build()


# ==================================================
def create_representation_matrix(gid, tag):
    """Representation matrices from the harmonics (util_representation_matrix), as MultiPie."""
    rm = representation_matrix(tag)
    pdf = _pdf("representation_matrix", gid, tag, landscape=True)

    row, tbl = [], []
    for irrep, mat in rm["matrix"].items():
        row.append(TagIrrep.latex(irrep))
        tbl.append([TagSymmetryOperation.latex(t) + ":" + to_latex(m, "matrix") for t, m in zip(rm["tag"], mat)])

    col = [""] * 6
    cap = "Representation matrices"
    pdf.table(tbl, row, col, "Irrep.", hl=True, stretch=1.2, rmath=True, tmath=True, long=True, caption=cap)
    pdf.build()


# ==================================================
def _chunk(lst, n=24):
    return [lst[i : i + n] for i in range(0, len(lst), n)] if len(lst) > n else [lst]


# ==================================================
def _mapping_rows(no, first, mapping, n_empty):
    """
    Rows of a table with a long mapping split into several lines.

    Args:
        no (int): number (from 1).
        first (list): LaTeX cells of the first line (before mapping).
        mapping (list): SO numbers.
        n_empty (int): number of cells before mapping.

    Returns:
        - (list) -- row headers.
        - (list) -- table rows.
    """
    mlst = _chunk(mapping)
    if len(mlst) == 1:
        return [no], [first + [r"{\tt " + str(mlst[0]).replace(" ", "") + "}"]]
    row = [no] + [""] * (len(mlst) - 1)
    tbl = [first + [r"{\tt " + str(mlst[0])[:-1].replace(" ", "") + ",}"]]
    for mi in mlst[1:-1]:
        tbl.append([""] * n_empty + [r"{\tt " + str(mi)[1:-1].replace(" ", "") + ",}"])
    tbl.append([""] * n_empty + [r"{\tt " + str(mlst[-1])[1:].replace(" ", "") + "}"])
    return row, tbl


# ==================================================
def create_wyckoff_site(gid, tag):
    _, group = _data()
    tp = gid.split(":")[0]
    g = group[gid]
    pdf = _pdf("wyckoff_site", gid, tag)

    if tp == "SG":
        pdf.text(_plus_set(g["symmetry_operation"]))

    for w, val in g["wyckoff"]["site"].items():
        row, tbl = [], []
        for no, (pos, m) in enumerate(zip(val["conventional"], val["mapping"])):
            r, t = _mapping_rows(no + 1, ["$" + to_latex(pos, "vector") + "$"], m, 1)
            row += r
            tbl += t
        cap = r"Wyckoff site: {\tt " + str(w) + "}" + r", site symmetry: {\tt " + val["symmetry"] + "}"
        pdf.table(tbl, row, ["position", "mapping"], "No.", hl=True, stretch=1.2, long=True, caption=cap)

    pdf.build()


# ==================================================
def create_wyckoff_bond(gid, tag):
    _, group = _data()
    tp = gid.split(":")[0]
    g = group[gid]
    wp = g["wyckoff"]
    pdf = _pdf("wyckoff_bond", gid, tag)

    if tp == "SG":
        pdf.text("\n" + _plus_set(g["symmetry_operation"]))
    pdf.text("\n")

    for s_wp, val in wp["site"].items():
        if s_wp == "1o":  # bonds at the origin of point group (as MultiPie).
            continue
        pdf.text(r"* Wyckoff site: {\tt " + str(s_wp) + "}" + r", site symmetry: {\tt " + val["symmetry"] + "}")
        for b_wp in val["bond"]:
            bond = wp["bond"][b_wp]
            row, tbl = [], []
            for no, (b, m) in enumerate(zip(bond["conventional"], bond["mapping"])):
                v, c = b[0:3], b[3:6]
                first = ["$" + to_latex(v, "vector") + "$", "$" + to_latex(c, "vector") + "$"]
                r, t = _mapping_rows(no + 1, first, m, 2)
                row += r
                tbl += t
            cap = r"Wyckoff bond: {\tt " + str(b_wp) + "}" + r", bond symmetry: {\tt " + bond["symmetry"] + "}"
            pdf.table(tbl, row, ["vector", "center", "mapping"], "No.", hl=True, stretch=1.2, long=True, caption=cap)

    pdf.build()


def create_response_tensor(gid, tag, X):
    """Response tensors of type X up to rank 4 (response_tensor.pkl), as MultiPie."""
    d = _response_tensor()[gid]
    pdf = _pdf(f"response_tensor_{X}", gid, tag, landscape=True, suffix=f" [{X} tensor]")

    empty = True
    for (Xt, rank, opt), (m, ex) in d.items():
        if Xt != X:
            continue
        pdf.text(f"* Rank {rank} tensor." if opt == "" else f"* Rank {rank} tensor ({opt}).")
        if not np.all(m == sp.S(0)):
            empty = False
            pdf.equation(to_latex(m, "matrix"))
            pdf.equation([sp.latex(sp.Eq(i, j)) for i, j in ex.items()])
    if empty:
        pdf.text()
    pdf.build()


_creator = {
    "character_table": create_character_table,
    "symmetry_operation": create_symmetry_operation,
    "wyckoff_site": create_wyckoff_site,
    "wyckoff_bond": create_wyckoff_bond,
    "representation_matrix": create_representation_matrix,
    "response_tensor_Q": lambda gid, tag: create_response_tensor(gid, tag, "Q"),
    "response_tensor_G": lambda gid, tag: create_response_tensor(gid, tag, "G"),
    "response_tensor_T": lambda gid, tag: create_response_tensor(gid, tag, "T"),
    "response_tensor_M": lambda gid, tag: create_response_tensor(gid, tag, "M"),
}


def _worker(arg):
    kind, gid, tag = arg
    _creator[kind](gid, tag)
    return arg


# ==================================================
def create_pdf_group(group_types=("PG", "SG", "MPG", "MSG"), kind=None, serial=False):
    """
    Create PDFs of each group.

    Args:
        group_types (tuple, optional): group types, subset of ("PG", "SG", "MPG", "MSG").
        kind (list, optional): kinds of PDF (all kinds of each group type if None).
        serial (bool, optional): run in a single process ?
    """
    info, _ = _data()
    tasks = []
    for tp in group_types:
        ids = info["id_set"][tp]["all"] + (info["id_set"]["PG"]["complex"] if tp == "PG" else [])
        for k in kinds[tp]:
            if kind is None or k in kind:
                tasks += [(k, gid, info["tag"][gid]) for gid in ids if has_pdf(k, gid)]

    t0 = time.time()
    if serial:
        for t in tasks:
            _worker(t)
    else:
        with Pool(os.cpu_count()) as pool:
            for _ in pool.imap_unordered(_worker, tasks, chunksize=2):
                pass
    print(f"created {len(tasks)} PDFs in {time.time() - t0:.0f} s.")


# ==================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create PDFs of each group")
    parser.add_argument("--group", default="PG,SG,MPG,MSG", help="group types, e.g. PG,SG")
    parser.add_argument("--kind", default=None, help="kinds, e.g. character_table,wyckoff_site")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    args = parser.parse_args()
    create_pdf_group(args.group.split(","), None if args.kind is None else args.kind.split(","), args.serial)
