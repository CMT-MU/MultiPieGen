# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create LaTeX and PDF of overall information (based on MultiPie others/tools/create_PDF/create_pdf_group_info.py,
create_pdf_compatibility_relation.py, create_pdf_active_multipole.py).

    - multipiegen/pdf/info/{PG,SG,MPG,MSG}.{tex,pdf}: group information.
    - multipiegen/pdf/info/compatibility_relation.{tex,pdf}: compatibility relation of Oh and D6h subgroups, with
      the orbitals of each irrep (component), and the list of orbitals.
    - multipiegen/pdf/info/active_multipole.{tex,pdf}: active multipoles in magnetic point groups.

The data are taken from info.pkl and group.pkl.

Run:  python -m multipiegen.create_pdf.create_pdf_info
"""

import os
import random

import sympy as sp
from multipie.util.util import str_to_sympy, to_latex
from multipie.util.util_pdf_latex import PDFviaLaTeX
from multipie.util.util_tag import TagIrrep

from multipiegen import __top_dir__
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.data.data_top_harmonics_info import top_harmonics_symmetry
from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

pdfdir = os.path.join(__top_dir__, "multipiegen", "pdf", "info")
top_group = ["Oh", "D6h"]
component_name = "xyz"

_x, _y, _z = sp.symbols("x y z", real=True)


# ==================================================
def _changes(keys):
    """Row indices after which the key changes (horizontal lines of a table)."""
    return [i for i in range(len(keys) - 1) if keys[i] != keys[i + 1]]


# ==================================================
def create_group_info_pg(info, group):
    pdf = PDFviaLaTeX("PG", dir=pdfdir, pt=10, style="narrowest", english=True)

    ids = info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]
    col = ["tag", "", r"Sch\"onflies", "international", "setting", "crystal"]
    row, tbl = [], []
    for i in ids:
        g = group[i]["info"]
        row.append(r"\texttt{" + i.split(":")[1] + "}")
        tbl.append(
            [r"\texttt{" + g.tag + "}", r"\texttt{" + i + "}", f"${g.schoenflies}$", f"${g.international}$", g.setting, g.crystal]
        )
    hl = _changes([group[i]["info"].crystal for i in ids])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, "ID", caption="Point group information.", stretch=1.0, long=True, hl=hl, cpos="lllllll")
    pdf.text(r"\end{center}")
    pdf.build()


# ==================================================
def create_group_info_sg(info, group):
    pdf = PDFviaLaTeX("SG", dir=pdfdir, pt=10, style="narrowest", english=True)

    ids = info["id_set"]["SG"]["all"]
    col = ["tag", "", "", r"Sch\"onflies", "international", "setting", "crystal", "PG"]
    row, tbl = [], []
    for i in ids:
        g = group[i]["info"]
        no = i.split(":")[1]
        pg = "$" + group[g.PG]["info"].schoenflies + "$"
        row.append(r"\texttt{" + no + "}")
        tbl.append(
            [
                r"\texttt{" + g.tag + "}",
                r"\texttt{" + i + "}",
                r"\texttt{" + no + "}",
                f"${g.schoenflies}$",
                f"${g.international}$",
                g.setting,
                g.crystal,
                pg,
            ]
        )
    # lines between point groups (settings "-1" are not separated).
    hl = _changes([group[group[i]["info"].PG]["info"].tag.split("-")[0] for i in ids])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, "ID", caption="Space group information.", stretch=1.0, long=True, hl=hl, cpos="lllllllll")
    pdf.text(r"\end{center}")
    pdf.build()


# ==================================================
def create_group_info_mpg(info, group):
    pdf = PDFviaLaTeX("MPG", dir=pdfdir, pt=10, style="narrowest", english=True)

    ids = info["id_set"]["MPG"]["all"]
    col = ["tag", "international", "setting", "crystal", "PG", "type"]
    row, tbl = [], []
    for i in ids:
        g = group[i]["info"]
        pg = "$" + group[g.PG]["info"].schoenflies + "$"
        row.append(r"\texttt{" + i.split(":")[1] + "}")
        tbl.append([r"\texttt{MPG:" + g.tag + "}", f"${g.international}$", f"${g.setting}$", g.crystal, pg, g.type])
    hl = _changes([group[i]["info"].PG for i in ids])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, "ID", caption="Magnetic point group information.", stretch=1.0, long=True, hl=hl, cpos="lllllll")
    pdf.text(r"\end{center}")
    pdf.build()


# ==================================================
def create_group_info_msg(info, group):
    pdf = PDFviaLaTeX("MSG", dir=pdfdir, pt=10, style="narrowest", english=True)

    ids = info["id_set"]["MSG"]["all"]
    col = ["tag", "BNS", "OG", "crystal", "PG", "SG", "MPG", "type"]
    row, tbl = [], []
    for i in ids:
        g = group[i]["info"]
        pg = "$" + group[g.PG]["info"].schoenflies + "$"
        sg = "$" + group[g.SG]["info"].schoenflies + "$"
        mpg = r"\texttt{" + g.MPG.split(":")[1] + "} $(" + group[g.MPG]["info"].international + ")$"
        row.append(r"\texttt{" + i.split(":")[1] + "}")
        tbl.append([r"\texttt{MSG:" + g.tag + "}", f"${g.BNS}$", f"${g.OG}$", g.crystal, pg, sg, mpg, g.type])
    hl = _changes([group[i]["info"].SG for i in ids])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, "ID", caption="Magnetic space group information.", stretch=1.0, long=True, hl=hl, cpos="lllllllll")
    pdf.text(r"\end{center}")
    pdf.build()


# ==================================================
def _expression(e):
    """sympy expression of a top harmonics (string) in real x, y, z."""
    e = str_to_sympy(e) if isinstance(e, str) else e
    return sp.expand(e.subs({s: {"x": _x, "y": _y, "z": _z}[s.name] for s in e.free_symbols}))


# ==================================================
def orbital_top_key(info):
    """
    Top harmonics of each orbital.

    Args:
        info (BinaryManager): info.

    Returns:
        - (dict) -- {top group: {orbital: (top_key, c)}}, top harmonics = c * orbital, for orbitals in top harmonics.
    """
    bf = info["harmonics"]["basis_function"]
    names = list(info["harmonics"]["latex"]["orbital"].keys())
    rnd = random.Random(1)
    pts = [
        {_x: sp.Rational(rnd.randint(1, 97), 7), _y: sp.Rational(rnd.randint(1, 97), 11), _z: sp.Rational(rnd.randint(1, 97), 13)}
        for _ in range(3)
    ]

    dic = {}
    for top in top_group:
        expr = {k: _expression(e) for k, e in top_harmonics_symmetry[top].items()}
        dic[top] = {}
        for name in names:
            b = sp.expand(bf[name][0])
            rank = sp.Poly(b, _x, _y, _z).total_degree()
            bv = [b.subs(p) for p in pts]
            for k, e in expr.items():
                if k[0] != rank:
                    continue
                ratio = [sp.nsimplify(e.subs(p) / v) for p, v in zip(pts, bv)]
                if len(set(ratio)) == 1 and sp.expand(e - ratio[0] * b) == 0:
                    dic[top][name] = (k, ratio[0])
                    break
    return dic


# ==================================================
def _signed(c, s):
    """LaTeX math of s with the sign of c."""
    return "$" + ("-" if c < 0 else "") + s + "$"


# ==================================================
def _irrep_latex(s):
    """LaTeX of irrep with component, e.g. 'T1u(x)' or '-Eg(y)'."""
    sign = "-" if s.startswith("-") else ""
    s = s.lstrip("-")
    if s.count("("):
        ir, p = s.split("(")
        return "$" + sign + TagIrrep.latex(ir) + f"({p}" + "$"
    return "$" + sign + TagIrrep.latex(s) + "$"


# ==================================================
def _tag_latex(s):
    """LaTeX of point group tag, e.g. 'D4h' or 'D2d-1'."""
    if s.count("-"):
        t, p = s.split("-")
        return "$" + t[0] + r"_{\rm " + t[1:] + "}$(" + p + ")"
    return "$" + s[0] + r"_{\rm " + s[1:] + "}$"


# ==================================================
def _column_key(s):
    """(irrep, component) of compatibility-relation column, e.g. 'T1u(x)' -> ('T1u', 0), 'A1g' -> ('A1g', -1)."""
    if s.count("("):
        ir, p = s.split("(")
        return ir, component_name.index(p[0])
    return s, -1


# ==================================================
def _orbitals(X, top, irrep, comp, top_key, latex):
    """'/'-separated LaTeX of orbitals of X (Q or G) harmonics of (irrep, comp) in top group."""
    lst = []
    for key, (c, tk) in data_harmonics_mapping[top][X].items():
        if key[1] != irrep or key[3] != comp:
            continue
        for name, (k, r) in top_key.items():
            if k == tk:  # X(key) = c * top(tk) = c * r * orbital.
                lst.append((key[0], _signed(c * r, X + latex[name])))
    return r"{\footnotesize " + "/".join(s for _, s in sorted(lst, key=lambda x: x[0])) + "}"


# ==================================================
def create_compatibility_relation(info, group, top_key):
    cr = info["character"]["compatibility_relation"]
    latex = info["harmonics"]["latex"]["orbital"]
    pdf = PDFviaLaTeX("compatibility_relation", dir=pdfdir, pt=9, landscape=True, style="narrowest", english=True)
    pdf.title("Compatibility Relation for Point Groups")

    for top in top_group:
        cols = [i[0] for i in cr[(top, top)]]
        col = [_irrep_latex(i) for i in cols]
        row = ["Q", "G"]
        tbl = [[_orbitals(X, top, *_column_key(c), top_key[top], latex) for c in cols] for X in ["Q", "G"]]
        crystal = []
        for k, v in cr.items():
            if k[0] == top:
                row.append(_tag_latex(k[1]))
                tbl.append([_irrep_latex(i[1]) for i in v])
                crystal.append(group[info["id"][k[1]]]["info"].crystal)
        hl = [1] + [i + 2 for i in _changes(crystal)]

        # vertical line between g and u.
        parity = [_column_key(c)[0][-1] for c in cols]
        cpos = "c|" + "".join("c" + ("|" if i in _changes(parity) else "") for i in range(len(cols)))
        cap = _tag_latex(top) + " subgroup."
        pdf.table(tbl, row, col, _tag_latex(top), caption=cap, stretch=1.0, cmath=False, tmath=False, hl=hl, cpos=cpos)
        pdf.text(r"\newpage")

    bf = info["harmonics"]["basis_function"]
    names = sorted(latex.keys(), key=lambda n: sp.Poly(bf[n][0], _x, _y, _z).total_degree())
    orbital_table(pdf, info, names, "Orbitals (harmonics).")
    pdf.build()


# ==================================================
def orbital_table(pdf, info, names, caption):
    """Table of orbitals: rank (first row of each rank only), tag, symbol and expression (lines between ranks)."""
    bf = info["harmonics"]["basis_function"]
    latex = info["harmonics"]["latex"]["orbital"]
    rank = [sp.Poly(bf[n][0], _x, _y, _z).total_degree() for n in names]
    row = [str(r) if i == 0 or r != rank[i - 1] else "" for i, r in enumerate(rank)]
    tbl = [[r"\texttt{" + n + "}", r"$\mathbb{X}" + latex[n] + "$", "$" + to_latex(bf[n][0]) + "$"] for n in names]
    col = ["tag", "symbol", "expression"]
    pdf.table(tbl, row, col, "rank", cpos="cllll", long=True, stretch=1.6, caption=caption, hl=_changes(rank))


# ==================================================
def _active_list(active, order):
    return [r"$\bullet$" if i in active else "" for i in order]


# ==================================================
def create_active_multipole(info, group, top_key):
    basis = info["response_tensor"]["multipole"]
    pdf = PDFviaLaTeX("active_multipole", dir=pdfdir, landscape=True, style="narrow", pt=9, english=True)
    pdf.title("Active Multipoles in Magnetic Point Groups")

    names = list(dict.fromkeys(basis["cubic"] + basis["hexagonal"]))
    bf = info["harmonics"]["basis_function"]
    names.sort(key=lambda n: sp.Poly(bf[n][0], _x, _y, _z).total_degree())
    orbital_table(pdf, info, names, "Multipoles in response tensor.")

    for crystal, top, cap in [("cubic", "Oh", "cubic subgroups"), ("hexagonal", "D6h", "hexagonal subgroups")]:
        order = basis[crystal]
        latex = info["harmonics"]["latex"]["orbital"]
        lbl = ["MPG", "Type", "Crystal", "X"] + [r"$\mathbb{X}" + latex[i] + "$" for i in order]
        # vertical lines between (rank, irrep, multiplicity) of the top group.
        tk = [top_key[top][n][0][0:3] for n in order]
        cpos = "lllll|" + "".join("l" + ("|" if i in _changes(tk) else "") for i in range(len(order)))

        row, tbl, hl = [], [], []
        for id_s in info["id_set"]["MPG"]["all"]:
            g = group[id_s]
            gi = g["info"]
            hexagonal = gi.crystal in ["trigonal", "hexagonal"]
            if hexagonal != (crystal == "hexagonal"):
                continue
            active = {X: [i[1:] for i in g["active_multipole"] if i[0] == X] for X in "QGTM"}
            row += [id_s.split(":")[1], "(" + gi.PG.split(":")[1] + ")", "", ""]
            ss = "$" + group[gi.PG]["info"].schoenflies + "$"
            tbl.append([r"\texttt{$" + gi.international + "$}", gi.type, gi.crystal, "Q"] + _active_list(active["Q"], order))
            tbl.append([ss, "", "", "G"] + _active_list(active["G"], order))
            tbl.append(["", "", "", "T"] + _active_list(active["T"], order))
            tbl.append(["", "", "", "M"] + _active_list(active["M"], order))
            hl.append(len(tbl) - 1)

        pdf.table(tbl, row, lbl, r"\#", cpos=cpos, long=True, stretch=1.6, hl=hl[:-1], caption=f"Active Multipole ({cap})")

    pdf.build()


# ==================================================
def create_pdf_info():
    """
    Create all PDFs of overall information in multipiegen/pdf/info.
    """
    os.makedirs(pdfdir, exist_ok=True)
    info = load_info()
    group = load_group()
    top_key = orbital_top_key(info)

    create_group_info_pg(info, group)
    create_group_info_sg(info, group)
    create_group_info_mpg(info, group)
    create_group_info_msg(info, group)
    create_compatibility_relation(info, group, top_key)
    create_active_multipole(info, group, top_key)


# ==================================================
if __name__ == "__main__":
    create_pdf_info()
