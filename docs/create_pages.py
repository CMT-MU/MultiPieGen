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
    - docs/src/*_ja.md: the Japanese versions of the pages above (lang="ja"), as orphan pages.
      The fixed prose is translated (the strings are kept side by side with the English ones below); text taken from
      the data (comments, docstrings, headers of the binaries) is left in English.

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
crystal_name_ja = {
    "triclinic": "三斜晶",
    "monoclinic": "単斜晶",
    "orthorhombic": "直方晶",
    "tetragonal": "正方晶",
    "trigonal": "三方晶",
    "hexagonal": "六方晶",
    "cubic": "立方晶",
}
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
label_ja = {
    "harmonics": "調和関数",
    "character_table": "指標",
    "symmetry_operation": "対称操作",
    "wyckoff_site": "Wyckoff サイト",
    "wyckoff_bond": "Wyckoff ボンド",
    "representation_matrix": "表現行列",
    "response_tensor_Q": "Q テンソル",
    "atomic_multipole_lg": "原子 lg",
    "atomic_multipole_lgs": "原子 lgs",
    "atomic_multipole_jml": "原子 jml",
    "response_tensor_G": "G テンソル",
    "response_tensor_T": "T テンソル",
    "response_tensor_M": "M テンソル",
}  # the others (s=1 q, ...) as in English.
column_ja = {"No.": "番号", "international": "国際記号", "tag": "タグ", "type": "タイプ"}  # the others as in English.

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
head_ja = {  # (title, body).
    "PG": (
        "点群",
        """すべての 47 点群です：結晶点群 32、別の設定（"-1"）5、複素の群（"-c"）10。

- 調和関数：ランク 11 までの調和関数。極性（Q）と軸性（G）の記号を並べて示します（式は極性の調和関数のものです）。
- 指標：指標表、極性と軸性の変換、対称積と反対称積。
- 対称操作、Wyckoff サイトとボンド、表現行列：名前のとおりです。
- s=1..3 q/g：内部ランク s で、極性（q）または軸性（g）の内部基底をもつ点群多極子調和関数（ランク <= 3）。
- 原子 lg/lgs/jml：(L,gamma)、(L,gamma,s)、(J,M;L) 基底での原子 SAMB（s, p, d, f 軌道の bra-ket ブロックの行列）。
""",
    ),
    "SG": (
        "空間群",
        """すべての 230 空間群の対称操作、Wyckoff サイト、Wyckoff ボンドです。

番号、Hermann-Mauguin 記号、対称操作、Wyckoff 位置は、International Tables for Crystallography, Vol. A（ITA）および Bilbao Crystallographic Server と同じです。各空間群について、標準の設定（Bilbao Crystallographic Server の既定の設定）だけを示します：

- 単斜晶の群：unique axis b、cell choice 1。
- 原点の選び方が 2 つある群：origin choice 2（反転中心が原点）。
- 菱面体晶の群：六方晶の軸。
""",
    ),
    "MPG": (
        "磁気点群",
        """すべての 157 磁気点群の対称操作、Wyckoff サイト、応答テンソル（Q、G、T、M、ランク 4 まで）です。タイプ II（灰色）の群では、磁気的なテンソル（T、M）は示しません。

ID "a.b.c"：

- a：点群の番号（PG の列）。結晶点群は 1-32、別の設定は 33-37 です。
- b：点群の中での番号。タイプ I（点群そのもの）は 1、タイプ II（時間反転 1' を含む灰色の群）は 2、タイプ III（白黒の群）は 3 以上です。
- c：磁気点群の通し番号、1-157。

タイプ：I（通常）、II（灰色）、III（白黒。一部の操作が時間反転と組み合わさる）。
""",
    ),
    "MSG": (
        "磁気空間群",
        """すべての 1651 磁気空間群の対称操作と Wyckoff サイトです。各磁気空間群について、BNS 表記の標準の設定だけを示します。

ID "a.b" は BNS 番号です：

- a：空間群の番号。
- b：各晶族の中での通し番号。三斜晶、単斜晶、直方晶、正方晶、三方晶と六方晶、立方晶の群で、それぞれ 1 から始まります。

BNS と OG の列は、Belov-Neronova-Smirnova 記号と Opechowski-Guccione 記号です。タイプ：I（通常）、II（灰色）、III（白黒。一部の操作が時間反転と組み合わさる）、IV（白黒。時間反転が並進と組み合わさる）。
""",
    ),
}


# ==================================================
def _ja_head(title, en_file, front=None):
    """
    Beginning of a Japanese page: front matter (orphan, and the given keys), title, and the link to the English page.
    """
    return ["---"] + (front or []) + ["orphan: true", "---", "", f"# {title}", "", f"英語版は[こちら]({en_file})です。"]


def _page_name(name, lang):
    """File name of a page, e.g. "data.md" -> "data_ja.md"."""
    return name if lang == "en" else name[:-3] + f"_{lang}.md"


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
def create_group_pages(lang="en"):
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
        if lang == "en":
            cols += [label[k] for k in names]
            lines = [head[tp]]
        else:
            cols = [column_ja.get(c, c) for c in cols] + [label_ja.get(k, label[k]) for k in names]
            title, body = head_ja[tp]
            lines = ["\n".join(_ja_head(title, f"{tp.lower()}_pdf.md")) + "\n\n" + body]
        for crystal in crystal_name:
            lst = [i for i in ids if group[i]["info"].crystal == crystal]
            lines += [
                f"\n## {crystal.capitalize() if lang == 'en' else crystal_name_ja[crystal]}\n",
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
        write(_page_name(f"{tp.lower()}_pdf.md", lang), lines)


# ==================================================
example_desc_ja = {  # Japanese descriptions of the examples (the English ones are in docs/examples/README.md).
    "cluster_basis.py": "Oh/D6h の正規直交な（スカラーの）クラスター基底（cluster_scalar_basis.pkl、キーは PGMultipoleType）"
    "を読み込み、クラスター上の関数を展開します。",
    "cluster_vector_basis.py": "Oh/D6h のクラスターの正規直交なベクトルクラスター基底（内部ランク s=1、極性・軸性の内部基底、"
    "キーは PGMultipoleType）を読み込み、ベクトル場を展開します。",
    "bond_cluster_basis.py": "Oh/D6h のクラスターの正規直交なボンドクラスター基底（ベクトルクラスター基底をボンドベクトルと縮約し、"
    "正規直交化し直したもの。キーは MultiPie Ver.2 の bond_a と同じく PGMultipoleType）を読み込みます。",
    "cluster_samb.py": "実の点群と空間群の、各 Wyckoff サイトとボンドのクラスター SAMB（サイト；ボンドは対称な Q と反対称な T/M）"
    "です（cluster_samb.pkl）。",
    "group.py": "各群のデータ（group.pkl）です：対称操作、積、Oh/D6h への対応、指標、Wyckoff サイトとボンド、活性な多極子。",
    "representation_matrix.py": "ランク 11 までの調和関数から求めた点群の表現行列と、ある既約表現のすべての調和関数が"
    "同じ行列を与えることの確認です。",
    "response_tensor.py": "磁気点群の応答テンソル（ランク 0-4）を活性な多極子で表したものです（response_tensor.pkl）。",
    "harmonics.py": "点群の調和関数の式と u-matrix です（harmonics.pkl）。どのバイナリの多極子のキーからでも引けます。",
    "multipolar_harmonics.py": "data_basis の軌道を極性・軸性の内部基底とする点群多極子調和関数（内部ランク s=1,2,3、ランク l<=3）"
    "です（multipolar_harmonics.pkl）。",
    "atomic_multipole.py": "原子多極子の行列（atomic_multipole.pkl）と、jml/lgs/lg 基底での点群の原子 SAMB"
    "（キーは PGMultipoleType）です（atomic_multipole_group.pkl）。",
    "harmonics_symbol.py": "47 点群すべてについて、調和関数の極性と軸性の記号を並べたものです（調和関数の PDF と同じ）。",
    "info.py": "群の全体の情報です（info.pkl）：群の id/タグ/名前、軌道基底と基底関数、ルートクラスター。",
    "binary_data.py": "MultiPie 向けのバイナリデータ（binary_data）です：info.pkl と群ごとの 1 ファイル。BinaryManager だけで"
    "読めます。PG のファイルで SAMB のキーの調和関数を引く方法も示します。",
}


def create_examples(lang="en"):
    """
    The page of the examples, as a notebook (MyST-NB) executed when the docs are built: the code of each example after
    its docstring, split into cells at the "# ===..." lines, with the output below each cell.
    """
    front = [
        "jupytext:",
        "  text_representation:",
        "    format_name: myst",
        "kernelspec:",
        "  display_name: Python 3",
        "  language: python",
        "  name: python3",
    ]
    if lang == "en":
        lines = ["---"] + front + ["---", "", "# Examples"]
        lines += [
            "",
            "Sample code for each feature, in `docs/examples`. Run each file from anywhere after `pip install -e .`, e.g.",
            "`python docs/examples/group.py`. The outputs below are made when the documentation is built.",
        ]
    else:
        lines = _ja_head("例", "examples.md", front)
        lines += [
            "",
            "機能ごとのサンプルコードです（`docs/examples`）。`pip install -e .` の後、各ファイルはどこからでも実行できます。"
            "例：`python docs/examples/group.py`。下の出力は、ドキュメントのビルド時に作られます。",
        ]
    readme = open(os.path.join(doc_dir, "examples", "README.md")).read()
    for name, desc in re.findall(r"\| \[(\w+\.py)\]\(\S+\) \| (.+) \|", readme):
        if lang == "en":
            desc = f"{desc[0].upper() + desc[1:]}."
        else:
            assert name in example_desc_ja, f"no Japanese description of {name} in docs/create_pages.py"
            desc = example_desc_ja[name]
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
        lines += ["", f"## {name}", "", desc]
        for cell in cells:
            cell = "\n".join(cell).strip("\n")
            if cell:
                lines += ["", "```{code-cell} ipython3", cell, "```"]
    write(_page_name("examples.md", lang), lines)


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
data_content_ja = {
    "data_basis.py": "軌道と基底の定義：原子軌道、スピン、LM/JM/LMs/Lgs 基底、bra-ket ブロック、LaTeX の名前、Wannier90 の"
    "別名、立方晶・六方晶の軌道（点群多極子調和関数の内部基底）、調和関数のキー、最大ランクとランクの名前",
    "data_character_table.py": "47 点群すべての指標表",
    "data_compatibility_relation.py": "Oh と D6h の部分群の適合関係",
    "data_group_pg.py": "点群の情報",
    "data_group_sg.py": "空間群の情報",
    "data_group_mpg.py": "磁気点群の情報",
    "data_group_msg.py": "磁気空間群の情報",
    "data_parent_pg.py": "各点群の親の点群",
    "data_parent_sg.py": "各空間群の親の空間群",
    "data_response_tensor.py": "応答テンソルに現れる多極子と、そのデカルト成分",
    "data_root_cluster.py": "代表サイトと代表ベクトルをもつ Oh と D6h のクラスター",
    "data_top_harmonics.py": "トップの群（Oh、D6h）の調和関数",
    "data_wyckoff_site.py": "点群と空間群の Wyckoff サイト",
    "data_wyckoff_bond_sym.py": "点群と空間群の Wyckoff ボンドと、ボンドの対称性",
    "data_wyckoff_mpg.py": "磁気点群の Wyckoff サイト",
    "data_wyckoff_msg.py": "磁気空間群の Wyckoff サイト",
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
data_var_note_ja = {  # the others (structures of the data) as in English.
    ("data_basis.py", "atomic_lgs_basis"): "原子軌道の標準のスピンあり基底の並び（Lgs）",
    ("data_basis.py", "cubic_basis"): "cubic_orbital_basis のスピンあり基底の並び",
    ("data_basis.py", "hexagonal_basis"): "hexagonal_orbital_basis のスピンあり基底の並び",
    ("data_top_harmonics.py", "top_harmonics_tag"): "トップの調和関数のタグ",
    ("data_top_harmonics.py", "top_harmonics"): "Oh と D6h のトップの調和関数（デカルト座標の式）",
    ("data_wyckoff_site.py", "wyckoff_site_data"): "{ PG/SG: { Wyckoff: (position, site-symmetry) } }（点群が先）",
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
binary_content_ja = {
    "harmonics.pkl": "多極子のキーで引く、全 PG の調和関数（ランク <= 11）の式と u-matrix",
    "info.pkl": "全群の全体の情報、調和関数、基底、ルートクラスター、応答テンソル",
    "group.pkl": "PG、SG、MPG、MSG それぞれのデータ：対称操作、指標、Wyckoff サイトとボンド",
    "response_tensor.pkl": "全 MPG の応答テンソル（ランク 0-4）",
    "multipolar_harmonics.pkl": "全 PG の点群多極子調和関数（s=1-3、l<=3）",
    "atomic_multipole.pkl": "LMs/JML/LM 基底での原子多極子の行列",
    "atomic_multipole_group.pkl": "JML/Lgs/Lg 基底での全 PG の原子 SAMB",
    "cluster_scalar_basis.pkl": "Oh と D6h のクラスターの正規直交な（スカラーの）クラスター基底",
    "cluster_vector_basis.pkl": "Oh と D6h のクラスターの正規直交なベクトルクラスター基底（s=1）",
    "cluster_bond_basis.pkl": "Oh と D6h のクラスターの正規直交なボンドクラスター基底",
    "cluster_samb.pkl": "実の PG と SG の各 Wyckoff サイトとボンドのクラスター SAMB",
}


def _sentence(text, lang):
    """A description (noun phrase) as a sentence."""
    return text[0].upper() + text[1:] + "." if lang == "en" else text + "です。"


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
data_text = {
    "en": {
        "origin": {True: "MultiPie Ver.2.4.0", False: "new"},
        "created": "Created by {}.",
        "var_head": "| variable | type (size) | description |",
        "no_note": "see the description below",
    },
    "ja": {
        "origin": {True: "MultiPie Ver.2.4.0", False: "新規"},
        "created": "{} で作ります。",
        "var_head": "| 変数 | 型（サイズ） | 説明 |",
        "no_note": "下の説明を参照",
    },
}
data_intro_ja = [
    "`multipiegen/data` のデータテーブルは、生成するすべてのデータの元です。データテーブルは手で管理しています。一部は、"
    "最初に `others/creation` のスクリプトで作ります（スクリプトは自分と同じ場所に書き出し、その結果を `multipiegen/data` に"
    "コピーします）。「由来」は、MultiPie Ver.2.4.0 から取ったテーブル（一部は変更あり）では「MultiPie Ver.2.4.0」、"
    "それ以外では「新規」です。MultiPie Ver.2.4.0 から取ったテーブルは一度だけコピーし、それ以降はここで管理しています。"
    "生成では MultiPie Ver.2 のデータを読むことはなく、その一般的な関数とユーティリティだけを使います。",
    "",
    "これらのテーブルから導けるデータは、保存せずに `multipiegen` の中で計算します：",
    "",
    "| データ | 計算する場所 | 元のデータ |",
    "|---|---|---|",
    "| `*_tag`, `*_id`, `*_id_set`, `group_name_list` | `util/util_group_data` | `data_group_pg/sg/mpg/msg` の `*_info` |",
    "| `top_harmonics_symmetry`, `top_harmonics_u_matrix` | `util/util_top_harmonics` | `data_top_harmonics`"
    "（方域の部分は点群多極子調和関数から） |",
    "| `data_harmonics_mapping` | `util/util_harmonics_mapping` | `data_compatibility_relation`, "
    "`data_top_harmonics.top_harmonics_tag` |",
    "| `op_mapping_d6h_oh_cluster`, `op_mapping_parent_sg` | `util/util_op_mapping` | `data_parent_pg`, `data_parent_sg` |",
    "| 極性と軸性の変換、対称積と反対称積 | `create_binary/create_group.character_plus` | `data_character_table` と対称操作 |",
    "",
    "| ファイル | 内容 | 由来 | 作成 | 使用 |",
    "|---|---|---|---|---|",
]


def create_data_page(lang="en"):
    creators, users = _creators(), _users()
    ddir = os.path.join(__top_dir__, "multipiegen", "data")
    files = sorted(f for f in os.listdir(ddir) if f.startswith("data_") and f.endswith(".py"))
    txt = data_text[lang]
    content = data_content if lang == "en" else data_content_ja
    var_note = data_var_note if lang == "en" else {**data_var_note, **data_var_note_ja}
    if lang == "ja":
        lines = _ja_head("データテーブル", "data.md") + [""] + data_intro_ja
    else:
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
        assert f in content, f"no description ({lang}) of {f} in docs/create_pages.py"
        created = ", ".join(f"`{c}`" for c in creators.get(f, [])) or "-"
        used = ", ".join(f"`{u}`" for u in users.get(name, [])) or "-"
        origin = txt["origin"][f in data_from_multipie]
        lines.append(f"| [`{f}`](#{name.replace('_', '-')}-py) | {content[f]} | {origin} | {created} | {used} |")

    for f in files:
        name = f[:-3]
        src = open(os.path.join(ddir, f)).read()
        tree = ast.parse(src)
        doc = ast.get_docstring(tree)
        mod = importlib.import_module(f"multipiegen.data.{name}")
        lines += ["", f"## {f}", "", _sentence(content[f], lang)]
        created = creators.get(f, [])
        if created:
            lines += ["", txt["created"].format(", ".join(f"`{c}`" for c in created))]
        lines += ["", txt["var_head"], "|---|---|---|"]
        seen = set()
        for node in tree.body:
            if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                v = node.targets[0].id
                if v in seen:
                    continue
                seen.add(v)
                note = var_note.get((f, v)) or _var_comment(src.splitlines(), node.lineno) or _doc_note(doc, v) or txt["no_note"]
                lines.append(f"| `{v}` | {_size(getattr(mod, v))} | {note} |")
        if doc:
            lines += ["", "```text", doc, "```"]
    write(_page_name("data.md", lang), lines)


# ==================================================
binary_text = {
    "en": {
        "keys": "{} keys: ",
        "create": "create",
        "load": "load",
        "top": "top-level keys",
        "size": "size",
        "header": "Description of the dictionary keys (header of the binary):",
    },
    "ja": {
        "keys": "全 {} 個、",
        "create": "作成",
        "load": "読み込み",
        "top": "最上位のキー",
        "size": "サイズ",
        "header": "辞書のキーの説明（バイナリのヘッダ）：",
    },
}
binary_intro_ja = [
    "`multipiegen/binary` のバイナリデータは、`multipiegen/create_binary` のモジュールがデータテーブルから作ります。"
    "MultiPie Ver.2 の `BinaryManager` の形式（ヘッダ付きの gzip 圧縮した pickle）なので、MultiPieGen なしで読めます。"
    '例：`BinaryManager("info", topdir=...)`。`str(data)` でヘッダ、つまり下の各節に示す辞書のキーの説明が表示されます。',
    "",
    "| ファイル | 内容 | 最上位のキー | 作成 | 読み込み |",
    "|---|---|---|---|---|",
]


def create_binary_page(lang="en"):
    txt = binary_text[lang]
    sep = ": " if lang == "en" else "："
    if lang == "ja":
        lines = _ja_head("バイナリデータ", "binary.md") + [""] + binary_intro_ja
    else:
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
        if lang == "ja":
            content = binary_content_ja[f]
        name = f[:-4]
        bm = BinaryManager(name, topdir=__bin_dir__)
        keys = list(bm.keys())
        if len(keys) <= 10:
            top = ", ".join(f"`{k}`" for k in keys)
        else:
            top = txt["keys"].format(len(keys)) + ", ".join(f"`{k}`" for k in keys[:2]) + f", ..., `{keys[-1]}`"
        anchor = name.replace("_", "-") + "-pkl"
        cmd = f"`python -m multipiegen.create_binary.{create}`"
        ld = f"`multipiegen.util.{load}()`"
        lines.append(f"| [`{f}`](#{anchor}) | {content} | {top} | {cmd} | {ld} |")
        size = os.path.getsize(os.path.join(__bin_dir__, f))
        sections += [
            "",
            f"## {f}",
            "",
            _sentence(content, lang),
            "",
            f"- {txt['create']}{sep}{cmd}",
            f"- {txt['load']}{sep}{ld}",
            f"- {txt['top']}{sep}{top}",
            f"- {txt['size']}{sep}{size / 1024:,.0f} KB",
            "",
            txt["header"],
            "",
            "```text",
            str(bm).strip().strip('"').strip(),
            "```",
        ]
    write(_page_name("binary.md", lang), lines + sections + binary_data_section(lang))


# ==================================================
binary_data_intro_ja = [
    "",
    "## MultiPie 向けのバイナリデータ（binary_data）",
    "",
    "`multipiegen/binary_data` には、MultiPie Ver.3 以降に組み込むバイナリデータを、MultiPie の `multipie/binary_data` と"
    "同じ構成で置いています：`info.pkl`（`binary/info.pkl` と同じ）と、群ごとのファイル `{PG,SG,MPG,MSG}/{id}.pkl` です。"
    "これらは、上のバイナリを群ごとに並べ直して（PG の表現行列を追加します）、"
    "`python -m multipiegen.create_binary.create_binary_data` で作ります。"
    'ファイルは `BinaryManager("PG/PG_32", topdir=...)`（ファイル名では id の `:` を `_` に置き換えます）か '
    '`multipiegen.util.util_binary_data.load_binary_data("PG:32")` で読み込みます。'
    "SAMB と点群多極子調和関数は、調和関数の式を持ちません。キーの `(X, l, Gamma, n)` で、PG のファイルの "
    '`"harmonics"` から引いてください（T -> Q、M -> G）。SG と MPG では、`info.PG` の PG のファイルを使います。',
    "",
    "| 群 | ファイル | キー | サイズ |",
    "|---|---|---|---|",
]
binary_data_text = {
    "en": {
        "complex": "The keys of the complex PG (`-c`) lack `cluster_samb`.",
        "files": "### {} files",
        "header": "Description of the dictionary keys (header of `{}.pkl`):",
    },
    "ja": {
        "complex": "複素の PG（`-c`）のキーには `cluster_samb` がありません。",
        "files": "### {} のファイル",
        "header": "辞書のキーの説明（`{}.pkl` のヘッダ）：",
    },
}


def binary_data_section(lang="en"):
    """Section of binary.md on the binary data for MultiPie (multipiegen/binary_data)."""
    from multipiegen.util.util_binary_data import binary_data_dir, group_file, group_types, load_binary_data

    txt = binary_data_text[lang]
    lines = (
        list(binary_data_intro_ja)
        if lang == "ja"
        else [
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
    )
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
    lines += ["", txt["complex"], ""]
    for tp, (gid, bm) in examples.items():
        lines += [
            txt["files"].format(tp),
            "",
            txt["header"].format(group_file(gid)),
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
    for lang in ["en", "ja"]:
        create_group_pages(lang)
        create_examples(lang)
        create_data_page(lang)
        create_binary_page(lang)
