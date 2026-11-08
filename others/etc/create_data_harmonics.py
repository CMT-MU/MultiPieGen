# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create binary data (information for harmonics in each groups).

Usage:
    python create_data_harmonics.py                       # create binary data (default).
    python create_data_harmonics.py --mapping             # print only mapping (no binary data needed).
    python create_data_harmonics.py --mapping -g Td D4h   # only selected point groups.
    python create_data_harmonics.py --mapping --renumber  # mapping after multiplicity renumbering.
    python create_data_harmonics.py --mapping -f json -o mapping.json
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import argparse
import json
import logging
import numpy as np
import sympy as sp
from collections import namedtuple

from multipie import PGMultipoleType, IrrepListType
from multipie.util.util import timer, str_to_sympy
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager

from multipiegen.data.data_top_harmonics import top_harmonics, top_harmonics_tag
from multipiegen.data.data_compatibility_relation import compatibility_relation

# ==================================================
h_harmonics = """
* Harmonics for all PG.
- PG_id (str): (dict) harmonics data for each PG.
  - "harmonics" (str): (Dict) polar (Q) and axial (G) harmonics (s=k=0,x="q").
    - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,sympy), ndarray(2l+1,dim,sympy), ndarray(dim,sympy)) [cartesian ex.], u-matrix, <m|gamma>, and [tesseral ex.] for each component.
  - "irrep" (str): (Dict) Q and G harmonics grouped by irrep.
    - (X,Gamma) (str,str): ([(int,int)]) list of (l,n).
  - "index" (str): (dict) harmonics index by (X,expression).
    - (X,ex) (str,str): ((X,l,Gamma,n,p,s,k,x),comp).
NOTE:
  - u-matrix is in order of m=[l, l-1, ..., -l][gamma].
  - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
  - for complex expresson, use point group tag with "-c".
  - multipoles are sorted in order of [Gamma, l, Q/G, n, p].
"""

# polar (Q) and axial (G).
HARMONICS_TYPES = ("Q", "G")


# ==================================================
# top harmonics.
# ==================================================
def create_top_harmonics():
    """
    Create top harmonics (polar up to rank 11 for Oh and D6h).

    Returns:
        - (dict) -- top harmonics data.

    Notes:
        - "Oh/D6h" (dict): data.
          {("Oh/D6h",l,irrep.,multiplicity,component): (tesseral ex., cartesian ex., u-matrix)} ({(str,int,str,int,int): (sympy,sympy,[sympy(2l+1)])}).
        - u-matrix <m|gamma> is in order of m=[l, l-1, ..., -l].
        - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
        - complex expressions (postfixed with a/b in irrep) are added for creating subgroup harmonics.
        - axial irrep is obtained by interchanging "g" and "u".
        - representation matrices of symmetry operations are common upto sign and order for the same irrep.
    """
    info = BinaryManager("info", topdir=BIN_DIR)
    x, y, z = info["harmonics"]["variable"]
    sub = {"x": x, "y": y, "z": z}

    dic_harm = {}
    for tag, (t_ex, c_ex, u_mat) in top_harmonics.items():
        t_ex = str_to_sympy(t_ex, subs=sub)
        c_ex = str_to_sympy(c_ex, subs=sub)
        u_mat = str_to_sympy(u_mat)
        dic_harm[tag] = (t_ex, c_ex, u_mat)
    dic = {"data": dic_harm, "tag": top_harmonics_tag}

    return dic


# ==================================================
# mapping (top group harmonics -> subgroup harmonics).
# ==================================================
def _top_irrep_comp(irrep, comp, sub_irrep):
    """
    Resolve (irrep, comp) of the top harmonics to be referred.

    If the irrep. of the subgroup is a complex one (contains "a"/"b"), the complex expression of the top group
    (irrep+"a" / irrep+"b", single component) is used instead of the real one.

    Args:
        irrep (str): irrep. of top group.
        comp (int): component of top group.
        sub_irrep (str): irrep. of subgroup (polar or axial).

    Returns:
        - (str, int) -- irrep. and component of top harmonics to be referred.
    """
    if "a" in sub_irrep or "b" in sub_irrep:
        if comp == 0:
            return irrep + "a", -1
        if comp == 1:
            return irrep + "b", -1
    return irrep, comp


# ==================================================
def create_mapping(group_tag):
    """
    Create mapping from top group harmonics to subgroup harmonics (no sympy data are needed).

    Args:
        group_tag (str): point group tag (with "-c" for complex expression).

    Returns:
        - (dict) -- mapping, {"Q": {...}, "G": {...}}.

    Notes:
        - mapping[X] = {(l, irrep_sub, comp_sub): [(coeff, (top, l, irrep_top, mul_top, comp_top)), ...]}
          ({(int,str,int): [(int, (str,int,str,int,int))]}).
        - the key is the harmonics of the subgroup, and the list is top group harmonics belonging to it.
        - the list order follows top_harmonics_tag.
        - multiplicity of the subgroup is not yet assigned (see renumber_multiplicity).
    """
    top, cr = compatibility_relation[group_tag]

    mapping = {X: {} for X in HARMONICS_TYPES}
    for rank, irrep, mul, comp in top_harmonics_tag[top]:
        polar, axial = cr[(irrep, comp)]
        for X, (sub_irrep, sub_comp, coeff) in zip(HARMONICS_TYPES, (polar, axial)):
            top_irrep, top_comp = _top_irrep_comp(irrep, comp, sub_irrep)
            key = (rank, sub_irrep, sub_comp)
            mapping[X].setdefault(key, []).append((coeff, (top, rank, top_irrep, mul, top_comp)))

    return mapping


# ==================================================
def create_mappings(group_tags=None):
    """
    Create mappings for several point groups.

    Args:
        group_tags (list, optional): point group tags. If None, all groups in compatibility_relation.

    Returns:
        - (dict) -- {group_tag: mapping}.
    """
    if group_tags is None:
        group_tags = list(compatibility_relation)
    return {tag: create_mapping(tag) for tag in group_tags}


# ==================================================
def renumber_multiplicity(mapping):
    """
    Assign multiplicity of the subgroup harmonics.

    Args:
        mapping (dict): mapping created by create_mapping.

    Returns:
        - (dict) -- {"Q": {...}, "G": {...}}, each is {(l, irrep, n, comp): (coeff, top_key)}.

    Notes:
        - n = -1 if only one top harmonics belongs to the key, otherwise n = 1, 2, ... in the list order.
    """
    renumbered = {X: {} for X in HARMONICS_TYPES}
    for X, mapping_X in mapping.items():
        for (rank, irrep, comp), data in mapping_X.items():
            if len(data) == 1:
                renumbered[X][(rank, irrep, -1, comp)] = data[0]
            else:
                for n, v in enumerate(data, start=1):
                    renumbered[X][(rank, irrep, n, comp)] = v
    return renumbered


# ==================================================
# create harmonics of subgroup from mapping.
# ==================================================
def _apply_coefficients(renumbered, top_data):
    """
    Multiply coefficients to top harmonics.

    Args:
        renumbered (dict): output of renumber_multiplicity.
        top_data (dict): top_harmonics["data"].

    Returns:
        - (dict) -- {X: {(l, irrep, n, comp): (tesseral ex., cartesian ex., u-matrix)}}.
    """
    harmonics_data = {X: {} for X in HARMONICS_TYPES}
    for X, data_X in renumbered.items():
        for key, (coeff, top_key) in data_X.items():
            t, c, u = top_data[top_key]
            harmonics_data[X][key] = (coeff * t, coeff * c, coeff * u)
    return harmonics_data


# ==================================================
def _make_irrep_list(harmonics_data):
    """
    Create irrep list, {(X, irrep): [(l, n)]}.
    """
    irrep_list = Dict(IrrepListType)
    for X, data_X in harmonics_data.items():
        for rank, irrep, mul, comp in data_X.keys():
            if comp == -1 or comp == 0:  # count once for each irrep.
                irrep_list[(X, irrep)] = irrep_list.get((X, irrep), []) + [(rank, mul)]
    return irrep_list.sort(("X", list(HARMONICS_TYPES)))


# ==================================================
def _make_harmonics(harmonics_data, group_tag):
    """
    Create harmonics, {(X, l, Gamma, n, p, s, k, x): (cartesian ex., u-matrix, tesseral ex.)}.
    """
    tmp = Dict(namedtuple("tmp", ["l", "X", "irrep", "n", "comp"]))
    for X in HARMONICS_TYPES:
        for (rank, irrep, mul, comp), v in harmonics_data[X].items():
            tmp[(rank, X, irrep, mul, comp)] = v
    tmp = tmp.sort("l", ("X", list(HARMONICS_TYPES)), "irrep", "n", "comp")

    # gather components.
    grouped = Dict(namedtuple("Harmonics", ["l", "X", "irrep", "n"]))
    for (l, X, irrep, mul, comp), v in tmp.items():
        grouped[(l, X, irrep, mul)] = grouped.get((l, X, irrep, mul), []) + [v]

    harmonics = Dict(PGMultipoleType)
    for (l, X, Gamma, n), comps in grouped.items():
        tset = np.asarray([t for t, _, _ in comps])
        cset = [c for _, c, _ in comps]
        uset = np.asarray([u for _, _, u in comps]).T
        if "-c" not in group_tag:  # to avoid factorization for complex point group due to extremely slow.
            cset = np.vectorize(sp.factor)(np.asarray(cset))
        harmonics[(X, l, Gamma, n, -1, 0, 0, "q")] = (cset, uset, tset)

    return harmonics.sort("Gamma", "l", ("X", list(HARMONICS_TYPES)))


# ==================================================
def expand_component(basis, vector):
    ex_basis = {}
    for idx, v in basis.items():
        if len(v[0]) == 1:
            if vector:
                ex_basis[(idx, -1)] = (v[0][0], v[1][0], v[2][0])
            else:
                ex_basis[(idx, -1)] = (v[0][0], v[1][0])
        else:
            if vector:
                for no, (b, ex, ex1) in enumerate(zip(v[0], v[1], v[2])):
                    ex_basis[(idx, no)] = (b, ex, ex1)
            else:
                for no, (b, ex) in enumerate(zip(v[0], v[1])):
                    ex_basis[(idx, no)] = (b, ex)

    return ex_basis


# ==================================================
def _make_index(harmonics):
    """
    Create reverse dict, {(X, expression): (X, l, Gamma, n, p, s, k, x), comp}.
    """
    index = {}
    for X in HARMONICS_TYPES:
        q_harm = expand_component(harmonics.select(X=X), vector=False)
        for idx, val in q_harm.items():
            v = str(sp.expand(val[0])).replace(" ", "")
            index[(X, v)] = idx
    return index


# ==================================================
def create_group_harmonics(group_tag, top_harm=None):
    """
    Create harmonics data of a point group.

    Args:
        group_tag (str): point group tag.
        top_harm (dict, optional): output of create_top_harmonics (created if None; pass it to avoid recreation).

    Returns:
        - (dict) -- {"harmonics": ..., "irrep": ..., "index": ...}.
    """
    if top_harm is None:
        top_harm = create_top_harmonics()

    mapping = create_mapping(group_tag)
    renumbered = renumber_multiplicity(mapping)
    harmonics_data = _apply_coefficients(renumbered, top_harm["data"])

    irrep_list = _make_irrep_list(harmonics_data)
    harmonics = _make_harmonics(harmonics_data, group_tag)
    index = _make_index(harmonics)

    return {"harmonics": harmonics, "irrep": irrep_list, "index": index}


# ==================================================
@timer
def create_harmonics():
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)
    top_harm = create_top_harmonics()  # create only once.

    harmonics = BinaryManager(verbose=True, topdir=BIN_DIR)
    for no in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        g_info = group[no]
        group_tag = g_info["info"].tag
        print("creating", group_tag, flush=True)
        harmonics[no] = create_group_harmonics(group_tag, top_harm)

    harmonics.add_comment(h_harmonics)
    harmonics.save_binary("harmonics")


# ==================================================
# mapping output.
# ==================================================
def _sources(val):
    """
    Normalize to list of (coeff, top_key) (raw mapping: list, renumbered mapping: single tuple).
    """
    return val if isinstance(val, list) else [val]


# ==================================================
def format_mapping(group_tag, mapping, renumbered=False):
    """
    Format mapping as text.

    Args:
        group_tag (str): point group tag.
        mapping (dict): output of create_mapping (or renumber_multiplicity if renumbered=True).
        renumbered (bool, optional): mapping is renumbered one?

    Returns:
        - (str) -- text.
    """
    key_name = "(l, irrep, n, comp)" if renumbered else "(l, irrep, comp)"
    lines = [f"[{group_tag}] top={compatibility_relation[group_tag][0]}, key={key_name} <- coeff * (top, l, irrep, n, comp)"]
    for X in HARMONICS_TYPES:
        lines.append(f"  {X}:")
        for key, val in mapping[X].items():
            srcs = [f"{c:+d} * {tk}" for c, tk in _sources(val)]
            lines.append(f"    {key} <- " + ", ".join(srcs))
    return "\n".join(lines)


# ==================================================
def mapping_to_jsonable(mapping):
    """
    Convert mapping to JSON-serializable dict (tuple keys are converted to "a,b,c" strings).
    """
    return {
        X: {",".join(map(str, key)): [{"coeff": c, "top": list(tk)} for c, tk in _sources(val)] for key, val in mapping_X.items()}
        for X, mapping_X in mapping.items()
    }


# ==================================================
def output_mapping(group_tags=None, renumber=False, fmt="text", output=None):
    """
    Output only mapping (binary data are not touched).

    Args:
        group_tags (list, optional): point group tags (all if None).
        renumber (bool, optional): output mapping after multiplicity renumbering?
        fmt (str, optional): "text" or "json".
        output (str, optional): output file (stdout if None).
    """
    if group_tags is None:
        group_tags = list(compatibility_relation)
    unknown = [tag for tag in group_tags if tag not in compatibility_relation]
    if unknown:
        raise SystemExit(f"unknown group tag: {unknown}\navailable: {list(compatibility_relation)}")

    mappings = create_mappings(group_tags)
    if renumber:
        mappings = {tag: renumber_multiplicity(m) for tag, m in mappings.items()}

    if fmt == "json":
        text = json.dumps({tag: mapping_to_jsonable(m) for tag, m in mappings.items()}, indent=2)
    else:
        text = "\n\n".join(format_mapping(tag, m, renumber) for tag, m in mappings.items())

    if output is None:
        print(text)
    else:
        with open(output, "w", encoding="utf-8") as f:
            f.write(text + "\n")


# ================================================== main
def main():
    parser = argparse.ArgumentParser(description="Create harmonics binary data, or output only mapping.")
    parser.add_argument("--mapping", action="store_true", help="output only mapping (binary data are not created).")
    parser.add_argument("-g", "--group", nargs="+", metavar="TAG", help="point group tags for --mapping (default: all).")
    parser.add_argument("--renumber", action="store_true", help="output mapping after multiplicity renumbering.")
    parser.add_argument("-f", "--format", choices=("text", "json"), default="text", help="output format for --mapping.")
    parser.add_argument("-o", "--output", help="output file for --mapping (default: stdout).")
    args = parser.parse_args()

    logging.basicConfig(format="%(message)s", level=logging.INFO)
    if args.mapping:
        output_mapping(args.group, args.renumber, args.format, args.output)
    else:
        create_harmonics()


if __name__ == "__main__":
    main()
