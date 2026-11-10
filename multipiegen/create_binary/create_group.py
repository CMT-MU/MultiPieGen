# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create group.pkl (symmetry operation, character, Wyckoff site/bond for all PG, SG, MPG and MSG).

    python -m multipiegen.create_binary.create_group [--serial]

Based on MultiPie (others/tools/create_binary/create_data_group.py, create_data_group_mapping.py,
create_data_group_bond.py), which are merged here. Differences from MultiPie:

    - the representative site/vector are those of the new root cluster (info.pkl).
    - "mapping" of symmetry operation is the mapping to the Oh/D6h cluster (data_op_mapping), for PG and SG.
    - Wyckoff bond has "symmetry" (bond symmetry).
    - complex PG ("-c") has all data, copied from the real PG except for the character.
    - all sympy values are stored in a simple form (simplify_value).
"""

import argparse
import logging
import os
import time
from multiprocessing import Pool

import numpy as np
import sympy as sp
from multipie import SOMatrixType
from multipie.util.util import normalize_vector, str_to_sympy
from multipie.util.util_binary import BinaryManager
from multipie.util.util_crystal import (
    convert_to_cartesian_hexagonal,
    convert_to_fractional_hexagonal,
    convert_to_fractional_hexagonal_matrix,
    convert_to_primitive,
    convert_to_primitive_matrix,
    t_dict,
)
from multipie.util.util_dict import Dict
from multipie.util.util_tag import TagIrrep

from multipiegen import __bin_dir__
from multipiegen.data.data_character_table import character, character_plus
from multipiegen.data.data_group_mpg import mpg_info
from multipiegen.data.data_group_msg import msg_info
from multipiegen.data.data_group_pg import pg_info
from multipiegen.data.data_group_sg import sg_info
from multipiegen.data.data_op_mapping import op_mapping_d6h_oh_cluster
from multipiegen.data.data_wyckoff_bond_sym import wyckoff_bond_data
from multipiegen.data.data_wyckoff_mpg import wyckoff_mpg_data
from multipiegen.data.data_wyckoff_msg import wyckoff_msg_data
from multipiegen.data.data_wyckoff_site import wyckoff_site_data
from multipiegen.util.util_group import group_file
from multipiegen.util.util_info import load_info
from multipiegen.util.util_symmetry_operation import (
    point_tag,
    simplify_array,
    simplify_value,
    symmetry_operation_matrix,
    symmetry_operation_matrix_multipole,
)

TOL = 1e-8

# ==================================================
h_group = """
* Information for all PG ("PG").
- PG_id (str): (dict) id, "PG:{ID}".
  - "info" (str): (namedtuple) info.
    - tag: Schoenflies in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: = id.
    - SG: associated SG, 1st SG in the same PGs.
    - MPG: associated MPG, 1st MPG with type II in the same PGs.
    - MSG: associated MSG, PG -> MPG -> MSG.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "generator" (str): ([str]) list of generator SO tag.
    - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
    - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
      - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
    - "product" (str): (dict) product table.
      - (SO_tag1,SO_tag2) (str,str): (str) SO_tag of product, SO1 SO2.
    - "mapping" (str): (dict) mapping of SO to the Oh/D6h cluster (data_op_mapping).
      - "top" (str): (str) top group, "Oh" or "D6h".
      - "mapping" (str): ([[int]]) mapping[coset][i] = SO number (from 1) of top group for i-th SO of group.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SOs of SG and PG (without translational part) are in the same order.
      - SO matrices are in the same order of "tag".
      - monopole (s=0) and dipole (s=1) basis are [1] and [x,y,z].
      - mapping[0] is the group itself, and mapping[k] (k>0) is the left coset g_k * group.
  - "wyckoff" (str): (dict) Wyckoff site/bond of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
        - "bond" (str): ([str]) list of bw_tag.
    - "bond" (str): (dict) info. of Wyckoff bond.
      - bw_tag (str): (dict) bond Wyckoff tag.
        - "symmetry" (str): (str) bond symmetry.
        - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional).
        - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
    NOTE:
      - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
      - Wyckoff sites are in the same order of ITA (BCS).
      - representative site/vector are scaled site0/vector0 of the root cluster (info["root_cluster"]).
  - "character" (str): (dict) character table of group.
    - "conjugacy" (str): ([[str]]) list of SO_tag in each conjugacy class.
    - "dimension" (str): (dict) dimension of irrep.
      - irrep_tag (str): (int) dimension.
    - "table" (str): (dict) character table (first SO in each conjugacy class).
      - irrep_tag (str): (ndarray(n,sympy) list of character.
    - "table_full" (str): (dict) character table for all SO.
      - irrep_tag (str): (ndarray(n,sympy) list of character.
    - "polar_axial_conversion" (str): (dict) conversion between polar and axial.
      - irrep_tag (str): (str) converted irrep_tag.
    - "symmetric_product" (str): (dict) symmetric product decomposition.
      - (irrep_tag,irrep_tag) (str,str): [(int,str)] list of (n, irrep_tag).
    - "anti_symmetric_product" (str): (dict) anti-symmetric product decomposition.
      - irrep_tag (str): ([(int,str)]) list of (n,irrep_tag).
    NOTE:
      - wp = exp(+2pi i/3), wm = exp(-2pi i/3) in "table", and they are substituted in "table_full".
      - for complex character, use PG_tag with "-c", otherwise without "-c", real version is given.
  NOTE:
    - complex PG ("-c") has the same "symmetry_operation" and "wyckoff" as the real PG.

* Information for all SG ("SG").
- SG_id (str): (dict) id, "SG:{ID}".
  - "info" (str): (namedtuple) info.
    - tag: Schoenflies in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting comment.
    - PG: associated PG, unique.
    - SG: = id.
    - MPG: associated MPG, SG -> PG -> MPG.
    - MSG: associated MSG, 1st MSG with type II in the same SGs.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "generator" (str): ([str]) list of generator SO tag.
    - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
    - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "cartesian_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, primitive).
    - "fractional_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, primitive).
    - "plus_set" (str): (ndarray(n,3,sympy)) list of plus set vector (fractional, conventional).
    - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
      - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
    - "mapping" (str): (dict) mapping of SO to the Oh/D6h cluster (that of the PG, data_op_mapping).
      - "top" (str): (str) top group, "Oh" or "D6h".
      - "mapping" (str): ([[int]]) mapping[coset][i] = SO number (from 1) of top group for i-th SO of group.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SOs of SG and PG (without translational part) are in the same order.
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site/bond of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
        - "primitive" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, primitive).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
        - "bond" (str): ([str]) list of bw_tag.
    - "bond" (str): (dict) info. of Wyckoff bond.
      - bw_tag (str): (dict) bond Wyckoff tag.
        - "symmetry" (str): (str) bond symmetry.
        - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, conventional, no plus set).
        - "primitive" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, primitive).
        - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
    NOTE:
      - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
      - Wyckoff sites are in the same order of ITA (BCS).

* Information for all MPG ("MPG").
- MPG_id (str): (dict) id, "MPG:{PG}.{no}.{ID}".
  - "info" (str): (namedtuple) info.
    - tag: MPG_id (PG.no.ID).
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, MPG -> MSG -> SG.
    - MPG: = id.
    - MSG: associated MSG, 1st MSG in the same MPGs.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III.
    - SO: symmetry operations.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
    - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
  - "active_multipole" (list): active MPs list, [str] (multipoles in info["response_tensor"]["multipole"]).

* Information for all MSG ("MSG").
- MSG_id (str): (dict) id, "MSG:{SG}.{no}".
  - "info" (str): (namedtuple) info.
    - tag: BNS_id (SG.no).
    - BNS: Belov-Neronova-Smirnova (international) notation in LaTeX.
    - OG: Opechowski-Guccione notation in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, unique.
    - MPG: associated MPG, unique.
    - MSG: = id.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III, IV.
    - SO: symmetry operations.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
    - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.

* all sympy values are expanded with rationalized denominators.
"""

_wp_wm = {sp.Symbol("wp"): -sp.S(1) / 2 + sp.sqrt(3) * sp.I / 2, sp.Symbol("wm"): -sp.S(1) / 2 - sp.sqrt(3) * sp.I / 2}


# ================================================== representative site/vector, basis functions (per process).
_shared = {}


def _info_data():
    """
    Data taken from info.pkl (loaded once per process).

    Returns:
        - (dict) -- "site"/"vector": {hexagonal(bool): ndarray(3)} (fractional), "multipole": {hexagonal: [str]},
          "basis_function": {name: sympy}.
    """
    if not _shared:
        info = load_info()
        rc = info["root_cluster"]
        _shared["site"] = {False: rc["oh_scaled_site0_frac"], True: rc["d6h_scaled_site0_frac"]}
        _shared["vector"] = {False: rc["oh_vector0_frac"], True: rc["d6h_vector0_frac"]}
        mp = info["response_tensor"]["multipole"]
        _shared["multipole"] = {False: list(mp["cubic"]), True: list(mp["hexagonal"])}
        bf = info["harmonics"]["basis_function"]
        _shared["basis_function"] = {name: bf[name][0] for name in mp["cubic"] + mp["hexagonal"]}
    return _shared


# ==================================================
def _sub(expr, s):
    """Substitute x,y,z (or X,Y,Z) in a string expression by the components of s."""
    for c, v in zip("xyz", s):
        expr = expr.replace(c, f"({v})")
    return expr


# ==================================================
def _is_zero_mod1(d):
    """Is float array d zero modulo 1 (all components) ?"""
    return bool(np.all(np.abs(d - np.round(d)) < TOL))


# ==================================================
def _real_tag(tag):
    """Real PG tag of complex PG ("-c")."""
    return tag.replace("-c", "")


# ==================================================
def _so_mapping(pg_tag, so):
    """
    Mapping of SO to the Oh/D6h cluster, through the point group (op_mapping_d6h_oh_cluster).

    Args:
        pg_tag (str): real PG tag (of the SG for SG).
        so (list): SO tags of the group (SG has the translational part, ignored).

    Returns:
        - (dict) -- {"top": "Oh"/"D6h", "mapping": [[int]]}.
    """
    top, mapping = op_mapping_d6h_oh_cluster[pg_tag]
    pg_so = pg_info[[k for k, v in pg_info.items() if v.tag == pg_tag][0]].SO
    idx = [pg_so.index(t.split(":")[0]) for t in so]
    return {"top": top, "mapping": [[c[i] for i in idx] for c in mapping]}


# ==================================================
def _product_table(tag, frac):
    """
    Product table of point group.

    Args:
        tag (list): SO tags.
        frac (ndarray): SO matrices (fractional, integer entries).

    Returns:
        - (dict) -- {(SO_tag1,SO_tag2): SO_tag of SO1 SO2}.
    """
    m = np.asarray(frac, dtype=float).round().astype(int)
    key = {mi.tobytes(): t for mi, t in zip(m, tag)}
    return {(t1, t2): key[(m1 @ m2).tobytes()] for t1, m1 in zip(tag, m) for t2, m2 in zip(tag, m)}


# ==================================================
def _matrices(so, hexagonal):
    """
    Cartesian and fractional SO matrices.

    Args:
        so (list): SO tags.
        hexagonal (bool): hexagonal ?

    Returns:
        - (ndarray) -- cartesian SO matrices.
        - (ndarray) -- fractional SO matrices.
        - (ndarray) -- determinants.
    """
    cart = np.asarray([symmetry_operation_matrix(i, hexagonal) for i in so])
    frac = simplify_array(convert_to_fractional_hexagonal_matrix(cart)) if hexagonal else cart.copy()
    det = simplify_array([sp.Matrix(i[0:3, 0:3]).det() for i in cart])
    return cart, frac, det


# ==================================================
def create_symmetry_operation(g_info, tag, pg_tag):
    """
    Symmetry operation of PG or SG.

    Args:
        g_info (namedtuple): PGInfoType or SGInfoType.
        tag (str): group tag (SG has "^").
        pg_tag (str): real PG tag (of the SG for SG), key of op_mapping_d6h_oh_cluster.

    Returns:
        - (dict) -- symmetry operation.
    """
    hexagonal = g_info.hexagonal_g
    so = g_info.SO
    cart, frac, det = _matrices(so, hexagonal)

    dic = {"tag": so, "generator": g_info.SO_gen, "cartesian": cart, "fractional": frac, "det": det}
    if tag.count("^") > 0:
        lattice = g_info.lattice
        dic["cartesian_primitive"] = simplify_array(convert_to_primitive_matrix(lattice, cart))
        dic["fractional_primitive"] = simplify_array(convert_to_primitive_matrix(lattice, frac))
        dic["plus_set"] = simplify_array(t_dict[lattice])
    else:
        dic["product"] = _product_table(so, frac)

    so_mat = Dict(SOMatrixType)
    for s in range(2):
        so_mat[("Q", s)] = [symmetry_operation_matrix_multipole(s, i, hexagonal, axial=False, tesseral=True) for i in so]
        so_mat[("G", s)] = [symmetry_operation_matrix_multipole(s, i, hexagonal, axial=True, tesseral=True) for i in so]
    dic["so_matrix"] = so_mat
    dic["mapping"] = _so_mapping(pg_tag, so)

    return dic


# ==================================================
def create_m_symmetry_operation(g_info):
    """
    Symmetry operation of MPG or MSG.

    Args:
        g_info (namedtuple): MPGInfoType or MSGInfoType.

    Returns:
        - (dict) -- symmetry operation.
    """
    so = g_info.SO
    cart, frac, det = _matrices(so, g_info.hexagonal_g)
    tr = np.asarray([-sp.S(1) if i.count("'") > 0 else sp.S(1) for i in so], dtype=object)
    return {"tag": so, "tr_sign": tr, "cartesian": cart, "fractional": frac, "det": det}


# ==================================================
def create_character(tag, so_tag):
    """
    Character of PG.

    Args:
        tag (str): PG tag (with "-c" for complex PG).
        so_tag (list): SO tags.

    Returns:
        - (dict) -- character.
    """
    conjugacy, table = character[tag]
    plus = character_plus[tag]

    cls = {s: i for i, c in enumerate(conjugacy) for s in c}
    table = {irrep: str_to_sympy(v) for irrep, v in table.items()}
    table_full = {}
    for irrep, ch in table.items():
        v = np.asarray([ch[cls[s]] for s in so_tag], dtype=object)
        table_full[irrep] = simplify_array([sp.S(i).subs(_wp_wm) for i in v])

    return {
        "conjugacy": conjugacy,
        "dimension": {i: TagIrrep.parse(i)["dimension"] for i in table.keys()},
        "table": table,
        "table_full": table_full,
        "polar_axial_conversion": plus["polar_axial"],
        "symmetric_product": plus["symmetric_product"],
        "anti_symmetric_product": plus["anti_symmetric_product"],
    }


# ==================================================
def create_wyckoff_site(g_info, wyckoff, sym_op, space_group, plus_set=True):
    """
    Wyckoff site.

    Args:
        g_info (namedtuple): group info.
        wyckoff (dict): Wyckoff data of group, {wp: (expression, symmetry)}.
        sym_op (dict): symmetry operation.
        space_group (bool): space group (translation part) ?
        plus_set (bool, optional): add plus set for space group (SG) ?

    Returns:
        - (dict) -- Wyckoff site.
    """
    hexagonal = g_info.hexagonal_g
    lattice = g_info.lattice
    s = _info_data()["site"][hexagonal]
    so = sym_op["fractional"]
    so_f = np.asarray(so, dtype=float)
    t_lst = sym_op["plus_set"] if space_group and plus_set else np.zeros((1, 3), dtype=object)
    t_f = np.asarray(t_lst, dtype=float)

    dic = {}
    for wp, (wp_data, sym) in wyckoff.items():
        pos = str_to_sympy(_sub(wp_data, s))
        pos_f = np.asarray(pos, dtype=float)

        # mapping by applying SO to the first site.
        if space_group:
            pw = (so_f @ np.append(pos_f[0], 1.0))[:, 0:3]
            mapping = [[no + 1 for no, q in enumerate(pw) if any(_is_zero_mod1(p - q - t) for t in t_f)] for p in pos_f]
            pos = np.mod(np.concatenate([pos + t for t in t_lst]), 1)
        else:
            pw = so_f @ pos_f[0]
            mapping = [[no + 1 for no, q in enumerate(pw) if np.all(np.abs(p - q) < TOL)] for p in pos_f]

            # normalize site in cartesian coordinate.
            cart = convert_to_cartesian_hexagonal(pos) if hexagonal else pos
            cart = normalize_vector(simplify_array(cart))
            pos = convert_to_fractional_hexagonal(simplify_array(cart)) if hexagonal else cart

        pos_ex = simplify_array(str_to_sympy(wp_data))
        d = {"symmetry": sym, "conventional": pos_ex}
        if space_group and plus_set:
            d["primitive"] = simplify_array(convert_to_primitive(lattice, pos_ex, shift=True))
        d["reference"] = simplify_array(pos)
        d["mapping"] = mapping
        dic[wp] = d

    return dic


# ==================================================
def create_wyckoff_bond(g_info, bond_data, site, sym_op, space_group):
    """
    Wyckoff bond (and bond list of each Wyckoff site, added to site).

    Args:
        g_info (namedtuple): group info.
        bond_data (dict): bond data of group, {bw: (vector, center, mapping, bond_symmetry)}.
        site (dict): Wyckoff site (bond list is added).
        sym_op (dict): symmetry operation.
        space_group (bool): space group ?

    Returns:
        - (dict) -- Wyckoff bond.
    """
    hexagonal = g_info.hexagonal_g
    lattice = g_info.lattice
    s = _info_data()["vector"][hexagonal]

    for v in site.values():
        v["bond"] = []

    dic = {}
    for bw, (vector, center, mapping, b_sym) in bond_data.items():
        s_wp = bw.split("@")[1]
        site[s_wp]["bond"].append(bw)
        n_site = len(site[s_wp]["conventional"])
        site_mp = {tuple(c): r for c, r in zip(site[s_wp]["conventional"], site[s_wp]["reference"][:n_site])}

        center = simplify_array(str_to_sympy(center))
        bond_v = simplify_array(str_to_sympy(_sub(vector.replace("X", "x").replace("Y", "y").replace("Z", "z"), s)))
        bond_c = np.asarray([site_mp[tuple(i)] for i in center], dtype=object)
        vector = simplify_array(str_to_sympy(vector))

        if space_group:
            plus_set = sym_op["plus_set"]
            bond_v = np.tile(bond_v, (len(plus_set), 1))
            bond_c = np.mod(np.concatenate([bond_c + t for t in plus_set]), 1)

        position = np.hstack((vector, center))
        reference = simplify_array(np.hstack((bond_v, bond_c)))

        d = {"symmetry": b_sym, "conventional": position}
        if space_group:
            prim_c = convert_to_primitive(lattice, center, shift=True)
            d["primitive"] = simplify_array(np.hstack((vector, prim_c)))
        d["reference"] = reference
        d["mapping"] = mapping
        dic[bw] = d

    return dic


# ==================================================
_bf_cache = {}


def _transformed_basis(name, ptag, hexagonal):
    """Basis function of multipole transformed by SO (cached)."""
    key = (name, ptag, hexagonal)
    if key not in _bf_cache:
        x, y, z = sp.symbols("x y z", real=True)
        bf = _info_data()["basis_function"][name]
        r = symmetry_operation_matrix(ptag, hexagonal).T @ np.array([x, y, z])
        _bf_cache[key] = sp.expand(bf.subs({x: r[0], y: r[1], z: r[2]}, simultaneous=True))
    return _bf_cache[key]


# ==================================================
def create_active_multipole(g_info, sym_op):
    """
    Active multipoles of MPG (invariant under all SOs).

    Args:
        g_info (namedtuple): MPGInfoType.
        sym_op (dict): symmetry operation.

    Returns:
        - (list) -- active multipoles, ["Q..."] + ["G..."] + ["T..."] + ["M..."].
    """
    hexagonal = g_info.hexagonal_g
    mp_list = _info_data()["multipole"][hexagonal]
    ptags = [point_tag(i) for i in sym_op["tag"]]
    det = [int(i) for i in sym_op["det"]]
    tr = [int(i) for i in sym_op["tr_sign"]]

    signs = {
        "Q": [1] * len(ptags),
        "G": det,
        "T": tr,
        "M": [d * t for d, t in zip(det, tr)],
    }
    active = {X: [] for X in signs}
    for name in mp_list:
        bf = sp.expand(_info_data()["basis_function"][name])
        tf = [_transformed_basis(name, p, hexagonal) for p in ptags]
        for X, sg in signs.items():
            if all(sp.expand(c * f - bf) == 0 for c, f in zip(sg, tf)):
                active[X].append(X + name)

    return active["Q"] + active["G"] + active["T"] + active["M"]


# ==================================================
def create_pg(pg_id, tag):
    """
    Data of PG.

    Args:
        pg_id (str): PG id.
        tag (str): PG tag.

    Returns:
        - (dict) -- data of PG.
    """
    rtag = _real_tag(tag)
    g_info = pg_info[pg_id]
    so = create_symmetry_operation(g_info, rtag, rtag)
    site = create_wyckoff_site(g_info, wyckoff_site_data[rtag], so, False)
    bond = create_wyckoff_bond(g_info, wyckoff_bond_data[rtag], site, so, False)
    return {
        "info": g_info,
        "symmetry_operation": so,
        "wyckoff": {"site": site, "bond": bond},
        "character": create_character(tag, so["tag"]),
    }


# ==================================================
def create_sg(sg_id, tag):
    """
    Data of SG.

    Args:
        sg_id (str): SG id.
        tag (str): SG tag.

    Returns:
        - (dict) -- data of SG.
    """
    g_info = sg_info[sg_id]
    so = create_symmetry_operation(g_info, tag, pg_info[g_info.PG].tag)
    site = create_wyckoff_site(g_info, wyckoff_site_data[tag], so, True)
    bond = create_wyckoff_bond(g_info, wyckoff_bond_data[tag], site, so, True)
    return {"info": g_info, "symmetry_operation": so, "wyckoff": {"site": site, "bond": bond}}


# ==================================================
def create_mpg(mpg_id, tag):
    """
    Data of MPG.

    Args:
        mpg_id (str): MPG id.
        tag (str): MPG tag.

    Returns:
        - (dict) -- data of MPG.
    """
    g_info = mpg_info[mpg_id]
    so = create_m_symmetry_operation(g_info)
    site = create_wyckoff_site(g_info, wyckoff_mpg_data[mpg_id.split(":")[1]], so, False)
    return {
        "info": g_info,
        "symmetry_operation": so,
        "wyckoff": {"site": site},
        "active_multipole": create_active_multipole(g_info, so),
    }


# ==================================================
def create_msg(msg_id, tag):
    """
    Data of MSG.

    Args:
        msg_id (str): MSG id.
        tag (str): MSG tag.

    Returns:
        - (dict) -- data of MSG.
    """
    g_info = msg_info[msg_id]
    so = create_m_symmetry_operation(g_info)
    site = create_wyckoff_site(g_info, wyckoff_msg_data[msg_id.split(":")[1]], so, True, plus_set=False)
    return {"info": g_info, "symmetry_operation": so, "wyckoff": {"site": site}}


_creator = {"PG": create_pg, "SG": create_sg, "MPG": create_mpg, "MSG": create_msg}


# ==================================================
def create_group_data(group_id, tag=None):
    """
    Data of a group.

    Args:
        group_id (str): group id, "PG:32", "SG:221", "MPG:...", "MSG:...".
        tag (str, optional): group tag (taken from info.pkl if None).

    Returns:
        - (dict) -- data of the group.
    """
    if tag is None:
        tag = load_info()["tag"][group_id]
    return _creator[group_id.split(":")[0]](group_id, tag)


def _worker(arg):
    return arg[0], create_group_data(*arg)


# ==================================================
def group_ids(info=None):
    """
    All group ids in group.pkl.

    Args:
        info (BinaryManager, optional): info.

    Returns:
        - (list) -- [(group id, tag)].
    """
    if info is None:
        info = load_info()
    ids = info["id_set"]
    lst = ids["PG"]["all"] + ids["PG"]["complex"] + ids["SG"]["all"] + ids["MPG"]["all"] + ids["MSG"]["all"]
    return [(i, info["tag"][i]) for i in lst]


# ==================================================
def create_group(ids=None, serial=False, verbose=True):
    """
    Create group (BinaryManager).

    Args:
        ids (list, optional): [(group id, tag)], all groups if None.
        serial (bool, optional): run in a single process ?
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- group.
    """
    if ids is None:
        ids = group_ids()

    group = BinaryManager(verbose=verbose, topdir=__bin_dir__)
    if serial:
        results = map(_worker, ids)
    else:
        pool = Pool(os.cpu_count())
        results = pool.imap(_worker, ids, chunksize=4)

    for gid, data in results:
        if verbose:
            logging.info(f"created {gid}")
        group[gid] = data

    if not serial:
        pool.close()
        pool.join()

    group.add_comment(h_group)
    return group


# ==================================================
def save_group(serial=False):
    """
    Create and save group.pkl in multipiegen/binary.

    Args:
        serial (bool, optional): run in a single process ?
    """
    t0 = time.time()
    create_group(serial=serial, verbose=False).save_binary(group_file)
    logging.info(f"group.pkl: {time.time() - t0:.1f} s")


# ==================================================
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    parser = argparse.ArgumentParser(description="create group.pkl")
    parser.add_argument("--serial", action="store_true", help="run in a single process")
    save_group(parser.parse_args().serial)
