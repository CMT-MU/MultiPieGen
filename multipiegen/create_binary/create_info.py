# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create binary data of overall information of groups (info.pkl).
    - multipiegen/binary/info.pkl (load by multipiegen.util.util_info.load_info).

Based on MultiPie's others/tools/create_binary/create_data_info.py, with
    - data_definition replaced by data_basis (orbital tags renamed), data_harmonics, data_root_cluster,
    - basis functions taken from the top harmonics (data_top_harmonics_info) through data_basis.harmonics_key,
    - orbital decomposition (orbital_decomposition) removed,
    - group name list (data_group_name_list) added.

Run:  python -m multipiegen.create_binary.create_info
"""

import numpy as np
import sympy as sp
from multipie.util.util import normalize_vector, str_to_sympy
from multipie.util.util_binary import BinaryManager
from multipiegen import __bin_dir__
from multipiegen.data.data_basis import (
    atomic_orbital_basis,
    atomic_spin_basis,
    atomic_lm_basis,
    atomic_jm_basis,
    atomic_lms_basis,
    atomic_lgs_basis,
    atomic_rank_block,
    atomic_latex_name,
    wannier90_orbital,
    cubic_orbital_basis,
    hexagonal_orbital_basis,
    cubic_basis,
    hexagonal_basis,
    orbital_latex_name,
    harmonics_key,
)
from multipiegen.data.data_compatibility_relation import compatibility_relation
from multipiegen.data.data_group_mpg import mpg_id_set, mpg_tag, mpg_id
from multipiegen.data.data_group_msg import msg_id_set, msg_tag, msg_id
from multipiegen.data.data_group_name_list import group_name_list
from multipiegen.data.data_group_pg import pg_id_set, pg_tag, pg_id
from multipiegen.data.data_group_sg import sg_id_set, sg_tag, sg_id
from multipiegen.data.data_harmonics import max_rank, max_s_rank, rank_name
from multipiegen.data.data_response_tensor import rt_cartesian_multipole, rt_multipole_list
from multipiegen.data.data_root_cluster import (
    cluster_scale,
    oh_cluster_cart_site0,
    d6h_cluster_cart_site0,
    oh_cluster_cart_vector0,
    d6h_cluster_cart_vector0,
)
from multipiegen.data.data_top_harmonics_info import top_harmonics_symmetry, top_harmonics_u_matrix

info_file = "info"

# ==================================================
h_info = """
* Information of all groups.
- "id_set" (str): (dict) info. of group's IDs.
  - "PG/SG/MPG/MSG" (str):  (dict) ID set.
    - PG: all/crystal/complex/irrep.
    - SG: all/crystal/PG.
    - MPG: all/crystal/PG/type.
    - MSG: all/crystal/PG/SG/MPG/type.
- "tag" (str): (dict) from id to tag.
- "id" (str):  (dict) from tag to id.
- "group_name" (str): (dict) name list of groups (data_group_name_list).
  - "PG/SG/MPG/MSG" (str): (dict) name info. of each group.
    - id (str): (dict) "dir" (directory name), "md" (markdown file name), "no" (number), "tag" (tag),
                       "name" (name in LaTeX), "PG/SG/MPG/MSG" (associated group id).
- "character" (str): (dict) info. of character.
  - "alias" (str): (dict) alias.
    - "wp/wm" (str): (sympy) definition of wp and wm, exp(2 pi i/3) or exp(-2 pi i/3).
  - "compatibility_relation" (str): (dict) compatibility relation.
    - (PG_tag1,PG_tag2) (str,str): ([(str,str)]) list of compatibility relation, (irrep1,irrep2).
- "root_cluster" (str): (dict) info. of root cluster.
  - "point_group" (str): ([str]) list of PG_tag for root cluster, ["Oh", "D6h"].
  - "scale" (str): (int) scale factor, scaled site = site/scale.
  - "oh_site0_cart" (str): (ndarray(3,sympy)) representative site of Oh cluster (cartesian = fractional for cubic, r=1).
  - "oh_scaled_site0_cart" (str): (ndarray(3,sympy)) scaled representative site of Oh cluster (cartesian, r=1/scale).
  - "oh_vector0_cart" (str): (ndarray(3,sympy)) representative vector of Oh cluster (cartesian = fractional for cubic, v=1).
  - "oh_site0_frac" (str): (ndarray(3,sympy)) representative site of Oh cluster (fractional for cubic, same as cartesian).
  - "oh_scaled_site0_frac" (str): (ndarray(3,sympy)) scaled representative site of Oh cluster (fractional for cubic, same as cartesian).
  - "oh_vector0_frac" (str): (ndarray(3,sympy)) representative vector of Oh cluster (fractional for cubic, same as cartesian).
  - "d6h_site0_cart" (str): (ndarray(3,sympy)) representative site of D6h cluster (cartesian, r=1).
  - "d6h_scaled_site0_cart" (str): (ndarray(3,sympy)) scaled representative site of D6h cluster (cartesian, r=1/scale).
  - "d6h_vector0_cart" (str): (ndarray(3,sympy)) representative vector of D6h cluster (cartesian, v=1).
  - "d6h_site0_frac" (str): (ndarray(3,sympy)) representative site of D6h cluster (fractional for hexagonal).
  - "d6h_scaled_site0_frac" (str): (ndarray(3,sympy)) scaled representative site of D6h cluster (fractional for hexagonal).
  - "d6h_vector0_frac" (str): (ndarray(3,sympy)) representative vector of D6h cluster (fractional for hexagonal).
- "harmonics" (str): (dict) info. of harmonics.
  - "max_rank" (str): (int) max rank of harmonics.
  - "max_s_rank" (str): (int) max internal rank.
  - "variable_site" (str): (ndarray(3,sympy)) variable of site, [x,y,z].
  - "variable_vector" (str): (ndarray(3,sympy)) variable of (bond) vector, [X,Y,Z].
  - "internal_basis" (str): (dict) internal basis of multipolar harmonics.
    - "cubic/hexagonal" (str): (dict) basis (data_basis.cubic_orbital_basis/hexagonal_orbital_basis).
      - (s,x) (int,str): ([sympy]) list of basis, x = "q" (polar)/"g" (axial).
  - "atomic_basis" (str): (dict) atomic basis.
    - "spinless" (str): (dict) spinless basis.
      - "lm/lg/cubic/hexagonal" (str): (dict) basis.
        - L (int): ([str]) list of basis, M or gamma.
    - "spinful" (str): (dict) spinful basis.
      - "jml/lms/lgs/cubic/hexagonal" (str): (dict) basis.
        - L (int): ([str]) list of basis, (J,M) or (M,s) or (gamma,s).
    - "spin" (str): ([str]) spin basis, ["u","d"].
    - "rank_block" (str): ([(int,int)]) bra-ket rank block (L1,L2) (upper triangular block only).
  - "basis_function" (str): (dict) orbital basis function.
    - name (str): ((sympy,ndarray(2L+1,sympy))) cartesian expression and u-matrix, <m|gamma>.T (m: descending).
  - "harmonics_key" (str): (dict) key of top harmonics of orbital.
    - name (str): ((str,tuple)) ("tesseral", (l,m)) or ("Oh"/"D6h", (rank, irrep, mul, comp)), m<0 for S.
  - "wannier90" (str): (dict) alias for Wannier90.
    - tag (str): (str) name corresponding to Wannier90 tag.
  - "tesseral" (str): (dict) tesseral tag.
    - name (str): ((int,int,str)) tesseral tag (L,M,"c/s").
  - "rank_name" (str): (dict) conversion between rank and name.
    - name/rank (str/int): (int/str) rank/name.
  - "latex" (str): (dict) LaTeX name.
    - "atomic" (str): (dict) name -> LaTeX of atomic orbital and spin.
    - "orbital" (str): (dict) name -> LaTeX suffix of internal basis (add prefix Q/G/T/M).
  NOTE:
    - see "multipiegen/data/data_basis.py" for detail.
- "response_tensor" (str): (dict) info. of response tensor.
  - "multipole" (str): (dict) multipole info.
    - "cubic/hexagonal" (str): ([str]) list of multipoles.
  - "cartesian_multipole" (str): (dict) linear combination info.
    - "cubic/hexagonal" (str): { comp (str): ([(int,str)]) }, relation between cartesian component and multipole.
"""

_xyz = sp.symbols("x y z", real=True)
_XYZ = sp.symbols("X Y Z", real=True)


# ==================================================
def create_group_info(info):
    info["id_set"] = {"PG": pg_id_set, "SG": sg_id_set, "MPG": mpg_id_set, "MSG": msg_id_set}
    info["tag"] = pg_tag | sg_tag | mpg_tag | msg_tag
    info["id"] = pg_id | sg_id | mpg_id | msg_id


# ==================================================
def create_group_name_info(info):
    fields = ["dir", "md", "no", "tag", "name", "PG", "SG", "MPG", "MSG"]
    info["group_name"] = {
        g_type: {id_s: dict(zip(fields, value)) for id_s, value in lst.items()} for g_type, lst in group_name_list.items()
    }


# ==================================================
def create_character_info(info):
    # complex character.
    subs = {"wp": sp.exp(2 * sp.pi * sp.I / 3), "wm": sp.exp(-2 * sp.pi * sp.I / 3)}

    # compatibility table.
    sconv = {1: "", -1: "-"}
    cconv = {-1: "", 0: "(x)", 1: "(y)", 2: "(z)"}

    ct = {}
    for tag in ["Oh", "D6h"]:
        dp = {}
        for pg, (top, dic) in compatibility_relation.items():
            if top != tag:
                continue
            dp[pg] = [sconv[vp] + ip + cconv[cp] for ((ip, cp, vp), (ia, ca, va)) in dic.values()]
        for pg1, tbl1 in dp.items():
            for pg2, tbl2 in dp.items():
                ct[(pg1, pg2)] = [(i, j) for i, j in zip(tbl1, tbl2)]

    info["character"] = {"alias": subs, "compatibility_relation": ct}


# ==================================================
def _hex_fractional(v):
    x, y, z = v
    return np.vectorize(sp.simplify)(np.array([x + y / sp.sqrt(3), 2 * y / sp.sqrt(3), z], dtype=object))


# ==================================================
def create_root_cluster_info(info):
    simplify = np.vectorize(sp.simplify)

    site_c = normalize_vector(str_to_sympy(oh_cluster_cart_site0))
    vector_c = normalize_vector(str_to_sympy(oh_cluster_cart_vector0))
    site_h = normalize_vector(str_to_sympy(d6h_cluster_cart_site0))
    vector_h = normalize_vector(str_to_sympy(d6h_cluster_cart_vector0))

    info["root_cluster"] = {
        "point_group": ["Oh", "D6h"],
        "scale": cluster_scale,
        "oh_site0_cart": site_c,
        "oh_scaled_site0_cart": simplify(site_c / cluster_scale),
        "oh_vector0_cart": vector_c,
        "oh_site0_frac": site_c.copy(),  # fractional = cartesian for cubic.
        "oh_scaled_site0_frac": simplify(site_c / cluster_scale),
        "oh_vector0_frac": vector_c.copy(),
        "d6h_site0_cart": site_h,
        "d6h_scaled_site0_cart": simplify(site_h / cluster_scale),
        "d6h_vector0_cart": vector_h,
        "d6h_site0_frac": _hex_fractional(site_h),
        "d6h_scaled_site0_frac": _hex_fractional(simplify(site_h / cluster_scale)),
        "d6h_vector0_frac": _hex_fractional(vector_h),
    }


# ==================================================
def tesseral_tag():
    """
    Tesseral tag of regular atomic orbitals, name -> (L, M, "c/s").
    The order of each rank is C(l,l), S(l,l), C(l,l-1), S(l,l-1), ..., C(l,0).
    """
    dic = {}
    for l, names in enumerate(atomic_orbital_basis):
        ms = sum([[(m, "c"), (m, "s")] for m in range(l, 0, -1)], []) + [(0, "c")]
        for name, (m, cs) in zip(names, ms):
            dic[name] = (l, m, cs)
    return dic


# ==================================================
def basis_function():
    """
    Orbital basis function, name -> (cartesian expression, <m|gamma>.T), from top harmonics.
    """
    subs = dict(zip(sp.symbols("x y z"), _xyz))
    dic = {}
    for name, (top, key) in harmonics_key.items():
        expr = str_to_sympy(top_harmonics_symmetry[top][key]).subs(subs)
        u = str_to_sympy(top_harmonics_u_matrix[top][key])
        dic[name] = (expr, u)
    return dic


# ==================================================
def internal_basis():
    """
    Internal basis of multipolar harmonics, {"cubic/hexagonal": {(s,x): [symbol]}}.
    """
    sets = {"cubic": cubic_orbital_basis, "hexagonal": hexagonal_orbital_basis}
    dic = {}
    for b_type, basis in sets.items():
        dic[b_type] = {}
        for X, x in [("Q", "q"), ("G", "g")]:
            for s in range(max_s_rank + 1):
                dic[b_type][(s, x)] = [sp.Symbol(X + orbital_latex_name[name], real=True) for name in basis[s]]
    return dic


# ==================================================
def create_harmonics_info(info):
    def by_rank(lst):
        return {l: list(names) for l, names in enumerate(lst)}

    atomic_basis = {
        "spinless": {
            "lm": by_rank(atomic_lm_basis),
            "lg": by_rank(atomic_orbital_basis),
            "cubic": by_rank(cubic_orbital_basis),
            "hexagonal": by_rank(hexagonal_orbital_basis),
        },
        "spinful": {
            "jml": by_rank(atomic_jm_basis),
            "lms": by_rank(atomic_lms_basis),
            "lgs": by_rank(atomic_lgs_basis),
            "cubic": by_rank(cubic_basis),
            "hexagonal": by_rank(hexagonal_basis),
        },
        "spin": list(atomic_spin_basis),
        "rank_block": list(atomic_rank_block),
    }

    info["harmonics"] = {
        "max_rank": max_rank,
        "max_s_rank": max_s_rank,
        "variable_site": np.array(_xyz, dtype=object),
        "variable_vector": np.array(_XYZ, dtype=object),
        "internal_basis": internal_basis(),
        "atomic_basis": atomic_basis,
        "basis_function": basis_function(),
        "harmonics_key": dict(harmonics_key),
        "wannier90": dict(wannier90_orbital),
        "tesseral": tesseral_tag(),
        "rank_name": dict(rank_name),
        "latex": {"atomic": dict(atomic_latex_name), "orbital": dict(orbital_latex_name)},
    }


# ==================================================
def create_response_tensor_info(info):
    info["response_tensor"] = {"multipole": rt_multipole_list, "cartesian_multipole": rt_cartesian_multipole}


# ==================================================
def create_info(verbose=True):
    """
    Create info (BinaryManager).

    Args:
        verbose (bool, optional): verbose ?

    Returns:
        - (BinaryManager) -- info.
    """
    info = BinaryManager(verbose=verbose, topdir=__bin_dir__)

    create_group_info(info)
    create_group_name_info(info)
    create_character_info(info)
    create_root_cluster_info(info)
    create_harmonics_info(info)
    create_response_tensor_info(info)

    info.add_comment(h_info)
    return info


# ==================================================
def save_info():
    """
    Create and save info.pkl in multipiegen/binary.
    """
    create_info().save_binary(info_file)


# ==================================================
if __name__ == "__main__":
    save_info()
