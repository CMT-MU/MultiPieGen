# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Cluster SAMB of each Wyckoff site and bond of the real point groups and the space groups.

The same data as MultiPie Ver.2's cluster_samb ("site", "bond_s", "bond_a"), created in a different way:

- General point. The cluster of the general point of a group H (as many points as the order of H, point i =
  (operation i of H) site0) is a subset of the Oh or D6h cluster: the operations of H are those of the top group
  listed in coset 0 of
  data_op_mapping.op_mapping_d6h_oh_cluster[PG tag]. The basis of the top cluster restricted to these points gives
  the basis of H: cluster_scalar_basis.pkl for "site" and "bond_s", and cluster_bond_basis.pkl (x = "q") for "bond_a"
  (on the general bond, symmetric and anti-symmetric bonds span the same space). The key of the top group is kept,
  but the irrep and the component are replaced by those of H through the top harmonics T, and the value is multiplied by the
  coefficients of data_harmonics_mapping (key of the top group = c_top T, key of H = c T: factor c c_top).
- Special point. The SAMB of a Wyckoff site (bond) is the average of the general-point SAMB, f_j = sum_{i in M_j}
  f_i, with the mapping M_j of the Wyckoff site (bond) in group.pkl (operations i mapping the representative to
  the site or bond j). For "bond_a", an operation reversing the bond (negative in the mapping) gives -f_i.
- Orthonormalization. Gram-Schmidt within each irrep of H (component by component) in the order of [Gamma, l, n, p] for
  "site" and "bond_s" (as MultiPie Ver.2), and [Gamma, l, k, Q/G, n, p] for "bond_a" (k of the source vector cluster
  basis); vanishing and linearly dependent ones are removed (selected numerically, then orthonormalized exactly).
  Each component is normalized to 1 on the N sites (bonds) of the Wyckoff position.
- Key. PGMultipoleType(X, l, Gamma, n, p, s=0, k=0, x): "site" and "bond_s" have X = "Q" and x = "q" (as
  cluster_scalar_basis); "bond_a" has the keys of cluster_bond_basis (X = Q/G, x = q); n is the multiplicity of
  the harmonics (X, l, Gamma) of H (its xyz polynomial) and p is that of the top key.
- Space group. The operations (without the plus set) of a space group have the same point parts as those of its
  point group (in the same setting), so each operation is mapped to the top group through the point group,
  op_mapping_d6h_oh_cluster[PG tag] (as "mapping" of the SG in group.pkl). (data_op_mapping.op_mapping_parent_sg
  maps to a parent space group with a change of axes, and is not used.) The harmonics are those of the point group.
- Order. Gamma (character table), l, (k, Q/G for bond_a), n, p.

- cluster_samb (dict)::

    { PG_id/SG_id: { "site"/"bond_s"/"bond_a": { wyckoff: { PGMultipoleType: samb } } } }
        samb: ndarray (dim, N) (sympy), value on the N sites (bonds) of the Wyckoff position for each component.
        The expression of the harmonics (X, l, Gamma, n) is in harmonics.pkl (util_harmonics.harmonics_expression).
"""

import numpy as np
import sympy as sp
from multipie import PGMultipoleType
from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__
from multipiegen.data.data_character_table import character
from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
from multipiegen.data.data_op_mapping import op_mapping_d6h_oh_cluster
from multipiegen.util.util_bond_cluster_basis import load_bond_cluster_basis
from multipiegen.util.util_cluster_basis import load_cluster_basis

# ==================================================
cluster_samb_file = "cluster_samb"
TOL_DEPENDENT = 1e-8
_shared = {}


# ==================================================
def _simplify(e):
    return sp.radsimp(sp.expand(e))


# ==================================================
def _dot(u, v):
    return _simplify(sum(a * b for a, b in zip(u, v)))


# ==================================================
def _top_data():
    """Top cluster bases {"site": {top: basis}, "bond": {top: (basis, source)}} (loaded once per process)."""
    if "site" not in _shared:
        cb, bb = load_cluster_basis(), load_bond_cluster_basis()
        _shared["site"] = {t: cb[t]["basis"] for t in ["Oh", "D6h"]}
        _shared["bond"] = {t: (bb[t]["basis"], bb[t]["source"]) for t in ["Oh", "D6h"]}
    return _shared


# ==================================================
def top_index(pg, pg_ops, ops=None):
    """
    Index of the top cluster point of each operation of a group (coset 0 of op_mapping_d6h_oh_cluster[pg]).

    Args:
        pg (str): point group tag.
        pg_ops (list): symmetry operation tags of the point group (group.pkl).
        ops (list, optional): symmetry operation tags of the group (space group), [pg_ops] if None.

    Returns:
        - (str) -- top group, "Oh" or "D6h".
        - (list) -- 0-based index of the top cluster for each operation of the group.
    """
    top, cosets = op_mapping_d6h_oh_cluster[pg]
    to_top = dict(zip(pg_ops, cosets[0]))
    ops = pg_ops if ops is None else [t.split(":")[0] for t in ops]
    return top, [to_top[t] - 1 for t in ops]


# ==================================================
def general_point_basis(pg, top, idx, kind):
    """
    Basis of the general point of a group (top basis restricted to coset 0, keys replaced).

    Args:
        pg (str): point group tag (of the space group for SG).
        top (str): top group, "Oh" or "D6h".
        idx (list): 0-based index of the top cluster for each operation of the group (top_index).
        kind (str): "site" (cluster_scalar_basis) or "bond" (cluster_bond_basis).

    Returns:
        - (list) -- [(order, (X, l, Gamma, n, p, x), comp, dim, ndarray (order of H,) (sympy))] sorted by order.
    """
    mapping = data_harmonics_mapping[pg]
    assert mapping["top"] == top, pg
    irreps = list(character[pg][1].keys())
    inv = {X: {tk: (key, sp.S(c)) for key, (c, tk) in mapping[X].items()} for X in ["Q", "G"]}
    top_map = data_harmonics_mapping[top]  # key of the top group -> top harmonics (axial irreps for G).

    data = _top_data()
    if kind == "site":
        basis, source = data["site"][top], {}
    else:
        basis, source = data["bond"][top]

    res = []
    for key, val in basis.items():
        if key.x != "q":  # cluster_bond_basis of polar (q) vector basis only (a complete basis by itself).
            continue
        val = np.asarray(val, dtype=object)
        dim = len(val)
        for comp, row in enumerate(val):
            c_top, tk = top_map[key.X][(key.l, key.Gamma, key.n, -1 if dim == 1 else comp)]
            (l, irrep, n, c_pg), c = inv[key.X][tk]
            c = c * sp.S(c_top)
            k = source[key].k if key in source else 0
            if kind == "site":
                order = (irreps.index(irrep), l, n, key.p)
            else:
                order = (irreps.index(irrep), l, k, ["Q", "G"].index(key.X), n, key.p)
            vec = np.array([_simplify(c * row[i]) for i in idx], dtype=object)
            res.append((order, (key.X, l, irrep, n, key.p, key.x), max(c_pg, 0), dim, vec))
    res.sort(key=lambda t: (t[0], t[2]))
    return res


# ==================================================
def wyckoff_samb(general, mapping, sign=False):
    """
    SAMB of a Wyckoff site (bond) averaged from the general-point basis.

    Args:
        general (list): general-point basis (general_point_basis).
        mapping (list): mapping of the Wyckoff site (bond), [[operation no. (1-based, negative for reversed bond)]].
        sign (bool, optional): anti-symmetric bond (reversed bond gives minus sign) ?

    Returns:
        - (dict) -- {(X, l, Gamma, n, p, x): [(comp, ndarray (N,) (sympy))]} orthonormal, in the order of general.
    """
    P = np.zeros((len(mapping), len(general[0][4])))
    for j, m in enumerate(mapping):
        for i in m:
            P[j, abs(i) - 1] += -1 if (sign and i < 0) else 1

    # numerical selection within each irrep (components of the same irrep are not orthogonal in general, e.g. the
    # real 2D irreps of complex type, E of C3).
    nums = [P @ np.asarray(v, dtype=float) for _, _, _, _, v in general]
    scale = max(np.linalg.norm(v) for v in nums)
    chains = {}
    for no, (_, key, comp, _, _) in enumerate(general):
        chains.setdefault(key[2], []).append(no)
    kept = set()
    for lst in chains.values():
        B = []
        for no in lst:
            v = nums[no].copy()
            v0 = np.linalg.norm(v)
            if v0 < TOL_DEPENDENT * scale:
                continue
            for b in B:
                v -= (b @ v) * b
            if np.linalg.norm(v) / v0 > TOL_DEPENDENT:
                B.append(v / np.linalg.norm(v))
                kept.add(no)

    # all components of a multiplet are kept or dropped together.
    status = {}
    for no, (_, key, _, _, _) in enumerate(general):
        status.setdefault(key, set()).add(no in kept)
    bad = [k for k, s in status.items() if len(s) != 1]
    if bad:
        raise Exception(f"components are partly dropped: {bad}.")

    # exact orthonormalization of the kept vectors.
    res = {}
    for lst in chains.values():
        ortho = []
        for no in lst:
            if no not in kept:
                continue
            _, key, comp, _, vec = general[no]
            u = np.array(
                [_simplify(sum(int(P[j, i]) * vec[i] for i in np.nonzero(P[j])[0])) for j in range(len(P))], dtype=object
            )
            for w, ww in ortho:
                c = _simplify(_dot(w, u) / ww)
                if c != 0:
                    u = np.array([_simplify(a - c * b) for a, b in zip(u, w)], dtype=object)
            uu = _dot(u, u)
            ortho.append((u, uu))
            nrm = sp.sqrt(uu)
            res.setdefault(key, []).append((comp, np.array([_simplify(a / nrm) for a in u], dtype=object)))
    order = {}
    for no, (o, key, _, _, _) in enumerate(general):
        order.setdefault(key, o)
    return {key: sorted(res[key], key=lambda t: t[0]) for key in sorted(res, key=lambda k: order[k])}


# ==================================================
def _to_samb(d):
    return {
        PGMultipoleType(X, l, irrep, n, p, 0, 0, x): np.array([v for _, v in comps], dtype=object)
        for (X, l, irrep, n, p, x), comps in d.items()
    }


# ==================================================
def create_cluster_samb(pg, pg_ops, ops, wyckoff):
    """
    Cluster SAMB of a group.

    Args:
        pg (str): point group tag (of the space group for SG).
        pg_ops (list): symmetry operation tags of the point group (group.pkl).
        ops (list): symmetry operation tags of the group (group.pkl).
        wyckoff (dict): group["wyckoff"] of group.pkl.

    Returns:
        - (dict) -- {"site"/"bond_s"/"bond_a": {wyckoff: {PGMultipoleType: samb}}}.
    """
    top, idx = top_index(pg, pg_ops, ops)
    site = general_point_basis(pg, top, idx, "site")
    bond = general_point_basis(pg, top, idx, "bond")

    res = {"site": {}, "bond_s": {}, "bond_a": {}}
    for w, val in wyckoff["site"].items():
        res["site"][w] = _to_samb(wyckoff_samb(site, val["mapping"]))
    for w, val in wyckoff["bond"].items():
        res["bond_s"][w] = _to_samb(wyckoff_samb(site, val["mapping"]))
        res["bond_a"][w] = _to_samb(wyckoff_samb(bond, val["mapping"], sign=True))
    return res


# ==================================================
def load_cluster_samb():
    """
    Load cluster SAMB (multipiegen/binary/cluster_samb.pkl).

    Returns:
        - (BinaryManager) -- {PG_id/SG_id: {"site"/"bond_s"/"bond_a": {wyckoff: {PGMultipoleType: samb}}}}.
    """
    return BinaryManager(cluster_samb_file, topdir=__bin_dir__)
