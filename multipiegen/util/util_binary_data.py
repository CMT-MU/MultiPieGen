# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Binary data for MultiPie (multipiegen/binary_data), in the layout of MultiPie's multipie/binary_data.

- info.pkl: information of all groups (= binary/info.pkl).
- {PG,SG,MPG,MSG}/{id}.pkl: data of each group, with ":" of the id replaced by "_" (git on Windows cannot check out a
  file name with ":"), e.g. PG/PG_32.pkl, SG/SG_221.pkl, MPG/MPG_32.1.118.pkl, MSG/MSG_221.92.pkl. Each file has its
  own header (str() of the loaded data).

The data are taken from the binaries in multipiegen/binary (nothing is recomputed, except the representation
matrices of PG by util_representation_matrix), and regrouped by group:

- PG: info, symmetry_operation, wyckoff, character (group.pkl), representation_matrix, harmonics (harmonics.pkl),
  multipolar_harmonics (multipolar_harmonics.pkl), atomic_samb (atomic_multipole_group.pkl), cluster_samb
  (cluster_samb.pkl, real PG only), and coset ("first" of coset.pkl, real PG with a parent only).
- SG: info, symmetry_operation, wyckoff (group.pkl), cluster_samb (cluster_samb.pkl), and coset ("first" of
  coset.pkl, SG with a parent only).

info.pkl has also "parent" and "coset" (the parents and the cosets in them, from coset.pkl).
- MPG: info, symmetry_operation, wyckoff, active_multipole (group.pkl), and response_tensor (response_tensor.pkl).
- MSG: info, symmetry_operation, wyckoff (group.pkl).

The SAMB and the multipolar harmonics hold no expression of the harmonics: look it up by (X, l, Gamma, n) of the key
in "harmonics" of the PG file (T -> Q, M -> G); for SG and MPG, the PG file is that of info.PG.
"""

import os

from multipie.util.util_binary import BinaryManager

from multipiegen import __top_dir__
from multipiegen.util.util_atomic_multipole import load_atomic_multipole_group
from multipiegen.util.util_cluster_samb import load_cluster_samb
from multipiegen.util.util_coset import load_coset
from multipiegen.util.util_group import load_group
from multipiegen.util.util_harmonics import load_harmonics
from multipiegen.util.util_harmonics_multipole import load_multipolar_harmonics
from multipiegen.util.util_info import load_info
from multipiegen.util.util_representation_matrix import representation_matrix
from multipiegen.util.util_response_tensor import load_response_tensor

# ==================================================
binary_data_dir = os.path.join(__top_dir__, "multipiegen", "binary_data") + os.sep
group_types = ["PG", "SG", "MPG", "MSG"]

h_info_layout = """
* Layout of the binary data (multipiegen/binary_data, as MultiPie's multipie/binary_data).
- info.pkl: this file.
- {PG,SG,MPG,MSG}/{id}.pkl: data of each group, with ":" of the id replaced by "_", e.g. "PG/PG_32" (id "PG:32"),
  "SG/SG_221", "MPG/MPG_32.1.118", "MSG/MSG_221.92".
  str() of each file shows the description of its keys.
"""

h_representation_matrix = """
- "representation_matrix" (str): (dict) representation matrices of the irreps.
  - "tag" (str): ([str]) list of SO tag.
  - "matrix" (str): (dict) representation matrices.
    - irrep_tag (str): (ndarray(n,dim,dim,sympy)) matrix of each SO, in the order of "tag".
  NOTE:
    - obtained from the u-matrices of the harmonics (util_representation_matrix).
"""

h_coset = """
- "coset" (str): (dict) first coset g H of the group H in its direct parent (the first one other than H).
  (only for a group with a parent; not for "-c" PG.)
  - "parent" (str): (str) direct parent.
  - "R" (str): (Matrix(sympy)) x_parent = R x_H.
  - "index" (str): (int) number of cosets [P:H].
  - "ops" (str): ([str]) SO tags of the parent in the coset.
  - "rep" (str): (str) SO tag of the representative g (first operation).
  - "rep_matrix" (str): (Matrix(sympy)) g in the coordinates of H (4x4 for SG, 3x3 for PG).
  - "order" (str): (int) smallest n with g^n in H (modulo the lattice of H).
  - "kind" (str): (str) "rotation" or "inversion".
  - "frame" (str): (Matrix(sympy)) frame of the image system (g in the coordinates of H) for "rotation", None otherwise.
  - "inversion" (str): (dict) for "inversion", None otherwise.
    - "image" (str): (str) image group H' = g H g^-1.
    - "map" (str): (Matrix(sympy)) M: coordinates of H -> those of H'.
    - "shift" (str): (Matrix(sympy)) origin shift k (coordinates of the parent) with k H' k^-1 = g H g^-1 (SG only).
    - "wyckoff" (str): (dict) W: (W', {parameter of W': expression of the parameters of W}).
    - "site"/"bond" (str): (dict) image of each site/bond Wyckoff position.
      - W (str): (dict) "image": (H', W'), "perm": row i of W -> row perm[i] of W', "sign": -1 for a reversed bond,
        "diag": the cluster SAMB are mapped label by label, "eps": {PGMultipoleType: [+-1 for each component]} (diag)
        or "U": {block: (rows of W, rows of W', U)} (otherwise), rows = [(PGMultipoleType, component)],
        T z_a = sum_b U_ab z'_b within each block (irrep, Q and T/M separately for bond).
  - "parent_limit" (str): (dict) parent structures, W: [{"param": {parameter: value}, "position": representative
    position, "parent_wyckoff": Wyckoff position of the parent, "preserved_by": indices of the cosets (in the order of
    the coset images) that keep it}].
"""

h_info_coset = """
* Parents and cosets of the real PG and SG (from coset.pkl).
- "parent" (str): (dict) parents of each group.
  - PG_id/SG_id (str): ([str]) parents along the chain of representative parents, the direct parent first.
- "coset" (str): (dict) left cosets g H in each parent.
  - PG_id/SG_id (str): (dict)
    - parent_id (str): ([[str]]) operations of each coset, its representative first, without H itself and the cosets
      of translations only (PG: SO tag of the parent; SG: "point part of the SO tag of the parent:[translation in the
      coordinates of H modulo its lattice]").
"""

_section_title = {
    "harmonics": "(dict) harmonics (rank l <= 11): cartesian expression and u-matrix for each multipole key.",
    "multipolar_harmonics": "(dict) multipolar harmonics (s = 1..3, rank l <= 3) with polar (q) and axial (g) internal basis.",
    "atomic_samb": "(dict) matrix of atomic multipole SAMB (rank l <= 7, s = 0, 1, orbitals of L = 0..3).",
    "cluster_samb": "(dict) cluster SAMB of each Wyckoff site and bond.",
    "response_tensor": "(dict) response tensors (up to rank 4).",
}

_shared = {}


# ==================================================
def _sources():
    """Source binaries in multipiegen/binary (loaded once per process)."""
    if not _shared:
        _shared.update(
            {
                "info": load_info(),
                "group": load_group(),
                "harmonics": load_harmonics(),
                "multipolar_harmonics": load_multipolar_harmonics(),
                "atomic_samb": load_atomic_multipole_group(),
                "cluster_samb": load_cluster_samb(),
                "coset": load_coset(),
                "response_tensor": load_response_tensor(),
            }
        )
    return _shared


# ==================================================
def _group_section(header, tp):
    """
    Description of the group data of a type in the header of group.pkl, one level up.

    Args:
        header (str): header of group.pkl.
        tp (str): "PG", "SG", "MPG" or "MSG".

    Returns:
        - (str) -- description of the keys "info", "symmetry_operation", ... of a group file.
    """
    lines = header.splitlines()
    i = next(n for n, l in enumerate(lines) if l.startswith(f'* Information for all {tp} ("{tp}")'))
    j = next((n for n in range(i + 1, len(lines)) if lines[n].startswith("* ")), len(lines))
    body = lines[i + 2 : j]  # skip the title and "- {tp}_id (str): ..." lines.
    return "\n".join(l[2:] for l in body).strip("\n")


# ==================================================
def _source_section(header, name):
    """
    Description of a source binary as a key of a group file.

    The first "- PG_id (str): (dict) text" line becomes '- "name" (str): (dict) text', the lines below it are kept, and
    the top-level NOTE is moved under the key. References to the source binaries are replaced by the group files.

    Args:
        header (str): header of the source binary.
        name (str): key in the group file.

    Returns:
        - (str) -- description.
    """
    lines = header.strip("\n").splitlines()
    i = next(n for n, l in enumerate(lines) if l.startswith("- "))
    first = f'- "{name}" (str): {_section_title[name]}'
    out, in_note = [first], False
    for l in lines[i + 1 :]:
        if l.startswith("NOTE"):
            in_note = True
            out.append("  " + l)
        else:
            out.append("  " + l if in_note else l)
    s = "\n".join(out)
    s = s.replace("of a key of the other binaries", "of a key of the other data")
    s = s.replace("is in harmonics.pkl", 'is in "harmonics" of the PG file (of info.PG for SG and MPG)')
    s = s.replace("cluster_scalar_basis (", "the scalar cluster basis (").replace(
        "cluster_bond_basis (", "the bond cluster basis ("
    )
    return s


# ==================================================
def group_header(tp):
    """
    Header of a group file.

    Args:
        tp (str): "PG", "SG", "MPG" or "MSG".

    Returns:
        - (str) -- header.
    """
    src = _sources()
    title = {
        "PG": "* Data of a PG.",
        "SG": "* Data of an SG.",
        "MPG": "* Data of an MPG.",
        "MSG": "* Data of an MSG.",
    }[tp]
    parts = [title, _group_section(src["group"].comment, tp)]
    if tp == "PG":
        parts += [
            h_representation_matrix.strip("\n"),
            _source_section(src["harmonics"].comment, "harmonics"),
            _source_section(src["multipolar_harmonics"].comment, "multipolar_harmonics"),
            _source_section(src["atomic_samb"].comment, "atomic_samb"),
            _source_section(src["cluster_samb"].comment, "cluster_samb") + '\n    - only for real PG (not for "-c").',
            h_coset.strip("\n"),
        ]
    elif tp == "SG":
        parts.append(_source_section(src["cluster_samb"].comment, "cluster_samb"))
        parts.append(h_coset.strip("\n"))
    elif tp == "MPG":
        parts.append(_source_section(src["response_tensor"].comment, "response_tensor"))
    parts.append("* all sympy values are expanded with rationalized denominators.")
    return "\n".join(parts) + "\n"


# ==================================================
def group_file(gid):
    """
    File name of a group in binary_data.

    Args:
        gid (str): group id, e.g. "PG:32".

    Returns:
        - (str) -- file name, e.g. "PG/PG_32".

    Note:
        - ":" of the id is replaced by "_", as a file name with ":" cannot be checked out by git on Windows.
    """
    return f"{gid.split(':')[0]}/{gid.replace(':', '_')}"


# ==================================================
def create_group_data(gid):
    """
    Data of a group file.

    Args:
        gid (str): group id, e.g. "PG:32".

    Returns:
        - (BinaryManager) -- group data with its header.
    """
    src = _sources()
    tp = gid.split(":")[0]
    bm = BinaryManager(topdir=binary_data_dir, comment=group_header(tp))
    bm.update(src["group"][gid])
    if tp == "PG":
        bm["representation_matrix"] = representation_matrix(src["info"]["tag"][gid])
        bm["harmonics"] = src["harmonics"][gid]
        bm["multipolar_harmonics"] = src["multipolar_harmonics"][gid]
        bm["atomic_samb"] = src["atomic_samb"][gid]
        if gid in src["cluster_samb"]:
            bm["cluster_samb"] = src["cluster_samb"][gid]
        if gid in src["coset"]:
            bm["coset"] = src["coset"][gid]["first"]
    elif tp == "SG":
        bm["cluster_samb"] = src["cluster_samb"][gid]
        if gid in src["coset"]:
            bm["coset"] = src["coset"][gid]["first"]
    elif tp == "MPG":
        bm["response_tensor"] = src["response_tensor"][gid]
    return bm


# ==================================================
def create_info_data():
    """
    Data of info.pkl.

    Returns:
        - (BinaryManager) -- info with its header.
    """
    info, coset = _sources()["info"], _sources()["coset"]
    bm = BinaryManager(topdir=binary_data_dir, comment=info.comment)
    bm.add_comment(h_info_coset)
    bm.add_comment(h_info_layout)
    bm.update(info)
    bm["parent"] = {gid: d["parent"] for gid, d in coset.items()}
    bm["coset"] = {gid: d["coset"] for gid, d in coset.items()}
    return bm


# ==================================================
def group_ids():
    """
    Ids of all groups, in the order PG, SG, MPG, MSG.

    Returns:
        - (list) -- group ids.
    """
    info = _sources()["info"]
    return [gid for tp in group_types for gid in info["tag"] if gid.split(":")[0] == tp]


# ==================================================
def load_binary_data(name):
    """
    Load a file of binary_data.

    Args:
        name (str): "info" or a group id, e.g. "PG:32".

    Returns:
        - (BinaryManager) -- data, str() gives the description of the keys.
    """
    return BinaryManager(name if name == "info" else group_file(name), topdir=binary_data_dir)
