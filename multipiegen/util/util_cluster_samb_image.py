# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Cluster SAMB under the inversion coset: the reference sites of the cluster SAMB (A) and the transport matrices (C).

For every group H (real PG and SG) whose first coset in its direct parent has the inversion as its representative
(util_coset.first_coset_*, kind = "inversion"), with the image H' = I H I^-1 and the map M of the coset:

- Orbit correspondence. Every site/bond Wyckoff position W of H is paired with the Wyckoff position W' of H' onto
  which M maps its orbit: row i of W -> row perm[i] of W', the bond reversed if sign[i] = -1. Only W' with the same
  orbit size and the same number of free parameters are taken, and the lattice and centering translations and the
  reversal of a bond are allowed.
- (A) Reference. A pair with W != W' (or H != H'): if the cluster SAMB of the pair satisfy T z(label) = +- z'(label)
  already, nothing is changed; otherwise the side with the larger (group, Wyckoff order) is averaged from the site
  (bond) perm[0], the image of the representative of the other side, instead of its row 0
  (util_cluster_samb.create_cluster_samb(..., reference=...)), which gives T z = +- z' label by label.
- (C) Transport matrix. A Wyckoff position mapped onto itself (W = W' in H = H') keeps its data, and the matrix U,
  T z_a = sum_b U_ab z'_b, is recorded within each block of the same irrep (Q and T/M separately for bonds); no single
  reference gives +-1 in general (48f of SG 210).

The transport T moves the value on site (bond) i to site (bond) perm[i], with the sign of the bond reversal for T/M.
"""

import itertools

import numpy as np
import sympy as sp

from multipiegen.util import util_coset

TOL = 1e-8


# ==================================================
def _num(a):
    return np.vectorize(lambda v: complex(sp.N(v)), otypes=[complex])(np.asarray(a, dtype=object))


# ==================================================
def _mat4(s):
    m = np.array(sp.Matrix(sp.sympify(s)).evalf(), dtype=float)
    if m.shape == (3, 3):
        m4 = np.eye(4)
        m4[:3, :3] = m
        return m4
    return m


# ==================================================
_aff = {}


def _affine(arr):
    """rows of a Wyckoff position as affine functions of the parameters: C (rows, cols, nsym), d (rows, cols)."""
    key = str(arr)
    if key not in _aff:
        B = sp.Matrix(arr)
        syms = sorted(B.free_symbols, key=str)
        C = np.array([[[float(sp.diff(B[i, j], s)) for s in syms] for j in range(B.cols)] for i in range(B.rows)])
        C = C.reshape(B.rows, B.cols, len(syms))
        d = np.array([[float(B[i, j].subs({s: 0 for s in syms})) for j in range(B.cols)] for i in range(B.rows)])
        _aff[key] = (C, d)
    return _aff[key]


def _rows_at(arr):
    C, d = _affine(arr)
    x = 0.1234 + 0.0517 * np.arange(C.shape[2])
    return C @ x + d


# ==================================================
def _match(rowsA, M, arrB, bond, cent, lattice):
    """yield (perm, sign): row i of A mapped by M is row perm[i] of B (bond reversed if sign = -1)."""
    MA = []
    for r in rowsA:
        if bond:
            MA.append(np.concatenate([M[:3, :3] @ r[:3], M[:3, :3] @ r[3:] + M[:3, 3]]))
        else:
            MA.append(M[:3, :3] @ r + M[:3, 3])

    def same_pos(d):
        if not lattice:
            return np.allclose(d, 0, atol=TOL)
        return any(np.allclose(d - c, np.rint(d - c), atol=TOL) for c in cent)

    C, d = _affine(arrB)
    nB = len(d)
    shifts = [np.array(m) + c for m in itertools.product((0, 1, -1), repeat=3) for c in cent] if lattice else [np.zeros(3)]
    seen = set()
    for sgn in (1, -1):
        for n in shifts:
            for i0 in range(nB):
                t = MA[0].copy()
                if bond:
                    if sgn == -1:
                        t[:3] = -t[:3]
                    t[3:] += n
                else:
                    t = t + n
                if C.shape[2]:
                    x = np.linalg.lstsq(C[i0], t - d[i0], rcond=None)[0]
                    if not np.allclose(C[i0] @ x + d[i0], t, atol=TOL):
                        continue
                    rowsB = C @ x + d
                elif np.allclose(d[i0], t, atol=TOL):
                    rowsB = d
                else:
                    continue
                perm, sign = [], []
                for r in MA:
                    hit = None
                    for k, q in enumerate(rowsB):
                        for s in (1, -1) if bond else (1,):
                            if bond:
                                if np.allclose(q[:3] - s * r[:3], 0, atol=TOL) and same_pos(q[3:] - r[3:]):
                                    hit = (k, s)
                            elif same_pos(q - r):
                                hit = (k, s)
                            if hit:
                                break
                        if hit:
                            break
                    if hit is None:
                        break
                    perm.append(hit[0])
                    sign.append(hit[1])
                if len(perm) == len(MA) and sorted(perm) == list(range(nB)):
                    key = (tuple(perm), tuple(sign))
                    if key not in seen:
                        seen.add(key)
                        yield perm, sign


# ==================================================
def _transport(Ca, perm, sign, odd):
    """image T z of the SAMB z (dim, N): column i -> column perm[i], sign for reversed bonds of T/M."""
    Cp = np.zeros_like(Ca)
    Cp[:, perm] = Ca * (np.array(sign, dtype=float) if odd else 1.0)
    return Cp


# ==================================================
def _exact(v):
    """+- sqrt(rational) as sympy if possible, otherwise the float (complex if not real)."""
    if abs(v.imag) > TOL:
        return complex(v)
    x = v.real
    if abs(x) < TOL:
        return sp.S(0)
    q = sp.nsimplify(x * x, rational=True, tolerance=1e-12)
    if q.q < 10**6 and abs(float(q) - x * x) < 1e-12:
        return sp.sqrt(q) * (1 if x > 0 else -1)
    return float(x)


# ==================================================
def transport_matrix(A, B, perm, sign, kind):
    """
    Transport matrix of the cluster SAMB of a Wyckoff position onto those of its image.

    Args:
        A (dict): cluster SAMB of W, {PGMultipoleType: samb (dim, N)}.
        B (dict): cluster SAMB of the image W'.
        perm (list): row i of W -> row perm[i] of W'.
        sign (list): +1, or -1 for a reversed bond.
        kind (str): "site" or "bond".

    Returns:
        - (dict) -- {block: (rows of A, rows of B, U)}, None if B is not the image of A.
        - (bool) -- U is +-1 label by label (and component by component) ?

    Note:
        - block is an irrep, with " (T/M)" for T/M of a bond; rows = [(PGMultipoleType, component)]; U is a complex
          ndarray with T z[a] = sum of U[a, b] z'[b] over b.
    """
    U, diag = {}, True
    blocks = []
    for k in A:
        b = (k.Gamma, k.X in ("T", "M"))
        if b not in blocks:
            blocks.append(b)
    for blk in blocks:
        G = blk[0] + (" (T/M)" if blk[1] else "")
        ka = [k for k in A if (k.Gamma, k.X in ("T", "M")) == blk]
        kb = [k for k in B if (k.Gamma, k.X in ("T", "M")) == blk]
        ra = [(k, c) for k in ka for c in range(np.asarray(A[k]).shape[0])]
        rb = [(k, c) for k in kb for c in range(np.asarray(B[k]).shape[0])]
        if len(ra) != len(rb):
            return None, False
        Va = np.array([_transport(_num(A[k])[c : c + 1], perm, sign, kind == "bond" and k.X in ("T", "M"))[0] for k, c in ra])
        Vb = np.array([_num(B[k])[c] for k, c in rb])
        Ua = Va @ Vb.conj().T
        if not np.allclose(Ua @ Ua.conj().T, np.eye(len(ra)), atol=1e-7):
            return None, False
        if not (ka == kb and np.allclose(np.abs(Ua), np.eye(len(ra)), atol=1e-7)):
            diag = False
        U[G] = (ra, rb, Ua)
    return U, diag


# ==================================================
def inversion_pairs(group, info):
    """
    Groups whose first coset has the inversion as its representative, with its image.

    Args:
        group (dict): group.pkl.
        info (dict): info.pkl.

    Returns:
        - (list) -- [(group id, image group id, M (4x4 float, coordinates of H -> of H'), "PG"/"SG")], "-c" excluded.
    """
    tag2id = {t: i for i, t in info["tag"].items() if i.startswith("PG:")}
    pairs = []
    for gid in group:
        tp = gid.split(":")[0]
        if tp not in ("PG", "SG") or (tp == "PG" and info["tag"][gid].endswith("-c")):
            continue
        rep = util_coset.first_coset_rep(gid)
        if rep is None or rep.split(":")[0] != "-1":
            continue
        if tp == "PG":
            tag = info["tag"][gid]
            if tag.endswith("-c") or util_coset.pg_real(tag) not in util_coset.data_representative_parent_pg:
                continue
            d = util_coset.first_coset_pg(tag)
            if d["image"] is None or d["kind"] != "inversion":
                continue
            pairs.append((gid, tag2id[d["image"]], _mat4(d["map"]), tp))
        elif tp == "SG":
            no = int(gid[3:])
            if no not in util_coset.data_parent_sg:
                continue
            d = util_coset.first_coset_sg(no)
            if d["image"] is None or d["kind"] != "inversion":
                continue
            pairs.append((gid, f"SG:{d['image']}", _mat4(d["map"]), tp))
    return pairs


# ==================================================
def _pair_correspondence(group, pair):
    ga, gb, M, kg = pair
    corr = {}
    cent = [np.array(c, dtype=float) for c in group[gb]["symmetry_operation"]["plus_set"]] if kg == "SG" else [np.zeros(3)]
    for kind in ("site", "bond"):
        for wa, va in group[ga]["wyckoff"][kind].items():
            rowsA = _rows_at(va["conventional"])
            npar = len(sp.Matrix(va["conventional"]).free_symbols)
            cands = []
            for wb, vb in group[gb]["wyckoff"][kind].items():
                if len(vb["mapping"]) != len(va["mapping"]) or len(sp.Matrix(vb["conventional"]).free_symbols) != npar:
                    continue
                for perm, sign in _match(rowsA, M, vb["conventional"], kind == "bond", cent, kg == "SG"):
                    cands.append((wb, perm, sign))
            corr[(ga, kind, wa)] = (gb, cands)
    return corr


_group = {}


def _pair_worker(pair):
    if "group" not in _group:
        from multipiegen.util.util_group import load_group

        _group["group"] = load_group()
    return _pair_correspondence(_group["group"], pair)


# ==================================================
def orbit_correspondence(group, pairs, pool=None):
    """
    Candidates of the image of every Wyckoff position of the groups of the inversion pairs.

    Args:
        group (dict): group.pkl.
        pairs (list): inversion_pairs.
        pool (Pool, optional): process pool (group.pkl is loaded in each process).

    Returns:
        - (dict) -- {(group id, "site"/"bond", W): (image group id, [(W', perm, sign)])}.
    """
    corr = {}
    results = pool.map(_pair_worker, pairs) if pool is not None else [_pair_correspondence(group, p) for p in pairs]
    for c in results:
        corr.update(c)
    return corr


# ==================================================
def cluster_samb_reference(group, cs, corr):
    """
    Reference sites (bonds) of the cluster SAMB (A).

    Args:
        group (dict): group.pkl.
        cs (dict): cluster SAMB with the representative (row 0) as the reference.
        corr (dict): orbit_correspondence.

    Returns:
        - (dict) -- {group id: {"site"/"bond": {W: r}}}, average W from its site (bond) r (only r != 0).
    """
    order = {gid: n for n, gid in enumerate(group.keys())}

    def wy_order(gid, kind, w):
        return list(group[gid]["wyckoff"][kind].keys()).index(w)

    reference = {}
    for (ga, kind, wa), (gb, cands) in corr.items():
        A = cs[ga][kind][wa]
        best = None
        for wb, perm, sign in cands:
            B = cs[gb][kind][wb]
            if list(A) != list(B):
                continue
            U, diag = transport_matrix(A, B, perm, sign, kind)
            if U is None:
                continue
            if diag:
                best = ("ok", wb, perm)
                break
            if best is None:
                best = ("mix", wb, perm)
        if best is None:
            raise Exception(f"image of {ga} {kind} {wa} is not found.")
        status, wb, perm = best
        if status == "ok" or (ga, wa) == (gb, wb):
            continue
        if (order[ga], wy_order(ga, kind, wa)) < (order[gb], wy_order(gb, kind, wb)):
            reference.setdefault(gb, {}).setdefault(kind, {})[wb] = perm[0]
    return reference


# ==================================================
def cluster_samb_transport(cs, corr):
    """
    Images of the cluster SAMB under the inversion coset (after the reference change).

    Args:
        cs (dict): cluster SAMB (cluster_samb.pkl).
        corr (dict): orbit_correspondence.

    Returns:
        - (dict) -- {group id: {"site"/"bond": {W: {"image": (image group id, W'), "perm": [int], "sign": [int],
          "diag": bool, "eps": {PGMultipoleType: [sympy for each component]} (diag) or "U": {block: (rows, rows', U)} (otherwise)}}}}.
          rows = [(PGMultipoleType, component)]; U is exact (+- sqrt(rational)) where possible.
    """
    out = {}
    for (ga, kind, wa), (gb, cands) in corr.items():
        A = cs[ga][kind][wa]
        best = None
        for wb, perm, sign in cands:
            B = cs[gb][kind][wb]
            if list(A) != list(B):
                continue
            U, diag = transport_matrix(A, B, perm, sign, kind)
            if U is None:
                continue
            if best is None or (diag and not best[4]):
                best = (wb, perm, sign, U, diag)
            if diag:
                break
        if best is None:
            raise Exception(f"image of {ga} {kind} {wa} with the same labels is not found.")
        wb, perm, sign, U, diag = best
        d = {"image": (gb, wb), "perm": [int(i) for i in perm], "sign": [int(s) for s in sign], "diag": diag}
        if diag:
            d["eps"] = {}
            for ra, rb, Ua in U.values():
                for i, (k, c) in enumerate(ra):
                    d["eps"].setdefault(k, []).append(_exact(Ua[i, i]))
        else:
            d["U"] = {G: (ra, rb, np.vectorize(_exact, otypes=[object])(Ua)) for G, (ra, rb, Ua) in U.items()}
        out.setdefault(ga, {}).setdefault(kind, {})[wa] = d
    return out
