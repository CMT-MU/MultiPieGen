# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Direct supergroups of the real PG and SG: the supergroups P in which a group H is a maximal subgroup, with the
embeddings of H in P, from the symmetry operations of group.pkl.

- Embedding: the affine map x_P = A x_H (A = [[R, s], [0, 1]], conventional coordinates) that maps every operation of H
  to one of P. R is right-handed (an enantiomorphic pair is not a sub/supergroup). For SG, R maps the lattice of H
  into that of P: N = L_P^-1 R L_H is an integer matrix (L: primitive vectors in the conventional coordinates), and
  [T_P : T_H] = det N. For PG, R is orthogonal in the Cartesian coordinates.

  - R: from the images V_i of the generators W_i of the point group of H (N W'_i = V'_i N in the primitive bases, a
    linear subspace of N), its integer points with small entries (abs(N_ij) <= 3 for a subspace of dimension <= 3, <= 2 for
    <= 5, and the signed permutations otherwise, i.e., triclinic H).
  - s: on the grid of 1/24 in the primitive coordinates of P.
  - Two embeddings are the same when they give the same set of the operations of P in the coordinates of H modulo
    the lattice of P (they differ by an element of the normalizer of P); different embeddings give different parent
    structures for the same atoms of H.

- Direct (maximal): by the theorem of Hermann, a maximal subgroup is a t-subgroup ([T_P : T_H] = 1, point-group index
  2, 3 or 4) or a k-subgroup (the same point-group order, [T_P : T_H] = 2, 3 or 4 here). An embedding is direct if no
  other embedding Q (of these candidates, isomorphic k-subgroups included) lies between: H < Q < P as the sets of
  the elements modulo the lattice of H. The isomorphic supergroups (P of the same type as H) are not listed.
- Index: [P : H] modulo the lattice of H ("index") and modulo the lattice of P ("index_parent").

NOTE: the "direct parent" of data_parent_pg/sg (coset.pkl) is the representative parent of MultiPieGen and is not
always a maximal supergroup (e.g., SG 152 < 150 < 164 with the lattice of 164 (c/3)).
"""

import functools
import itertools

import numpy as np
import sympy as sp
from multipie.util.util_crystal import P_dict

from multipiegen.util.util_group import load_group
from multipiegen.util.util_info import load_info

TOL = 1e-6
GRID = 24
MAX_LATTICE_INDEX = 4
_HEX = np.array([[1, -0.5, 0], [0, np.sqrt(3) / 2, 0], [0, 0, 1]])


# ==================================================
@functools.lru_cache(maxsize=1)
def _groups():
    return load_group()


# ==================================================
def _float(M):
    return np.array(sp.Matrix(M).evalf(), dtype=float)


# ==================================================
def _mat4(W, t=None):
    A = np.eye(4)
    A[:3, :3] = W
    if t is not None:
        A[:3, 3] = t
    return A


# ==================================================
def _wkey(W):
    return tuple(np.rint(np.asarray(W, dtype=float) * 1e6).astype(np.int64).ravel())


# ==================================================
def _reduce(v, L):
    """v modulo the lattice L (None: no lattice)."""
    v = np.asarray(v, dtype=float)
    if L is None:
        return np.where(np.abs(v) < 1e-9, 0.0, v)
    f = np.linalg.solve(L, v)
    f = f - np.floor(f + TOL)
    f = np.where(np.abs(f) < 1e-9, 0.0, f)
    return L @ f


# ==================================================
def okey(W, L):
    """key of an operation (4x4) modulo the lattice L."""
    return _wkey(W[:3, :3]), tuple(np.rint(_reduce(W[:3, 3], L) * 1e5).astype(np.int64))


# ==================================================
def _order(W):
    n, A = 1, np.asarray(W, dtype=float)
    while not np.allclose(A, np.eye(3)):
        A = A @ W
        n += 1
    return n


# ==================================================
@functools.lru_cache(maxsize=None)
def group_data(gid):
    """point operations, the operations (with the plus set), lattice (primitive vectors), of a group."""
    g = _groups()[gid]
    so = g["symmetry_operation"]
    pg = gid.startswith("PG")
    ops = [_float(W) for W in so["fractional"]]
    if pg:
        ops = [_mat4(W) for W in ops]
    L = None if pg else np.array(_float(P_dict[g["info"].lattice]))[:3, :3]
    plus = [np.zeros(3)] if pg else [np.array(c, dtype=float) for c in _float(so["plus_set"])]
    full = []
    for c in plus:
        for W in ops:
            W = W.copy()
            W[:3, 3] = W[:3, 3] + c
            full.append(W)
    point = {}
    for W, tag in zip(ops, so["tag"]):
        point.setdefault(_wkey(W[:3, :3]), (W[:3, :3], W[:3, 3], str(tag).split(":")[0]))
    return {"id": gid, "L": L, "full": full, "plus": plus, "point": point, "hex": bool(g["info"].hexagonal_g)}


# ==================================================
def _generators(mats):
    """a small generating set of a finite matrix group (list of 3x3)."""
    keys = {_wkey(W): W for W in mats}
    gens, span = [], {_wkey(np.eye(3)): np.eye(3)}
    for W in sorted(keys.values(), key=lambda W: (-_order(W), _wkey(W))):
        if _wkey(W) in span:
            continue
        gens.append(W)
        frontier = list(span.values())
        while frontier:
            new = []
            for A in frontier:
                for G in gens:
                    B = A @ G
                    k = _wkey(B)
                    if k not in span:
                        span[k] = B
                        new.append(B)
            frontier = new
        if len(span) == len(keys):
            break
    return gens


# ==================================================
def _nullspace(K):
    """rational basis (rows, RREF) and pivots of the null space of the integer matrix K."""
    ns = sp.Matrix(np.rint(K).astype(int)).nullspace()
    if not ns:
        return None, []
    B = sp.Matrix.hstack(*ns).T.rref()[0]
    piv = [next(j for j in range(B.cols) if B[i, j] != 0) for i in range(B.rows)]
    return B, piv


# ==================================================
_SIGNED_PERM = [
    np.array(P, dtype=float) * np.array(s, dtype=float)[:, None]
    for P in itertools.permutations(np.eye(3).tolist())
    for s in itertools.product([1, -1], repeat=3)
]


# ==================================================
def _point_candidates(Hd, Pd):
    """linear parts N (primitive bases) mapping the point group of H into that of P (right-handed)."""
    LH = np.eye(3) if Hd["L"] is None else Hd["L"]
    LP = np.eye(3) if Pd["L"] is None else Pd["L"]
    Hp = [np.linalg.inv(LH) @ W @ LH for W, _, _ in Hd["point"].values()]
    Pp = [np.linalg.inv(LP) @ V @ LP for V, _, _ in Pd["point"].values()]
    Pkeys = {_wkey(V) for V in Pp}
    gens = _generators(Hp)
    sig = lambda W: (round(np.linalg.det(W)), round(np.trace(W)), _order(W))
    cands = [[V for V in Pp if sig(V) == sig(W)] for W in gens]
    I3 = np.eye(3)
    sols = {}

    def rec(i, K):
        if K is not None and np.linalg.matrix_rank(K, tol=1e-8) >= 9:
            return
        if i == len(gens):
            B, piv = _nullspace(K)
            if B is not None:
                sols.setdefault(tuple(np.array(B.evalf(), dtype=float).ravel().round(8)), (B, piv))
            return
        for V in cands[i]:
            Ki = np.kron(I3, gens[i].T) - np.kron(V, I3)  # vec(N W - V N) = Ki vec(N) (row major).
            rec(i + 1, Ki if K is None else np.vstack([K, Ki]))

    pg = Hd["L"] is None
    Ns = []
    if not gens:  # point group {1}: signed permutations.
        Ns = list(_SIGNED_PERM)
    else:
        rec(0, None)
        for B, piv in sols.values():
            d = len(piv)
            Bf = np.array(B.evalf(), dtype=float)
            if d > 5:
                Ns += [N for N in _SIGNED_PERM if np.allclose(np.linalg.lstsq(Bf.T, N.ravel(), rcond=None)[0] @ Bf, N.ravel())]
                continue
            if pg:
                vals = [-1, -0.5, 0, 0.5, 1]
            else:
                vals = list(range(-3, 4)) if d <= 3 else list(range(-2, 3))
            for c in itertools.product(vals, repeat=d):
                Ns.append((np.array(c, dtype=float) @ Bf).reshape(3, 3))
    out = {}
    for N in Ns:
        if np.linalg.det(N) < 1e-6:
            continue
        if not pg and not np.allclose(N, np.rint(N), atol=1e-8):
            continue
        Ni = np.linalg.inv(N)
        if not all(_wkey(N @ W @ Ni) in Pkeys for W in Hp):
            continue
        if pg:
            CH = _HEX if Hd["hex"] else np.eye(3)
            CP = _HEX if Pd["hex"] else np.eye(3)
            O = CP @ N @ np.linalg.inv(CH)
            if not np.allclose(O @ O.T, np.eye(3), atol=1e-6):
                continue
        out.setdefault(_wkey(N), N)
    return list(out.values())


# ==================================================
_grid_cache = {}


def _grid():
    if GRID not in _grid_cache:
        r = np.arange(GRID) / GRID
        _grid_cache[GRID] = np.array(list(itertools.product(r, r, r)))
    return _grid_cache[GRID]


# ==================================================
def _shifts(Hd, Pd, R):
    """shifts s (conventional of P) with x_P = R x_H + s mapping H into P."""
    LP = Pd["L"]
    LPi = np.linalg.inv(LP)
    S = _grid()
    ok = np.ones(len(S), dtype=bool)
    Ri = np.linalg.inv(R)
    for W in Hd["full"]:
        V = R @ W[:3, :3] @ Ri
        t = Pd["point"].get(_wkey(V))
        if t is None:
            return []
        a = LPi @ (R @ W[:3, 3] - t[1])
        M = LPi @ (np.eye(3) - V) @ LP
        idx = np.nonzero(ok)[0]
        res = a[None, :] + S[idx] @ M.T
        good = np.all(np.abs(res - np.rint(res)) < 1e-6, axis=1)
        ok[idx[~good]] = False
        if not ok.any():
            return []
    return [LP @ s for s in S[ok]]


# ==================================================
@functools.lru_cache(maxsize=None)
def _normalizer_translations(h):
    """translations t (primitive coordinates of H, grid of 1/24) normalizing H: (I - W) t in the lattice for every W."""
    Hd = group_data(h)
    LH = Hd["L"]
    S = _grid()
    ok = np.ones(len(S), dtype=bool)
    for W, _, _ in Hd["point"].values():
        M = np.eye(3) - np.linalg.inv(LH) @ W @ LH
        res = S @ M.T
        ok &= np.all(np.abs(res - np.rint(res)) < 1e-6, axis=1)
    return S[ok]


# ==================================================
@functools.lru_cache(maxsize=None)
def _point_normalizer(h):
    """linear parts M (primitive coordinates of H) of the affine normalizer of H (within the entries searched)."""
    Hd = group_data(h)
    out = []
    LH = np.eye(3) if Hd["L"] is None else Hd["L"]
    for M in _point_candidates(Hd, Hd):
        if round(abs(np.linalg.det(M))) != 1:
            continue
        if Hd["L"] is not None and not _shifts(Hd, Hd, LH @ M @ np.linalg.inv(LH)):
            continue
        out.append(M)
    return out


# ==================================================
def _shift_orbits(Hd, Pd, R, shifts):
    """one shift for each orbit of the translations normalizing H (s ~ s + R t)."""
    T = _normalizer_translations(Hd["id"])
    LP = Pd["L"]
    D = (np.linalg.inv(LP) @ R @ Hd["L"] @ T.T).T  # primitive coordinates of P.
    D = D[np.all(np.abs(D * GRID - np.rint(D * GRID)) < 1e-6, axis=1)]
    Di = np.rint(D * GRID).astype(np.int64) % GRID
    idx = lambda v: (v[..., 0] * GRID + v[..., 1]) * GRID + v[..., 2]
    seen = np.zeros(GRID**3, dtype=bool)
    out = []
    for s in shifts:
        si = np.rint(np.linalg.solve(LP, s) * GRID).astype(np.int64) % GRID
        k = idx(si)
        if seen[k]:
            continue
        out.append(s)
        seen[idx((si[None, :] + Di) % GRID)] = True
    return out


# ==================================================
def _maximal(Hd, Pd, A):
    """is H a maximal subgroup of P (embedding A) ? (no g in P - H with H < <H, g> < P, modulo the lattice of H)."""
    L = Hd["L"]
    ops, LPh = _parent_ops(Pd, A)
    trans = [np.zeros(3)] if LPh is None else [L @ t for t in _lattice_points(Hd, LPh)]
    E = {}
    for W in ops:
        for t in trans:
            M = _mat4(W[:3, :3], _reduce(W[:3, 3] + t, L))
            E.setdefault(okey(M, L), M)
    Hk = {okey(W, L): W for W in Hd["full"]}
    if L is not None:
        Hk = {okey(_mat4(W[:3, :3], _reduce(W[:3, 3], L)), L): W for W in Hd["full"]}
    hgens = _sg_generators(list(Hk.values()), L)
    tested = set()
    for k, g in E.items():
        if k in Hk or k in tested:
            continue
        span = dict(Hk)
        span[k] = g
        frontier = [g]
        gens = hgens + [g]
        while frontier and len(span) < len(E):
            new = []
            for X in frontier:
                for G in gens:
                    for Y in (X @ G, G @ X):
                        Y[:3, 3] = _reduce(Y[:3, 3], L)
                        kk = okey(Y, L)
                        if kk not in span:
                            span[kk] = Y
                            new.append(Y)
            frontier = new
        if len(span) < len(E):
            return False
        tested.add(k)
    return True


# ==================================================
def _sg_generators(ops, L):
    """a generating set of a finite group of operations modulo L."""
    keys = {okey(W, L): W for W in ops}
    gens, span = [], {okey(np.eye(4), L)}
    for k, W in sorted(keys.items(), key=lambda kv: (-_order(kv[1][:3, :3]), kv[0])):
        if k in span:
            continue
        gens.append(W)
        span = {okey(np.eye(4), L)}
        frontier = [np.eye(4)]
        while frontier:
            new = []
            for X in frontier:
                for G in gens:
                    Y = X @ G
                    Y[:3, 3] = _reduce(Y[:3, 3], L)
                    kk = okey(Y, L)
                    if kk not in span:
                        span.add(kk)
                        new.append(Y)
            frontier = new
        if len(span) == len(keys):
            break
    return gens


# ==================================================
def _lattice_points(Hd, L):
    """points of the lattice L modulo the lattice of H (coordinates: primitive of H)."""
    gens = [np.linalg.solve(Hd["L"], L[:, i]) for i in range(3)]
    pts = {tuple(np.zeros(3)): np.zeros(3)}
    frontier = [np.zeros(3)]
    while frontier:
        new = []
        for t in frontier:
            for g in gens:
                u = t + g
                u = u - np.floor(u + TOL)
                k = tuple(np.round(u, 5) + 0.0)
                if k not in pts:
                    pts[k] = u
                    new.append(u)
        frontier = new
    return list(pts.values())


# ==================================================
def _parent_ops(Pd, A):
    Ai = np.linalg.inv(A)
    return [Ai @ W @ A for W in Pd["full"]], (None if Pd["L"] is None else Ai[:3, :3] @ Pd["L"])


# ==================================================
def class_key(Hd, Pd, A):
    """the operations of P in the coordinates of H modulo the lattice of P, and that lattice modulo that of H."""
    ops, L = _parent_ops(Pd, A)
    if L is None:
        return None, frozenset(okey(W, None) for W in ops)
    pts = _lattice_points(Hd, L)
    keys = set()
    for W in ops:
        c = np.linalg.solve(Hd["L"], W[:3, 3])
        cand = []
        for t in pts:
            u = c + t
            u = u - np.floor(u + TOL)
            cand.append(tuple(np.rint(u * 1e5).astype(np.int64)))
        keys.add((_wkey(W[:3, :3]), min(cand)))
    return frozenset(tuple(np.round(t, 5) + 0.0) for t in pts), frozenset(keys)


# ==================================================
def element_keys(Hd, Pd, A):
    """the elements of P modulo the lattice of H (coordinates of H)."""
    ops, L = _parent_ops(Pd, A)
    if L is None:
        return frozenset(okey(W, None) for W in ops)
    trans = [Hd["L"] @ t for t in _lattice_points(Hd, L)]
    return frozenset(okey(_mat4(W[:3, :3], W[:3, 3] + t), Hd["L"]) for W in ops for t in trans)


# ==================================================
def _cost(A):
    R, s = A[:3, :3], A[:3, 3]
    return (
        int(np.sum(np.diag(R) < 0)),
        float(np.sum(np.abs(R - np.diag(np.diag(R))))),
        float(np.sum(np.abs(R))),
        float(np.sum(np.minimum(s % 1, 1 - s % 1))),
        tuple(np.round(s % 1, 6)),
    )


# ==================================================
def embeddings(h, p, lattice_index=None):
    """
    Embeddings of H in P.

    Args:
        h (str): id of H.
        p (str): id of P.
        lattice_index (list, optional): allowed [T_P : T_H] (default: 1-4).

    Returns:
        - (list) -- [(A (4x4 float), index, index_parent)], the simplest first.
    """
    Hd, Pd = group_data(h), group_data(p)
    pg = Hd["L"] is None
    if (Pd["L"] is None) != pg:
        return []
    nH, nP = len(Hd["point"]), len(Pd["point"])
    if nP % nH:
        return []
    allowed = range(1, MAX_LATTICE_INDEX + 1) if lattice_index is None else lattice_index
    LH = np.eye(3) if pg else Hd["L"]
    LP = np.eye(3) if pg else Pd["L"]
    classes = {}
    Ms = [] if p == h else _point_normalizer(h)
    seen = set()
    cands = sorted(_point_candidates(Hd, Pd), key=lambda N: _cost(_mat4(LP @ N @ np.linalg.inv(LH))))
    for N in cands:
        if _wkey(N) in seen:
            continue
        for M in Ms:
            seen.add(_wkey(N @ M))
        lat = round(abs(np.linalg.det(N)))
        if not pg and lat not in allowed:
            continue
        R = LP @ N @ np.linalg.inv(LH)
        R = np.where(np.abs(R - np.rint(R * 6) / 6) < 1e-9, np.rint(R * 6) / 6, R) + 0.0
        if pg:
            shifts = [np.zeros(3)]
        else:
            if not all(
                np.allclose(np.linalg.solve(LP, R @ c), np.rint(np.linalg.solve(LP, R @ c)), atol=1e-6) for c in Hd["plus"]
            ):
                continue
            shifts = _shifts(Hd, Pd, R)
        if not pg:
            shifts = _shift_orbits(Hd, Pd, R, shifts)
        for s in shifts:
            A = _mat4(R, _reduce(s, Pd["L"]))
            k = class_key(Hd, Pd, A)
            if k not in classes or _cost(A) < _cost(classes[k][0]):
                classes[k] = (A, nP // nH * (1 if pg else lat), nP // nH)
    return sorted(classes.values(), key=lambda v: (v[1], _cost(v[0])))


# ==================================================
def direct_supergroups(h):
    """
    Direct supergroups of a group H.

    Args:
        h (str): id of H (real PG or SG).

    Returns:
        - (list) -- [{"parent": id, "A": 4x4 float (x_P = A x_H), "index", "index_parent"}], in order of the index and
          the group number.
    """
    info = load_info()
    tp = h.split(":")[0]
    Hd = group_data(h)
    nH = len(Hd["point"])
    pg = tp == "PG"
    nodes = []
    for p in info["id_set"][tp]["all"]:
        if p == h:
            continue  # isomorphic.
        nP = len(group_data(p)["point"])
        if nP % nH:
            continue
        r = nP // nH
        if r in (2, 3, 4):  # t-subgroup.
            lat = [1]
        elif r == 1 and not pg:  # k-subgroup.
            lat = [2, 3, 4]
        else:
            continue
        for A, n, n_p in embeddings(h, p, lat):
            nodes.append((p, A, n, n_p))
    out = []
    for p, A, n, n_p in nodes:
        if _maximal(Hd, group_data(p), A):
            out.append({"parent": p, "A": A, "index": n, "index_parent": n_p})
    out.sort(key=lambda d: (d["index"], int(d["parent"].split(":")[1]), _cost(d["A"])))
    return out


# ==================================================
supergroup_file = "supergroup"


def exact_matrix(A):
    """float 4x4 -> sympy Matrix of rationals."""
    return sp.Matrix(4, 4, lambda i, j: sp.nsimplify(round(float(A[i, j]), 10), rational=True))


# ==================================================
def load_supergroup():
    """
    Load the direct supergroups (multipiegen/binary/supergroup.pkl).

    Returns:
        - (BinaryManager) -- {PG_id/SG_id: [{"parent", "A", "index", "index_parent"}]}.
    """
    from multipie.util.util_binary import BinaryManager

    from multipiegen import __bin_dir__

    return BinaryManager(supergroup_file, topdir=__bin_dir__)
