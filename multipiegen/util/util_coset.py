# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Cosets of a group in its parents: images of the cosets, the first coset, the parent structures, and the cosets of
all parents along the chain of representative parents.

Groups and symmetry operations are taken from the tables of this repository only (util_group_data,
util_symmetry_operation.group_symmetry_operation, data_parent_pg/sg, data_wyckoff_site). Numbers are handled as
floats and converted to fractions (denominator <= 48) or to exact expressions with sqrt when written.

- Direct parent. For an SG H, the parent P of data_parent_sg[H] (the parent that restores the least symmetry); for a
  PG H, every parent in data_parent_pg[H], and the representative parent data_representative_parent_pg[H] (the
  "direct parent" used for the first coset). The "-c" point groups share the operations of the real group.
- Coset image (coset_image_sg, coset_image_pg). P is split into left cosets g H (modulo the lattice of H for SG), and
  the image g H g^-1 is looked for among the children of P (the groups that have P in data_parent_pg/sg) with the
  same axes (for SG, up to an origin shift k). The map M = R'^-1 k^-1 g R (coordinates of H -> those of the image)
  maps each Wyckoff position of H to one of the image, with a map of the parameters.
- First coset (first_coset). The second coset of data_parent_*[H][P] (the first one other than H), with its
  representative g (first operation), the order n of g modulo H (smallest n with g^n in H), and its image.
  g is a rotation, or the inversion (with a translation for SG): "kind" = "rotation" or "inversion".
- Parent structure (parent_limit_sg, parent_limit_pg). For each Wyckoff position of H with free parameters, the
  largest sets of parameter values at which the orbit of H is invariant under P with no coinciding sites, the
  Wyckoff position of P occupied there, and the cosets that keep it.
- Cosets of the parents (parent_chain, coset_list). The chain of representative parents up to the top group (Oh,
  D6h, SG 221-230 or 191-194), and the left cosets of H in each of them (tags of the operations of the parent; for SG
  with the translation in the coordinates of H modulo its lattice).
"""

import itertools
from fractions import Fraction

import numpy as np
import sympy as sp
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager

from multipiegen import __bin_dir__

from multipiegen.data.data_parent_pg import data_parent_pg, data_representative_parent_pg
from multipiegen.data.data_parent_sg import data_parent_sg, data_representative_parent_sg
from multipiegen.data.data_wyckoff_site import wyckoff_site_data
from multipiegen.util.util_group_data import pg_id, sg_tag
from multipiegen.util.util_symmetry_operation import group_symmetry_operation

# ==================================================
TOL = 1e-6
PG_TOP = ["Oh", "D6h"]
SG_TOP = [221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 191, 192, 193, 194]
CENTERING = {
    "P": [[0, 0, 0]],
    "A": [[0, 0, 0], [0, 1 / 2, 1 / 2]],
    "C": [[0, 0, 0], [1 / 2, 1 / 2, 0]],
    "I": [[0, 0, 0], [1 / 2, 1 / 2, 1 / 2]],
    "F": [[0, 0, 0], [0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]],
    "R": [[0, 0, 0], [2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]],
}
X, Y, Z = sp.symbols("x y z")
PARAMS = [X, Y, Z]
coset_file = "coset"
_cache = {}


# ==================================================
def frac(v, den=48):
    """float -> fraction string (denominator <= den)."""
    f = Fraction(float(v)).limit_denominator(den)
    assert abs(float(f) - float(v)) < 1e-6, v
    return str(f)


# ==================================================
def num_str(v):
    """float -> fraction string, or an exact expression with sqrt(2), sqrt(3), sqrt(6) when it is not rational."""
    f = Fraction(float(v)).limit_denominator(48)
    if abs(float(f) - float(v)) < 1e-9:
        return str(f)
    e = sp.nsimplify(float(v), [sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)], tolerance=1e-10)
    assert abs(float(e) - float(v)) < 1e-8, v
    return str(e)


# ==================================================
def frac_list(a, exact=False):
    a = np.asarray(a)
    f = num_str if exact else frac
    if a.ndim == 1:
        return "[" + ",".join(f(x) for x in a) + "]"
    return "[" + ",".join(frac_list(r, exact) for r in a) + "]"


# ==================================================
def to_float(s):
    """string (sympy expression or matrix) -> float ndarray."""
    return np.array(sp.Matrix(str_to_sympy(s)).evalf(), dtype=float)


# ==================================================
def mod1(v):
    v = np.asarray(v, dtype=float)
    r = v - np.floor(v + TOL)
    r[np.abs(r - 1) < TOL] = 0.0
    r[np.abs(r) < TOL] = 0.0
    return r


# ==================================================
def sg_number_tag(no):
    return sg_tag[f"SG:{no}"]


# ==================================================
def sg_ops(no):
    """affine 4x4 (fractional) operations of SG no. (point operations x centering)."""
    so = group_symmetry_operation(sg_number_tag(no))
    ops = [np.array(sp.Matrix(o).evalf(), dtype=float) for o in so["fractional"]]
    cent = CENTERING[so["lattice"]]
    out = []
    for c in cent:
        for o in ops:
            o = o.copy()
            o[:3, 3] = mod1(o[:3, 3] + np.array(c))
            out.append(o)
    return out, so["tag"], so["lattice"]


# ==================================================
def pg_ops(tag):
    """3x3 (fractional) operations of PG tag."""
    so = group_symmetry_operation(tag)
    return [np.array(sp.Matrix(o).evalf(), dtype=float) for o in so["fractional"]], so["tag"]


# ==================================================
def pg_real(tag):
    """the "-c" (complex) settings share the operations and the chain of the real group."""
    return tag[:-2] if tag.endswith("-c") else tag


# ==================================================
def op_key(W, Linv=None):
    """canonical key of an operation; translations reduced modulo the lattice with basis Linv^-1."""
    rot = tuple(np.round(W[:3, :3], 6).ravel()) if W.shape[0] == 4 else tuple(np.round(W, 6).ravel())
    if W.shape[0] == 3:
        return rot
    u = W[:3, 3] if Linv is None else Linv @ W[:3, 3]
    return rot, tuple(np.round(mod1(u), 5))


# ==================================================
def group_key(ops, Linv=None):
    return frozenset(op_key(W, Linv) for W in ops)


# ==================================================
def wyckoff_positions(group):
    """{letter: (list of sympy position vectors, site symmetry)} for a PG tag or SG no."""
    key = pg_real(group) if isinstance(group, str) else sg_number_tag(group)
    out = {}
    for w, (pos, ss) in wyckoff_site_data[key].items():
        lst = str_to_sympy(pos)
        out[w] = ([sp.Matrix(p) for p in lst], ss)
    return out


# ==================================================
def affine_parts(p):
    """position (sympy 3-vector in x, y, z) -> (A 3x3, b 3) with p = A (x,y,z) + b."""
    A = np.array([[float(sp.diff(pi, s)) for s in PARAMS] for pi in p], dtype=float)
    b = np.array([float(pi.subs({X: 0, Y: 0, Z: 0})) for pi in p], dtype=float)
    return A, b


# ==================================================
def params_of(p):
    return [s for s in PARAMS if any(pi.has(s) for pi in p)]


# ==================================================
def match_position(q_A, q_b, r_A, r_b, cent, periodic=True):
    """
    Find u = P x + d with r_A u + r_b + c + n = q_A x + q_b (c in cent, n integer if periodic).
    Returns (P, d) or None.
    """
    P, res, *_ = np.linalg.lstsq(r_A, q_A, rcond=None)
    if np.abs(r_A @ P - q_A).max() > TOL:
        return None
    for c in cent:
        v = q_b - r_b - np.asarray(c)
        if not periodic:
            d = np.linalg.lstsq(r_A, v, rcond=None)[0]
            if np.abs(r_A @ d - v).max() < TOL:
                return P, d
            continue
        for n in itertools.product((-2, -1, 0, 1, 2), repeat=3):
            w = v - np.array(n)
            d = np.linalg.lstsq(r_A, w, rcond=None)[0]
            if np.abs(r_A @ d - w).max() < TOL:
                return P, mod1(d) if np.abs(d).max() > 0 else d
    return None


# ==================================================
def param_map_string(P, d, params_out, params_in):
    """{'x': '-x', ...}: parameters of the image representative in terms of the parameters of the source."""
    out = {}
    for i, s in enumerate(PARAMS):
        if s not in params_out:
            continue
        expr = sum(sp.sympify(num_str(P[i, j])) * PARAMS[j] for j in range(3) if PARAMS[j] in params_in) + sp.sympify(
            num_str(d[i])
        )
        out[str(s)] = str(expr)
    return out


# ==================================================
def wyckoff_image(M, src, dst, cent_dst, periodic=True):
    """
    Image of every Wyckoff position of `src` under the coordinate map M (src -> dst coordinates).
    Returns {letter: (letter', {param of letter' rep: expression in params of letter rep})}.
    """
    wsrc, wdst = wyckoff_positions(src), wyckoff_positions(dst)
    out = {}
    for w, (pos, _) in wsrc.items():
        p = pos[0]
        A, b = affine_parts(p)
        qA = M[:3, :3] @ A
        qb = M[:3, :3] @ b + (M[:3, 3] if M.shape[0] == 4 else 0)
        mult = len(pos) * (len(cent_dst) if periodic else 1)
        found = None
        for w2, (pos2, _) in wdst.items():
            if len(pos2) * (len(cent_dst) if periodic else 1) != mult:
                continue
            for r in pos2:
                rA, rb = affine_parts(r)
                m = match_position(qA, qb, rA, rb, cent_dst if periodic else [[0, 0, 0]], periodic)
                if m is not None:
                    # express in the parameters of the representative r0 of w2: r = r0 under some op, so use r's params
                    found = (w2, param_map_string(m[0], m[1], params_of(r), params_of(p)), int(pos2.index(r)))
                    break
            if found:
                break
        out[w] = found
    return out


# ==================================================
def _sg_elements_mod_child(top, Rc):
    """elements of the top SG modulo the child lattice (top coordinates)."""
    ops, _, _ = sg_ops(top)
    Linv = np.linalg.inv(Rc[:3, :3])
    trans = {}
    for n in itertools.product(range(0, 4), repeat=3):
        u = mod1(Linv @ np.array(n, float))
        trans[tuple(np.round(u, 5))] = Rc[:3, :3] @ u
    els = {}
    for W in ops:
        for t in trans.values():
            Wt = W.copy()
            Wt[:3, 3] += t
            els.setdefault(op_key(Wt, Linv), Wt)
    return list(els.values())


# ==================================================
def _in_top(ops, R):
    Ri = np.linalg.inv(R)
    return [R @ W @ Ri for W in ops]


# ==================================================
def _find_shift(target_ops, cand_ops, Linv, Lc):
    """origin shifts k (top coordinates) with k H' k^-1 == target (all modulo the child lattice)."""
    tgt = {}
    for W in target_ops:
        tgt.setdefault(tuple(np.round(W[:3, :3], 6).ravel()), []).append(W[:3, 3])
    by_rot = {}
    for W in cand_ops:
        r = tuple(np.round(W[:3, :3], 6).ravel())
        if r not in tgt:
            return []
        by_rot.setdefault(r, W)
    if set(by_rot) != set(tgt):
        return []
    # grid of shifts in child coordinates (1/24), vectorized test on a few operations
    g = np.array(list(itertools.product(np.arange(24) / 24, repeat=3)))
    s_top = g @ Lc.T
    ok = np.ones(len(g), bool)
    for r, W in list(by_rot.items())[: min(4, len(by_rot))]:
        Wr = W[:3, :3]
        # k W k^-1 = (W, t + (1-W) s)
        new_t = W[:3, 3][None, :] + s_top @ (np.eye(3) - Wr).T
        hit = np.zeros(len(g), bool)
        for t in tgt[r]:
            d = (new_t - t[None, :]) @ Linv.T
            hit |= np.all(np.abs(d - np.rint(d)) < 1e-5, axis=1)
        ok &= hit
    out = []
    tkey = group_key(target_ops, Linv)
    for s in s_top[ok]:
        k = np.eye(4)
        k[:3, 3] = s
        ki = np.linalg.inv(k)
        if group_key([k @ W @ ki for W in cand_ops], Linv) == tkey:
            out.append(s)
    return out


# ==================================================
def _tag_of(top_ops, top_tags, W):
    for t, O in zip(top_tags, top_ops):
        if np.allclose(O[:3, :3], W[:3, :3], atol=TOL):
            return t.split(":")[0]
    return "?"


# ==================================================
def _children_sg(P):
    return sorted(c for c, v in data_parent_sg.items() if P in v)


# ==================================================
def _children_pg(P):
    return [c for c, v in data_parent_pg.items() if P in v]


# ==================================================
def coset_image_min_sg(no, P):
    R = to_float(data_parent_sg[no][P][1])
    Linv = np.linalg.inv(R[:3, :3])
    Hops_c, _, _ = sg_ops(no)
    Hops = _in_top(Hops_c, R)
    els = _sg_elements_mod_child(P, R)
    p_ops, p_tags, _ = sg_ops(P)
    cosets, seen = [], set()
    for g in sorted(els, key=lambda W: (0 if np.allclose(W[:3, :3], np.eye(3)) else 1, np.abs(W[:3, 3]).sum())):
        if op_key(g, Linv) in seen:
            continue
        seen |= {op_key(g @ h, Linv) for h in Hops}
        cosets.append(g)
    order = len(group_key(Hops, Linv))
    cands = []
    for s in _children_sg(P):
        Rs = to_float(data_parent_sg[s][P][1])
        M = np.linalg.inv(R[:3, :3]) @ Rs[:3, :3]
        if not (np.allclose(M, np.rint(M), atol=1e-8) and abs(abs(np.linalg.det(M)) - 1) < 1e-8):
            continue
        ops_t = _in_top(sg_ops(s)[0], Rs)
        if len(group_key(ops_t, Linv)) == order:
            cands.append((s, Rs, ops_t))
    result = []
    for g in cosets:
        gi = np.linalg.inv(g)
        target = [g @ h @ gi for h in Hops]
        best = None
        for s, Rs, ops_t in cands:
            shifts = _find_shift(target, ops_t, Linv, R[:3, :3])
            if not shifts:
                continue
            shifts = sorted(
                shifts, key=lambda v: (round(float(np.linalg.norm(mod1(Linv @ v))), 6), tuple(np.round(mod1(Linv @ v), 6)))
            )
            for sh in shifts[:24]:
                k = np.eye(4)
                k[:3, 3] = sh
                Mmap = np.linalg.inv(Rs) @ np.linalg.inv(k) @ g @ R
                Mmap[:3, 3] = mod1(Mmap[:3, 3])
                wy = wyckoff_image(Mmap, no, s, CENTERING[sg_ops(s)[2]])
                score = sum(1 for w, v in wy.items() if v is not None and v[0] == w)
                key = (-(s == no), -score, float(np.linalg.norm(mod1(Linv @ sh))))
                if best is None or key < best[0]:
                    best = (key, s, sh, Mmap, wy)
            if best is not None and best[1] == no:
                break
        rep = (_tag_of(p_ops, p_tags, g), frac_list(mod1(Linv @ g[:3, 3])))
        if best is None:
            result.append({"rep": rep, "image": None})
            continue
        _, s, sh, Mmap, wy = best
        trivial = s == no and op_key(Mmap) in {op_key(W) for W in Hops_c}
        result.append(
            {
                "rep": rep,
                "image": int(s),
                "shift": frac_list(sh),
                "map": frac_list(Mmap),
                "trivial": bool(trivial),
                "wyckoff": {w: (v[0], v[1]) if v is not None else None for w, v in wy.items()},
            }
        )
    return {"R": data_parent_sg[no][P][1], "coset": result}


# ==================================================
def _exact_ops(tag):
    so = group_symmetry_operation(tag)
    return [sp.Matrix(o) for o in so["fractional"]], so["tag"]


# ==================================================
def coset_image_min_pg(tag, P):
    real = pg_real(tag)
    Rs_ex = sp.Matrix(str_to_sympy(data_parent_pg[real][P][1]))
    R = np.array(Rs_ex.evalf(), dtype=float)
    Hc, _ = pg_ops(real)
    Hops = [R @ W @ np.linalg.inv(R) for W in Hc]
    p_ops, p_tags = pg_ops(P)
    p_ex, _ = _exact_ops(P)
    cosets, seen = [], set()
    for g, t, ge in zip(p_ops, p_tags, p_ex):
        if op_key(g) in seen:
            continue
        seen |= {op_key(g @ h) for h in Hops}
        cosets.append((g, t, ge))
    cands = []
    for s in _children_pg(P):
        Rp_ex = sp.Matrix(str_to_sympy(data_parent_pg[s][P][1]))
        Rp = np.array(Rp_ex.evalf(), dtype=float)
        ops = [Rp @ W @ np.linalg.inv(Rp) for W in pg_ops(s)[0]]
        cands.append((s, Rp, Rp_ex, group_key(ops)))
    complex_child = tag.endswith("-c")
    result = []
    for g, t, ge in cosets:
        gi = np.linalg.inv(g)
        target = group_key([g @ h @ gi for h in Hops])
        found = [c for c in cands if c[3] == target]
        if not found:
            result.append({"rep": t, "image": None})
            continue
        # prefer the child itself
        found.sort(key=lambda c: c[0] != real)
        s, Rp, Rp_ex, _ = found[0]
        M = np.linalg.inv(Rp) @ g @ R
        M_ex = sp.simplify(Rp_ex.inv() * ge * Rs_ex)
        wy = wyckoff_image(M, real, s, [[0, 0, 0]], periodic=False)
        img = s + "-c" if (complex_child and s + "-c" in pg_id) else s
        result.append(
            {
                "rep": t,
                "image": img,
                "map": str(M_ex.tolist()).replace(" ", ""),
                "wyckoff": {w: (v[0], v[1]) if v is not None else None for w, v in wy.items()},
            }
        )
    return {"R": data_parent_pg[real][P][1], "coset": result}


# ==================================================
def _standard(x0, N):
    """
    standard form of the affine subspace x0 + span(N) modulo integer shifts of the parameters:
    the free parameters (an invertible block of N, chosen in order) are set to 0 in x0, the others are reduced
    modulo 1, and the directions are written as x = x0 + T x_free.  Equal keys mean equal subspaces (duplicates that
    differ by T n_free may remain, which is harmless).
    """
    k = len(x0)
    if N.shape[1] == 0:
        return mod1(x0), N
    free = None
    for f in itertools.combinations(range(k), N.shape[1]):
        if abs(np.linalg.det(N[list(f), :])) > 1e-8:
            free = f
            break
    T = N @ np.linalg.inv(N[list(free), :])
    x = x0 - T @ x0[list(free)]
    dep = [i for i in range(k) if i not in free]
    x[dep] = mod1(x[dep])
    x[list(free)] = 0.0
    return x, T


# ==================================================
def _subspace_key(x0, N):
    x, T = _standard(x0, N)
    return (N.shape[1], tuple(np.round(x, 5)), tuple(np.round(T, 5).ravel()))


# ==================================================
def _solve(M, v):
    """affine solution set of M x = v (or None); x0 is the minimum-norm solution, N a basis of the null space."""
    k = M.shape[1]
    if k == 0:
        return (np.zeros(0), np.zeros((0, 0))) if np.abs(v).max() < 1e-7 else None
    u, s, vh = np.linalg.svd(M)
    tol = 1e-6 * max(1.0, float(s.max()) if s.size else 1.0)
    rank = int((s > tol).sum())
    x0 = vh[:rank].T @ ((u[:, :rank].T @ v) / s[:rank]) if rank else np.zeros(k)
    if np.abs(M @ x0 - v).max() > 1e-7:
        return None
    N = vh[rank:].T
    return x0, N


# ==================================================
def _canon(x0, N):
    """remove from x0 its component along the subspace directions."""
    if N.shape[1] == 0:
        return x0, N
    Q, _ = np.linalg.qr(N)
    return x0 - Q @ (Q.T @ x0), Q


# ==================================================
def _restrict(sub, C, Hops, A, b, periodic, full=True):
    """restrict the subspace sub=(x0,N) by the condition C p(x) in orbit_H(p(x))."""
    x0, N = sub
    out = []
    for W in Hops:
        if periodic:
            D = C[:3, :3] - W[:3, :3]
            const = D @ b + C[:3, 3] - W[:3, 3]
        else:
            D = C - W
            const = D @ b
        MA = D @ A  # 3 x k
        Mt = MA @ N
        base = -(MA @ x0 + const)
        if not periodic:
            r = _solve(Mt, base)
            if r is not None:
                out.append(_canon(x0 + N @ r[0], N @ r[1]))
            continue
        # integer n with Mt t = base + n: n = Mt t - base with the parameters t in a few unit cells, so around -base.
        # full = False searches around +base only (the first version of the tables, see _valid_limits).
        w = np.abs(Mt).sum(axis=1) * 2 + 1
        ns = []
        for c in (base, -base) if full else (base,):
            lo, hi = np.floor(c - w).astype(int), np.ceil(c + w).astype(int)
            ns += [n for n in itertools.product(*[range(l, h + 1) for l, h in zip(lo, hi)]) if n not in ns]
        for n in ns:
            r = _solve(Mt, base + np.array(n))
            if r is not None:
                out.append(_canon(x0 + N @ r[0], N @ r[1]))
    return out


# ==================================================
def _limits(Hops, gens, A, b, k, periodic, full=True):
    subs = [(np.zeros(k), np.eye(k))]
    for C in gens:
        subs = [_standard(*s) if periodic else s for s in subs]
        new = {}
        for sub in subs:
            for s in _restrict(sub, C, Hops, A, b, periodic, full):
                if periodic:
                    s = _standard(*s)
                new.setdefault(_subspace_key(*s), s)
        subs = list(new.values())
        if not subs:
            break
    # keep proper and maximal subspaces with representative points in [0,1)
    subs = [s for s in subs if s[1].shape[1] < k]
    out = []
    for s in sorted(subs, key=lambda s: -s[1].shape[1]):
        contained = False
        for t in out:
            if t[1].shape[1] > s[1].shape[1]:
                Q, _ = np.linalg.qr(t[1])
                d = s[0] - t[0]
                d = d - Q @ (Q.T @ d)
                if np.allclose(d - np.rint(d), 0, atol=1e-6) and np.allclose((np.eye(k) - Q @ Q.T) @ s[1], 0, atol=1e-6):
                    contained = True
        if not contained:
            out.append(s)
    return out


# ==================================================
def _invariant(Hops, Pops, A, b, sub, k, periodic, mult):
    """check at two generic points of sub that the orbit is invariant under all of Pops and has the multiplicity."""
    for trial in (0, 1):
        t = 0.1234 + 0.0731 * np.arange(sub[1].shape[1]) + 0.05 * trial
        x = sub[0] + (sub[1] @ t if sub[1].shape[1] else 0)
        q = A @ x + b
        orb = _orbit(Hops, q, [[0, 0, 0]], periodic)
        if len(orb) != mult:
            return False
        for W in Pops:
            qq = W[:3, :3] @ q + W[:3, 3] if periodic else W @ q
            if tuple(np.round(mod1(qq) if periodic else qq, 5)) not in orb:
                return False
    return True


# ==================================================
def _right_coset_reps(Hops, Pels, periodic):
    """representatives c of the right cosets H c of H in the parent (elements given in child coordinates)."""
    seen, reps = set(), []
    for c in Pels:
        if op_key(c) in seen:
            continue
        seen |= {op_key(h @ c) if periodic else op_key(h @ c) for h in Hops}
        reps.append(c)
    return reps


# ==================================================
def _valid_limits(Hops, Pels, A, b, k, periodic, mult):
    """
    limits for the parent with elements Pels (child coordinates, modulo the child lattice): the orbit of p(x) is
    invariant under the parent iff c p(x) lies on the orbit of H for every representative c of the right cosets H c
    (any element is h c, and h c p in H p <=> c p in H p), so the conditions on the representative point suffice.

    the limits found with the search window of the first version of the tables come first in their order; those
    found only with the full window follow, unless the orbit of H at their generic point is that of one listed before
    (the same parent structure).
    """
    reps = [c for c in _right_coset_reps(Hops, Pels, periodic) if op_key(c) not in {op_key(h) for h in Hops}]
    uniq = {}
    for full in (False, True):
        first = list(uniq.values())
        for sub in _limits(Hops, reps, A, b, k, periodic, full):
            key = _subspace_key(*sub)
            if key in uniq or not _invariant(Hops, Pels, A, b, sub, k, periodic, mult):
                continue
            if full and any(_same_orbits(Hops, A, b, sub, s, periodic) for s in uniq.values()):
                continue
            uniq[key] = sub
        if not periodic:
            break
    return list(uniq.values())


# ==================================================
def _same_orbits(Hops, A, b, s1, s2, periodic):
    """the orbit of H at a generic point of s1 has a point on s2 (s1 and s2 of the same dimension) ?"""
    if s1[1].shape[1] != s2[1].shape[1]:
        return False
    q = A @ _generic(s1) + b
    x0, N = s2
    AN = A @ N if N.shape[1] else np.zeros((3, 0))
    for W in Hops:
        r = W[:3, :3] @ q + (W[:3, 3] if periodic else 0)
        for n in itertools.product((-1, 0, 1), repeat=3) if periodic else [(0, 0, 0)]:
            v = r + np.array(n) - (A @ x0 + b)
            t = np.linalg.lstsq(AN, v, rcond=None)[0] if N.shape[1] else np.zeros(0)
            d = v - (AN @ t if N.shape[1] else 0)
            if np.abs(d - (np.rint(d) if periodic else 0)).max() < 1e-6:
                return True
    return False


# ==================================================
def _param_dict(sub, params):
    """parameter values: free parameters (pivots of N) are kept, the others expressed by them."""
    x0, N = sub
    k = len(params)
    if N.shape[1] == 0:
        return {str(p): num_str(mod1(x0)[i]) for i, p in enumerate(params)}
    # choose free parameters: columns of an invertible square block of N
    for free in itertools.combinations(range(k), N.shape[1]):
        B = N[list(free), :]
        if abs(np.linalg.det(B)) > 1e-8:
            break
    T = N @ np.linalg.inv(B)  # x = x0' + T x_free
    x0p = x0 - T @ x0[list(free)]
    out = {}
    for i, p in enumerate(params):
        if i in free:
            out[str(p)] = str(p)
        else:
            expr = sum(sp.sympify(num_str(T[i, j])) * params[f] for j, f in enumerate(free)) + sp.sympify(num_str(mod1(x0p)[i]))
            out[str(p)] = str(expr)
    return out, free


# ==================================================
def _orbit(ops, p, cent, periodic):
    pts = set()
    for W in ops:
        q = W[:3, :3] @ p + (W[:3, 3] if periodic else 0)
        for c in cent:
            pts.add(tuple(np.round(mod1(q + np.array(c)) if periodic else q, 5)))
    return frozenset(pts)


# ==================================================
def _parent_wyckoff(P, q, cent, periodic):
    for w, (pos, _) in wyckoff_positions(P).items():
        for r in pos:
            rA, rb = affine_parts(r)
            if match_position(np.zeros((3, 3)), q, rA, rb, cent, periodic) is not None:
                return w
    return None


# ==================================================
def _pos_str(p, pdict):
    sub = {sp.Symbol(k): sp.sympify(v) for k, v in pdict.items()}
    return "[" + ",".join(str(sp.nsimplify(sp.sympify(str(pi)).subs(sub))) for pi in p) + "]"


# ==================================================
def _generic(sub):
    t = 0.1234 + 0.0731 * np.arange(sub[1].shape[1])
    return sub[0] + (sub[1] @ t if sub[1].shape[1] else 0)


# ==================================================
def parent_limit_min_sg(no, P, cosets):
    R = to_float(data_parent_sg[no][P][1])
    Ri = np.linalg.inv(R)
    Hops, _, lat = sg_ops(no)
    centH = CENTERING[lat]
    Pops, _, latP = sg_ops(P)
    trans = {}
    for n in itertools.product(range(0, 4), repeat=3):
        u = mod1(Ri[:3, :3] @ np.array(n, float))
        trans[tuple(np.round(u, 5))] = u
    Pels = {}
    for W in [Ri @ W @ R for W in Pops]:
        for u in trans.values():
            Wt = W.copy()
            Wt[:3, 3] = mod1(Wt[:3, 3] + u)
            Pels.setdefault(op_key(Wt), Wt)
    Pels = list(Pels.values())
    out = {}
    for w, (pos, _) in wyckoff_positions(no).items():
        p = pos[0]
        params = params_of(p)
        if not params:
            continue
        A, b = affine_parts(p)
        A = A[:, [PARAMS.index(s) for s in params]]
        lims = []
        for sub in _valid_limits(Hops, Pels, A, b, len(params), True, len(pos) * len(centH)):
            try:
                pd = _param_dict(sub, params)
            except AssertionError:
                print(f"  warning: SG {no} -> {P} {w}: skipped a non-rational solution", flush=True)
                continue
            pdict = pd if isinstance(pd, dict) else pd[0]
            q = A @ _generic(sub) + b
            qP = R[:3, :3] @ q + R[:3, 3]
            orb = _orbit(Hops, q, [[0, 0, 0]], True)
            kept = []
            for ci, c in enumerate(cosets):
                if c.get("image") is None:
                    continue
                M = to_float(c["map"])
                if _orbit(sg_ops(c["image"])[0], M[:3, :3] @ q + M[:3, 3], [[0, 0, 0]], True) == orb:
                    kept.append(ci)
            lims.append(
                {
                    "param": pdict,
                    "position": _pos_str(p, pdict),
                    "parent_wyckoff": _parent_wyckoff(P, qP, CENTERING[latP], True),
                    "preserved_by": kept,
                }
            )
        if lims:
            out[w] = lims
    return out


# ==================================================
def parent_limit_min_pg(tag, P, cosets):
    real = pg_real(tag)
    R = np.array(sp.Matrix(str_to_sympy(data_parent_pg[real][P][1])).evalf(), dtype=float)
    Ri = np.linalg.inv(R)
    Hops, _ = pg_ops(real)
    Pels = [Ri @ W @ R for W in pg_ops(P)[0]]
    out = {}
    for w, (pos, _) in wyckoff_positions(real).items():
        p = pos[0]
        params = params_of(p)
        if not params:
            continue
        A, b = affine_parts(p)
        A = A[:, [PARAMS.index(s) for s in params]]
        lims = []
        for sub in _valid_limits(Hops, Pels, A, b, len(params), False, len(pos)):
            pd = _param_dict(sub, params)
            pdict = pd if isinstance(pd, dict) else pd[0]
            q = A @ _generic(sub) + b
            orb = _orbit(Hops, q, [[0, 0, 0]], False)
            kept = []
            for ci, c in enumerate(cosets):
                if c.get("image") is None:
                    continue
                M = np.array(sp.Matrix(str_to_sympy(c["map"])).evalf(), dtype=float)
                img = c["image"][:-2] if c["image"].endswith("-c") else c["image"]
                if _orbit(pg_ops(img)[0], M @ q, [[0, 0, 0]], False) == orb:
                    kept.append(ci)
            lims.append(
                {
                    "param": pdict,
                    "position": _pos_str(p, pdict),
                    "parent_wyckoff": _parent_wyckoff(P, R @ q, [[0, 0, 0]], False),
                    "preserved_by": kept,
                }
            )
        if lims:
            out[w] = lims
    return out


# ==================================================


# ==================================================
def _key(W, sg):
    """key of an operation (SG: translation modulo 1)."""
    if sg:
        return (tuple(np.rint(W[:3, :3] * 1e6).astype(int).ravel()), tuple(np.rint(mod1(W[:3, 3]) * 1e5).astype(int) % 100000))
    return tuple(np.rint(W * 1e6).astype(int).ravel())


# ==================================================
def _order(g, Hkeys, sg):
    """smallest n with g^n in H."""
    n, G = 1, g.copy()
    while _key(G, sg) not in Hkeys:
        G = G @ g
        n += 1
        assert n <= 12
    return n


# ==================================================
def _rep_matrix_sg(ent, P, R):
    """representative of a coset image entry (tag of the point part in P, translation in H) in the coordinates of H."""
    tag, t = ent["rep"]
    tv = np.array([float(sp.sympify(a)) for a in t.strip("[]").split(",")])
    for W, tg in zip(*sg_ops(P)[:2]):
        if tg.split(":")[0] == tag:
            gP = W.copy()
            gP[:3, 3] = R[:3, :3] @ tv
            return np.linalg.inv(R) @ gP @ R
    return None


# ==================================================
def coset_image_sg(no, P=None):
    """
    Images of the left cosets of SG no. in its direct parent (cached).

    Args:
        no (int): SG number.
        P (int, optional): parent SG number (the direct parent if None).

    Returns:
        - (dict) -- {"R": R, "coset": [coset]}, coset = {"rep": (tag, translation), "image": SG no. or None, "shift",
          "map", "trivial", "wyckoff": {W: (W', {param: expression})}}.
    """
    P = next(iter(data_parent_sg[no])) if P is None else P
    k = ("SG", no, P)
    if k not in _cache:
        _cache[k] = coset_image_min_sg(no, P)
    return _cache[k]


# ==================================================
def coset_image_pg(tag, P=None):
    """
    Images of the left cosets of a PG in a parent (cached).

    Args:
        tag (str): PG tag.
        P (str, optional): parent PG tag (the representative parent if None).

    Returns:
        - (dict) -- {"R": R, "coset": [coset]}, coset = {"rep": tag, "image": PG tag or None, "map", "wyckoff"}.
    """
    P = data_representative_parent_pg[pg_real(tag)] if P is None else P
    k = ("PG", tag, P)
    if k not in _cache:
        _cache[k] = coset_image_min_pg(tag, P)
    return _cache[k]


# ==================================================
def first_coset_sg(no):
    """
    First coset of SG no. in its direct parent.

    Args:
        no (int): SG number.

    Returns:
        - (dict) -- {"parent", "R", "index", "ops", "ops_child", "rep", "rep_child", "order", "kind", "image", "map",
          "shift", "wyckoff", "coset_index"} (matrices as float ndarray, "coset_index": index of the coset in
          coset_image_sg).
    """
    P = next(iter(data_parent_sg[no]))
    cosets, Rs = data_parent_sg[no][P]
    R = to_float(Rs)
    Ri = np.linalg.inv(R)
    Pops, Ptags, _ = sg_ops(P)
    Hkeys = {_key(W, True) for W in sg_ops(no)[0]}
    child = [Ri @ Pops[i - 1] @ R for i in cosets[1]]
    g = child[0]
    ci = coset_image_sg(no, P)
    idx = None
    for i, ent in enumerate(ci["coset"]):
        rW = _rep_matrix_sg(ent, P, R)
        if rW is not None and _key(np.linalg.inv(rW) @ g, True) in Hkeys:
            idx = i
            break
    img = ci["coset"][idx] if idx is not None and ci["coset"][idx]["image"] is not None else None
    return {
        "parent": P,
        "R": Rs,
        "index": len(cosets),
        "ops": [Ptags[i - 1] for i in cosets[1]],
        "ops_child": child,
        "rep": Ptags[cosets[1][0] - 1],
        "rep_child": g,
        "order": _order(g, Hkeys, True),
        "kind": "rotation" if np.linalg.det(g[:3, :3]) > 0 else "inversion",
        "image": None if img is None else img["image"],
        "map": None if img is None else img["map"],
        "shift": None if img is None else img["shift"],
        "wyckoff": None if img is None else img["wyckoff"],
        "coset_index": idx,
    }


# ==================================================
def first_coset_pg(tag):
    """
    First coset of a PG in its representative parent.

    Args:
        tag (str): PG tag.

    Returns:
        - (dict) -- as first_coset_sg (no "shift"), matrices 3x3.
    """
    real = pg_real(tag)
    P = data_representative_parent_pg[real]
    cosets, Rs = data_parent_pg[real][P]
    R = to_float(Rs)
    Ri = np.linalg.inv(R)
    Pops, Ptags = pg_ops(P)
    Pops = [np.array(W, float)[:3, :3] for W in Pops]
    Hkeys = {_key(np.array(W, float)[:3, :3], False) for W in pg_ops(real)[0]}
    child = [Ri @ Pops[i - 1] @ R for i in cosets[1]]
    g = child[0]
    ci = coset_image_pg(tag, P)
    idx = None
    for i, ent in enumerate(ci["coset"]):
        W = Pops[Ptags.index(ent["rep"])]
        if _key(np.linalg.inv(Ri @ W @ R) @ g, False) in Hkeys:
            idx = i
            break
    img = ci["coset"][idx] if idx is not None and ci["coset"][idx]["image"] is not None else None
    return {
        "parent": P,
        "R": Rs,
        "index": len(cosets),
        "ops": [Ptags[i - 1] for i in cosets[1]],
        "ops_child": child,
        "rep": Ptags[cosets[1][0] - 1],
        "rep_child": g,
        "order": _order(g, Hkeys, False),
        "kind": "rotation" if np.linalg.det(g) > 0 else "inversion",
        "image": None if img is None else img["image"],
        "map": None if img is None else img["map"],
        "wyckoff": None if img is None else img["wyckoff"],
        "coset_index": idx,
    }


# ==================================================
def parent_limit_sg(no, P=None):
    """
    Parent structures of SG no. for its direct parent.

    Args:
        no (int): SG number.
        P (int, optional): parent (the direct parent if None).

    Returns:
        - (dict) -- {W: [{"param", "position", "parent_wyckoff", "preserved_by"}]}.
    """
    P = next(iter(data_parent_sg[no])) if P is None else P
    return parent_limit_min_sg(no, P, coset_image_sg(no, P)["coset"])


# ==================================================
def parent_limit_pg(tag, P=None):
    """
    Parent structures of a PG for a parent.

    Args:
        tag (str): PG tag.
        P (str, optional): parent (the representative parent if None).

    Returns:
        - (dict) -- {W: [{"param", "position", "parent_wyckoff", "preserved_by"}]}.
    """
    P = data_representative_parent_pg[pg_real(tag)] if P is None else P
    return parent_limit_min_pg(tag, P, coset_image_pg(tag, P)["coset"])


# ==================================================
def parent_chain(gid):
    """
    Chain of representative parents of a group up to the top group, with the composed transformations.

    Args:
        gid (str): group id ("PG:n" of a real PG, or "SG:n").

    Returns:
        - (list) -- [(parent id, R (float ndarray, x_parent = R x_H))], the direct parent first ([] for a top group).
    """
    if gid.startswith("SG:"):
        no, R, out = int(gid[3:]), np.eye(4), []
        while no in data_parent_sg:
            p = next(iter(data_parent_sg[no]))
            R = to_float(data_parent_sg[no][p][1]) @ R
            out.append((f"SG:{p}", R))
            no = p
        return out
    tag, R, out = pg_real(_pg_tag(gid)), np.eye(3), []
    while tag in data_representative_parent_pg:
        p = data_representative_parent_pg[tag]
        R = to_float(data_parent_pg[tag][p][1]) @ R
        out.append((pg_id[p], R))
        tag = p
    return out


# ==================================================
def _pg_tag(gid):
    return next(t for t, i in pg_id.items() if i == gid)


# ==================================================
def coset_list(gid, R, parent):
    """
    Left cosets g H of a group in a parent, without H itself and the cosets of translations only.

    Args:
        gid (str): group id ("PG:n" of a real PG, or "SG:n").
        R (ndarray): x_parent = R x_H (parent_chain).
        parent (str): parent id.

    Returns:
        - (list) -- [[tag]] for each coset, its representative first; tag = SO tag of the parent for PG, and
          "point part of the SO tag of the parent:[translation in the coordinates of H modulo its lattice]" for SG.

    Note:
        - the cosets are ordered as in coset_image_sg (rotation parts of the identity first, then by the length of the
          translation), and the operations of a coset in the order of the operations of H.
    """
    if gid.startswith("SG:"):
        no, P = int(gid[3:]), int(parent[3:])
        Linv = np.linalg.inv(R[:3, :3])
        Hops = _in_top(sg_ops(no)[0], R)
        els = _sg_elements_mod_child(P, R)
        p_ops, p_tags, _ = sg_ops(P)
        res, seen = [], set()
        for g in sorted(els, key=lambda W: (0 if np.allclose(W[:3, :3], np.eye(3)) else 1, np.abs(W[:3, 3]).sum())):
            if op_key(g, Linv) in seen:
                continue
            members, ops = set(), []
            for h in Hops:
                W = g @ h
                k = op_key(W, Linv)
                if k not in members:
                    members.add(k)
                    ops.append(f"{_tag_of(p_ops, p_tags, W)}:{frac_list(mod1(Linv @ W[:3, 3]))}")
            seen |= members
            if not np.allclose(g[:3, :3], np.eye(3)):
                res.append(ops)
        return res
    tag, P = pg_real(_pg_tag(gid)), _pg_tag(parent)
    Hops = [R @ W @ np.linalg.inv(R) for W in pg_ops(tag)[0]]
    p_ops, p_tags = pg_ops(P)
    res, seen = [], set()
    for g, t in zip(p_ops, p_tags):
        if op_key(g) in seen:
            continue
        members = [op_key(g @ h) for h in Hops]
        seen |= set(members)
        if t != p_tags[0]:
            res.append([p_tags[[op_key(W) for W in p_ops].index(k)] for k in members])
    return res


# ==================================================
def first_coset_rep(gid):
    """
    Representative of the first coset of a group in its direct parent (without the coset images).

    Args:
        gid (str): group id ("PG:n" or "SG:n").

    Returns:
        - (str) -- SO tag of the parent, None if the group has no parent.
    """
    if gid.startswith("SG:"):
        no = int(gid[3:])
        if no not in data_parent_sg:
            return None
        P = next(iter(data_parent_sg[no]))
        return sg_ops(P)[1][data_parent_sg[no][P][0][1][0] - 1]
    tag = pg_real(_pg_tag(gid))
    if tag not in data_representative_parent_pg:
        return None
    P = data_representative_parent_pg[tag]
    return pg_ops(P)[1][data_parent_pg[tag][P][0][1][0] - 1]


# ==================================================
def load_coset():
    """
    Load coset data (multipiegen/binary/coset.pkl).

    Returns:
        - (BinaryManager) -- {PG_id/SG_id: {"parent", "coset", "first"}}.
    """
    return BinaryManager(coset_file, topdir=__bin_dir__)
