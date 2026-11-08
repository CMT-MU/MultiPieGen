"""Completeness check of the Wyckoff-bond data (MultiPie wyckoff_bond_data).

Independent of the file, for every group and every ordinary Wyckoff position W with site-symmetry group H:
  allowed bond types = H-classes of the nonzero subspaces S = Fix(K), K <= H x {+1,-1} acting on vectors as s*R
  (a vector is identified with its negative), i.e. the classes of vector "strata" of R^3.
and compares them with the entries of the file (vector family of each entry, mapped back to the site rep).
"""
import itertools
import pickle
import re
import sys
import collections

sys.path.insert(0, "/home/claude")
sys.path.insert(0, "/tmp/claude-0/-home-claude/ebc6dc5f-4f46-5d72-ab5d-2c96f63a86f3/scratchpad")
import sympy as sp
import spglib

import site_symmetry as S
import ita_data as D

data = pickle.load(open("data.pkl", "rb"))
x, y, z, X, Y, Z = S.x, S.y, S.z, S.X, S.Y, S.Z
I3 = sp.eye(3)

PG_START = {"C1": 1, "Ci": 2, "C2": 3, "Cs": 6, "C2h": 10, "D2": 16, "C2v": 25, "D2h": 47, "C4": 75, "S4": 81,
            "C4h": 83, "D4": 89, "C4v": 99, "D2d": 111, "D4h": 123, "C3": 143, "C3i": 147, "D3": 149, "C3v": 156,
            "D3d": 162, "C6": 168, "C3h": 174, "C6h": 175, "D6": 177, "C6v": 183, "D3h": 187, "D6h": 191, "T": 195,
            "Th": 200, "O": 207, "Td": 215, "Oh": 221}


def sg_number(name):
    pg, n = name.split("^")
    return PG_START[pg] + int(n) - 1


def lattice_letter(num):
    h = next(k for k in range(1, 531) if spglib.get_spacegroup_type(k).number == num)
    return spglib.get_spacegroup_type(h).international_short[0]


# ------------------------------------------------------------------ subspace utilities (canonical RREF rows)
def canon(vecs):
    vecs = [sp.Matrix(v) for v in vecs]
    if not vecs:
        return ()
    M = sp.Matrix.hstack(*vecs).T
    rref, piv = M.rref()
    return tuple(tuple(rref.row(i)) for i in range(len(piv)))


def cols(sub):
    return [sp.Matrix(r) for r in sub]


def fix_space(M):
    return canon((M - I3).nullspace())


def intersect(a, b):
    if not a or not b:
        return ()
    A, B = sp.Matrix(a).T, sp.Matrix(b).T
    ns = A.row_join(-B).nullspace()
    return canon([A * v[:A.shape[1], :] for v in ns])


def image(R, sub):
    return canon([R * v for v in cols(sub)])


FULL = canon([sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])])


def closure(subs):
    subs = set(subs)
    changed = True
    while changed:
        changed = False
        for a, b in itertools.combinations(list(subs), 2):
            c = intersect(a, b)
            if c not in subs:
                subs.add(c); changed = True
    return subs


_cache = {}


def vector_strata(Rs):
    """H-classes of nonzero fixed subspaces of H x {+-1} on R^3. -> list of sets of canonical subspaces."""
    key = tuple(sorted(tuple(R) for R in Rs))
    if key in _cache:
        return _cache[key]
    subs = {FULL}
    for R in Rs:
        for s in (1, -1):
            subs.add(fix_space(s * R))
    subs = {s for s in closure(subs) if s}
    return _cache.setdefault(key, orbit_classes(subs, Rs))


def point_strata(Rs):
    """Site strata of a point group: classes of fixed subspaces (zero subspace included)."""
    subs = {FULL} | {fix_space(R) for R in Rs}
    subs = closure(subs)
    out = []
    for cl in orbit_classes(subs, Rs):
        any_s = next(iter(cl))
        K = sum(1 for R in Rs if all((R - I3) * v == sp.zeros(3, 1) for v in cols(any_s)))
        out.append((cl, len(Rs) // K))
    return out


def orbit_classes(subs, Rs):
    left, classes = set(subs), []
    while left:
        s = min(left, key=str)
        cl = {image(R, s) for R in Rs}
        classes.append(cl)
        left -= cl
    return classes


def form(sub):
    """Readable vector family, e.g. ((1,0,0),(0,0,1)) -> [X,0,Z]."""
    if not sub:
        return "0"
    par = "XYZ"
    comp = [[] for _ in range(3)]
    for row in sub:
        piv = next(i for i, e in enumerate(row) if e != 0)
        for i, e in enumerate(row):
            if e != 0:
                comp[i].append(f"{e}*{par[piv]}")
    return "[" + ",".join("+".join(c).replace("1*", "").replace("+-", "-") if c else "0" for c in comp) + "]"


def entry_subspace(vec_expr):
    jac = sp.Matrix([[sp.diff(e, p) for p in (X, Y, Z)] for e in vec_expr])
    assert all(e.subs({X: 0, Y: 0, Z: 0}) == 0 for e in vec_expr), "constant part in vector"
    return canon([jac[:, j] for j in range(3)])


# ------------------------------------------------------------------ ops reconstruction
def rebuild_ops(d):
    """The general-position entry (first center == [x,y,z], all mapping singletons) gives all operations."""
    best = None
    for lab, (v, c, m) in d.items():
        cs = S._parse(c)
        if all(len(g) == 1 for g in m) and [str(e) for e in cs[0]] == ["x", "y", "z"]:
            if best is None or len(cs) > len(best[1]):
                best = (lab, cs, S._parse(v), m)
    if best is None:
        return None
    lab, cs, vs, m = best
    ops = []
    for j, (c, v) in enumerate(zip(cs, vs)):
        M = sp.eye(4)
        for i, e in enumerate(c):
            e = sp.expand(e)
            for k, p in enumerate((x, y, z)):
                M[i, k] = e.coeff(p)
            M[i, 3] = e.subs({x: 0, y: 0, z: 0})
        ok = list(M[:3, :3] * sp.Matrix([X, Y, Z])) == [sp.expand(e) for e in v] and m[j] == [j + 1]
        ops.append((M, ok))
    return lab, [o for o, _ in ops], all(ok for _, ok in ops)


# ------------------------------------------------------------------ main check
def check_group(name):
    d = data[name]
    res = dict(name=name, problems=[], per_W={})
    rb = rebuild_ops(d)
    if rb is None:
        res["problems"].append("no general-position entry: operations cannot be rebuilt")
        return res
    glab, ops, ops_ok = rb
    res["n_ops"] = len(ops)
    if not ops_ok:
        res["problems"].append("general entry: vectors/mapping inconsistent with its centers")
    space = "^" in name
    if space:
        num = sg_number(name)
        L = lattice_letter(num)
        three_z = any(S._op_type(M[:3, :3]) == "3" and S._axis(M[:3, :3], "3") == (0, 0, 1) for M, _ in [(o, 0) for o in ops])
        centering = L if L in "ABCIF" or (L == "R" and three_z) else None
        shifts = S._centering_shifts(centering)
        k = len(shifts) if shifts else 1
        res.update(number=num, centering=centering, k=k)
        # does the rebuilt operation list agree with the ITA numbering (first n ops of the full list)?
        ita = D.ops_for(num)
        n = len(ita) // (len(S._centering_shifts(L if L != "P" else None) or [0]))
        res["same_as_ita"] = len(ops) == n and all(a == b for a, b in zip(ops, ita[:n]))
        # ordinary Wyckoff positions: ITA list vs labels in the file
        wl = D._wl()[num]
        nW = len(wl)
        ita_W = {f"{len(c)}{D.LETTERS[nW - 1 - i]}" for i, c in enumerate(wl)}
        file_W = {lab.split("@")[1] for lab in d}
        res["W_missing"] = sorted(ita_W - file_W)
        res["W_extra"] = sorted(file_W - ita_W)
    else:
        centering, shifts, k = None, None, 1
        res.update(centering=None, k=1)

    Rs_all = [M[:3, :3] for M in ops]
    byW = collections.OrderedDict()
    for lab, ent in d.items():
        byW.setdefault(lab.split("@")[1], []).append((lab, ent))

    if not space:
        # independent list of site strata of the point group
        ps = point_strata(list({tuple(R): R for R in Rs_all}.values()))
        indep = sorted(m for _, m in ps)
        filem = sorted(int(re.match(r"\d+", w).group()) for w in byW)
        res["W_indep_mult"] = indep
        res["W_file_mult"] = filem
        if indep != filem:
            res["problems"].append(f"ordinary Wyckoff multiplicities: independent {indep} vs file {filem}")

    for W, ents in byW.items():
        lab0, (v0, c0, m0) = ents[0]
        c_ref = S._parse(c0)[0]
        cref_M = sp.Matrix(c_ref)
        r = S.site_symmetry(ops, [list(c_ref)], lattice_mod=space, centering=centering, oriented_errors_ok=True)
        H = [ops[i - 1][:3, :3] for i in r.stabilizer]
        H = list({tuple(R): R for R in H}.values())
        classes = vector_strata(H)
        hit = collections.defaultdict(list)
        unknown = []
        for lab, (v, c, m) in ents:
            ce = S._parse(c)[0]
            # relate the entry's first center to the rep center by an operation g
            g = None
            for M in ops:
                R, t = M[:3, :3], M[:3, 3]
                if S._key(R * cref_M + t, space, shifts) == S._key(sp.Matrix(ce), space, shifts):
                    g = R; break
            if g is None:
                res["problems"].append(f"{lab}: first center is not an image of the representative center")
                continue
            sub = image(g.inv(), entry_subspace(S._parse(v)[0]))
            for ci, cl in enumerate(classes):
                if sub in cl:
                    hit[ci].append(lab); break
            else:
                unknown.append((lab, form(sub)))
            # label multiplicity = number of listed bonds times the centering size
            nb = int(re.match(r"\d+", lab).group())
            if nb != len(m) * k:
                res["problems"].append(f"{lab}: bond multiplicity {nb} != {len(m)} listed bonds x {k}")
        missing = [form(min(cl, key=lambda s: (len(s), str(s)))) for ci, cl in enumerate(classes) if ci not in hit]
        dups = {ci: labs for ci, labs in hit.items() if len(labs) > 1}
        res["per_W"][W] = dict(n_strata=len(classes), n_file=len(ents), missing=missing,
                               dups=[(classes and labs) for labs in dups.values()], unknown=unknown, H_order=len(H))
    return res


if __name__ == "__main__":
    names = sys.argv[1:] or list(data)
    out = {}
    for n in names:
        out[n] = check_group(n)
    pickle.dump(out, open("complete_result.pkl", "wb"))
    nW = sum(len(r["per_W"]) for r in out.values())
    print(f"groups checked {len(out)}, ordinary Wyckoff positions {nW}")
