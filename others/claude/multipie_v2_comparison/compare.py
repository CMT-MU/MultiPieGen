"""
Compare the multipole-keyed data of MultiPieGen with MultiPie Ver.2 (verification only, not part of the package).

    python others/claude/multipie_v2_comparison/compare.py [section ...]

Sections: harmonics, mp_cache, multipolar, cluster_basis, cluster_samb, atomic, repmat, active (default: all, in this
order). Each section runs in its own process, one after another (no parallel workers), so that only one large binary
is in memory at a time. MultiPie Ver.2's harmonics_multipole.pkl alone needs ~8 GB and ~4 min to load; the "mp_cache"
section loads it once and writes the l <= 3 part, evaluated numerically, to cache/ (not in git), which "multipolar"
reads.

MultiPie Ver.2 is read from ../MultiPie (others/tools/binary_data, multipie/binary_data), or from $MULTIPIE_DIR.
Results are written to result/<section>.json.
"""

import collections
import json
import os
import pickle
import subprocess
import sys
import time

import numpy as np
import sympy as sp
from multipie.util.util_binary import BinaryManager

HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.environ.get("MULTIPIE_DIR", os.path.join(HERE, "../../../../MultiPie"))
TOOLS = os.path.join(MP, "others/tools/binary_data")
PKG = os.path.join(MP, "multipie/binary_data")
RESULT = os.path.join(HERE, "result")
CACHE = os.path.join(HERE, "cache")
SECTIONS = ["harmonics", "mp_cache", "multipolar", "cluster_basis", "cluster_samb", "atomic", "repmat", "active"]

_x, _y, _z = sp.symbols("x y z", real=True)
_P = [(0.31, 0.72, 1.13), (-0.57, 0.23, 0.41), (0.93, -0.38, -0.66)]  # sample points of x, y, z.


# ==================================================
def ex_val(e):
    """Values of a cartesian expression (in x, y, z) at the sample points."""
    e = sp.S(e)
    syms = {s.name: s for s in e.free_symbols}
    f = sp.lambdify([syms.get("x", _x), syms.get("y", _y), syms.get("z", _z)], e, "numpy")
    return np.array([complex(f(*p)) for p in _P])


def arr_val(a):
    """Values of an array of expressions, shape a.shape + (number of sample points,)."""
    a = np.asarray(a, dtype=object)
    return np.array([ex_val(e) for e in a.ravel()]).reshape(a.shape + (len(_P),))


def classify(a, b):
    if a.shape != b.shape:
        return "different"
    if np.allclose(a, b, atol=1e-9):
        return "same"
    if np.allclose(a, -b, atol=1e-9):
        return "sign"
    return "different"


def save(section, data):
    os.makedirs(RESULT, exist_ok=True)
    with open(os.path.join(RESULT, f"{section}.json"), "w") as f:
        json.dump(data, f, indent=1, default=str, ensure_ascii=False)


def _info():
    from multipiegen.util.util_info import load_info

    return load_info()


# ==================================================
def harmonics():
    """Harmonics of the 47 PG up to rank 11 (MultiPie Ver.2 harmonics.pkl)."""
    from multipiegen.util.util_harmonics_multipole import pg_harmonics

    info = _info()
    ref_all = BinaryManager("harmonics", topdir=TOOLS)
    total, groups = collections.Counter(), {}
    for gid in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        pg = info["tag"][gid]
        ref = {(k[0], k[1], k[2], k[3]): v for k, v in ref_all[gid]["harmonics"].items()}
        ours = pg_harmonics(pg, 11)
        res, diff = collections.Counter(), []
        for key in sorted(set(ref) | set(ours), key=str):
            if key not in ref or key not in ours:
                c = "only ours" if key not in ref else "only MultiPie"
            else:
                c = classify(arr_val(ours[key][0]), arr_val(ref[key][0]))
            res[c] += 1
            if c != "same":
                diff.append((c, list(key)))
        total.update(res)
        if diff:
            groups[pg] = {"count": dict(res), "diff": diff}
    save("harmonics", {"total": dict(total), "groups": groups})
    print("harmonics:", dict(total))


# ==================================================
def mp_cache():
    """Load MultiPie Ver.2 harmonics_multipole.pkl once, and cache the numerical values of l <= 3."""
    os.makedirs(CACHE, exist_ok=True)
    ref = BinaryManager("harmonics_multipole", topdir=TOOLS)
    data = {}
    for gid in ref.keys():
        data[gid] = [(tuple(k), arr_val(v[0]), arr_val(v[1])) for k, v in ref[gid].items() if k[1] <= 3]
    with open(os.path.join(CACHE, "harmonics_multipole_l3.pkl"), "wb") as f:
        pickle.dump(data, f)
    print("mp_cache:", sum(len(v) for v in data.values()), "entries")


def _tess_matrix(s, crystal):
    """T[a, b] = <tesseral b | orbital a>: coefficients of the orbitals (cubic/hexagonal) -> tesseral basis."""
    from multipiegen.data.data_basis import atomic_orbital_basis, cubic_orbital_basis, hexagonal_orbital_basis

    bf = _info()["harmonics"]["basis_function"]
    orb = (cubic_orbital_basis if crystal == "cubic" else hexagonal_orbital_basis)[s]
    Uo = np.array([np.array(bf[n][1], dtype=complex) for n in orb])
    Ut = np.array([np.array(bf[n][1], dtype=complex) for n in atomic_orbital_basis[s]])
    return Uo @ Ut.conj().T


def multipolar():
    """Multipolar harmonics (s = 1..3, l <= 3) vs the cache of MultiPie Ver.2 harmonics_multipole.pkl."""
    from multipiegen.data.data_harmonics_mapping import data_harmonics_mapping
    from multipiegen.util.util_harmonics_multipole import load_multipolar_harmonics

    info = _info()
    with open(os.path.join(CACHE, "harmonics_multipole_l3.pkl"), "rb") as f:
        cache = pickle.load(f)
    ours_all = load_multipolar_harmonics()
    total, groups, internal_mismatch = collections.Counter(), {}, collections.Counter()
    for gid in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        pg = info["tag"][gid]
        crystal = "cubic" if data_harmonics_mapping[pg]["top"] == "Oh" else "hexagonal"
        rg = collections.defaultdict(list)
        for k, ex, mh in cache.get(gid, []):
            rg[(k[0], k[1], k[2], k[5], k[6], k[7])].append((k, ex, mh))
        og = collections.defaultdict(list)
        for k, (ex, mh, mex) in ours_all[gid].items():
            T = _tess_matrix(k.s, crystal)
            m = np.einsum("iap,ab->ibp", arr_val(mh), T)
            og[(k.X, k.l, k.Gamma, k.s, k.k, k.x)].append((k, arr_val(ex), m, arr_val(mh)))
        res, diff = collections.Counter(), []
        for g in sorted(set(rg) | set(og), key=str):
            r, o = rg.get(g, []), og.get(g, [])
            if len(r) != len(o):
                res["(X,l,Gamma,s,k,x) with different number of multiplets"] += 1
                diff.append(("count", list(g), len(o), len(r)))
            for k, ex, m, mh in o:
                c = "not found"
                for kr, exr, mr in r:
                    cc = classify(m, mr)
                    if cc != "different":
                        c = "same" if cc == "same" else "sign"
                        res["n,p same" if (k.n, k.p) == (kr[3], kr[4]) else "n,p renumbered"] += 1
                        break
                if c == "not found":  # same harmonics (ex) but internal-basis coefficients differ ?
                    for kr, exr, mr in r:
                        if classify(ex, exr) != "different":
                            internal_mismatch[(crystal, k.s)] += 1
                            c = "same harmonics, different internal part"
                            break
                res[c] += 1
                if c not in ("same", "sign"):
                    diff.append((c, list(k)))
        total.update(res)
        groups[pg] = {"count": dict(res), "diff": diff[:30]}
    save(
        "multipolar",
        {
            "total": dict(total),
            "internal mismatch (crystal, s)": {str(k): v for k, v in internal_mismatch.items()},
            "groups": groups,
        },
    )
    print("multipolar:", dict(total), dict(internal_mismatch))


# ==================================================
def cluster_basis():
    """Cluster bases of Oh/D6h vs MultiPie Ver.2 root-cluster harmonics (raw harmonics, not orthonormalized)."""
    from multipiegen.util.util_bond_cluster_basis import load_bond_cluster_basis
    from multipiegen.util.util_cluster_basis import load_cluster_basis
    from multipiegen.util.util_cluster_vector_basis import load_cluster_vector_basis

    rc = BinaryManager("harmonics_root_cluster", topdir=TOOLS)
    cb, vb, bb = load_cluster_basis(), load_cluster_vector_basis(), load_bond_cluster_basis()
    out = {}
    for top, gid in [("Oh", "PG:32"), ("D6h", "PG:27")]:
        mono = {(k[0], k[1], k[2], k[3]) for k in rc[gid]["monopole"]}
        dip = {(k[0], k[1], k[2], k[3], k[6]) for k in rc[gid]["dipole"]}
        sc = {(k.X, k.l, k.Gamma, k.n) for k in cb[top]["basis"]}
        vq = {(k.X, k.l, k.Gamma, k.n, k.k) for k in vb[top]["basis"] if k.x == "q"}
        out[top] = {
            "MultiPie monopole multiplets (rank <= 11, raw)": len(mono),
            "MultiPie monopole points": int(np.shape(next(iter(rc[gid]["monopole"].values()))[0])[1]),
            "scalar basis multiplets (kept)": len(sc),
            "scalar basis keys not in MultiPie monopole": sorted(map(str, sc - mono)),
            "MultiPie dipole multiplets (rank <= 11, raw, x=q)": len(dip),
            "vector basis multiplets (x=q, kept)": len(vq),
            "vector basis multiplets (x=g, kept)": len([k for k in vb[top]["basis"] if k.x == "g"]),
            "vector basis keys (x=q) not in MultiPie dipole": sorted(map(str, vq - dip)),
            "bond basis multiplets (x=q)": len([k for k in bb[top]["basis"] if k.x == "q"]),
            "bond basis multiplets (x=g)": len([k for k in bb[top]["basis"] if k.x == "g"]),
        }
    save("cluster_basis", out)
    print("cluster_basis:", json.dumps(out, default=str)[:600])


# ==================================================
def cluster_samb():
    """Cluster SAMB of the real PG and SG vs MultiPie Ver.2 cluster_samb.pkl."""
    from multipiegen.util.util_cluster_samb import load_cluster_samb

    ours, ref = load_cluster_samb(), BinaryManager("cluster_samb", topdir=TOOLS)
    fl = lambda m: np.array(m, dtype=float)
    proj = lambda r: np.array(r).T @ np.linalg.pinv(np.array(r).T)
    total, bond_a = collections.Counter(), []
    pg_bond_a = collections.Counter()
    for gid in ours.keys():
        for kind in ["site", "bond_s", "bond_a"]:
            if list(ours[gid][kind]) != list(ref[gid][kind]):
                total[f"{kind}: Wyckoff list differs"] += 1
            for w, samb in ours[gid][kind].items():
                rr = ref[gid][kind][w]
                ga, gb = collections.defaultdict(list), collections.defaultdict(list)
                for k, (s, e) in rr.items():
                    ga[k[2]] += list(fl(s))
                for k, (s, e) in samb.items():
                    gb[k.Gamma] += list(fl(s))
                ok = set(ga) == set(gb) and all(np.allclose(proj(ga[G]), proj(gb[G]), atol=1e-8) for G in gb)
                total[f"{kind}: span per irrep {'same' if ok else 'differs'}"] += 1
                ka = collections.Counter((k[0], k[1], k[2], k[7]) for k in rr)
                kb = collections.Counter((k.X, k.l, k.Gamma, k.x) for k in samb)
                total[f"{kind}: keys (X,l,Gamma,x) {'same' if ka == kb else 'differ'}"] += 1
                kn = collections.Counter((k[0], k[1], k[2], k[3], k[4], k[7]) for k in rr)
                kbn = collections.Counter((k.X, k.l, k.Gamma, k.n, k.p, k.x) for k in samb)
                total[f"{kind}: full keys {'same' if kn == kbn else 'differ'}"] += 1
                if ka != kb and kind == "bond_a":
                    if gid.startswith("PG"):
                        pg_bond_a[gid] += 1
                    if len(bond_a) < 40:
                        bond_a.append((gid, w, sorted(map(list, (kb - ka).elements())), sorted(map(list, (ka - kb).elements()))))
    save("cluster_samb", {"total": dict(total), "bond_a PG": dict(pg_bond_a), "bond_a examples (ours, MultiPie)": bond_a})
    print("cluster_samb:", dict(total))


# ==================================================
def atomic():
    """Atomic multipoles vs MultiPie Ver.2 atomic_multipole.pkl and atomic_multipole_group.pkl."""
    from multipiegen.util.util_atomic_multipole import load_atomic_multipole, load_atomic_multipole_group

    info = _info()
    a, b = BinaryManager("atomic_multipole", topdir=TOOLS), load_atomic_multipole()
    keys_same = all([tuple(k) for k in a[bt][blk]] == [tuple(k) for k in b[bt][blk]] for bt in a.keys() for blk in a[bt])
    vals_same = all(
        np.allclose(np.array(a[bt][blk][k], dtype=complex), np.array(b[bt][blk][k], dtype=complex))
        for bt in a.keys()
        for blk in a[bt]
        for k in b[bt][blk]
    )
    del a, b
    a, b = BinaryManager("atomic_multipole_group", topdir=TOOLS), load_atomic_multipole_group()
    total, signs, nren = collections.Counter(), set(), collections.Counter()
    for gid in info["id_set"]["PG"]["all"]:
        for bt in ["jml", "lgs", "lg"]:
            for blk, d in b[gid][bt].items():
                A = collections.defaultdict(list)
                for k, (m, e) in a[gid][bt][blk].items():
                    A[(k[0], k[1], k[2], k[5], k[6])].append((k, np.array(m, dtype=complex)))
                for k, (m, e) in d.items():
                    m = np.array(m, dtype=complex)
                    c = "not found"
                    for kr, mr in A.get((k.X, k.l, k.Gamma, k.s, k.k), []):
                        cc = classify(m, mr)
                        if cc != "different":
                            c = cc
                            nren["n same" if kr[3] == k.n else "n differs"] += 1
                            break
                    total[c] += 1
                    if c == "sign":
                        signs.add(f"{info['tag'][gid]} l={k.l} {k.Gamma}")
    out = {
        "atomic_multipole: keys and order same": keys_same,
        "atomic_multipole: values same": vals_same,
        "atomic_multipole_group (37 real PG): matrices": dict(total),
        "atomic_multipole_group: n": dict(nren),
        "atomic_multipole_group: sign differences (PG, l, irrep)": sorted(signs),
    }
    save("atomic", out)
    print("atomic:", {k: v for k, v in out.items() if "sign differences" not in k})


# ==================================================
def repmat():
    """Representation matrices vs MultiPie Ver.2 representation_matrix.pkl."""
    from multipiegen.util.util_representation_matrix import representation_matrix

    info = _info()
    rm = BinaryManager("representation_matrix", topdir=TOOLS)
    total, diff = collections.Counter(), []
    for gid in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        pg = info["tag"][gid]
        o, r = representation_matrix(pg), rm[gid]
        if list(o["tag"]) != list(r["tag"]):
            total["SO order differs"] += 1
        for ir, m in o["matrix"].items():
            c = (
                classify(np.array(m, dtype=complex), np.array(r["matrix"][ir], dtype=complex))
                if ir in r["matrix"]
                else "irrep only ours"
            )
            total[c] += 1
            if c != "same":
                diff.append((pg, ir, c))
        total["irrep only MultiPie"] += len([ir for ir in r["matrix"].keys() if ir not in o["matrix"]])
    save("repmat", {"total": dict(total), "diff": diff})
    print("repmat:", dict(total), diff[:10])


# old (MultiPie Ver.2) -> new orbital names (data_basis), in the active multipoles "X + orbital".
_RENAME = {
    "f2": "fc",
    "f1": "fs",
    "f3x": "fx",
    "f3y": "fy",
    "ga": "gc",
    "gb": "gs",
    "gc": "gc4",
    "i": "i3",
    "l": "l3",
    "i0": "ixy",
    "j": "j3",
}


def _rename(tag):
    return tag[0] + _RENAME.get(tag[1:], tag[1:])


# ==================================================
def active():
    """Active multipoles of MPG (group.pkl) vs MultiPie Ver.2 MPG binaries."""
    from multipiegen.util.util_group import load_group

    info, group = _info(), load_group()
    total, diff = collections.Counter(), []
    for gid in info["id_set"]["MPG"]["all"]:
        r = BinaryManager(gid, topdir=os.path.join(PKG, "MPG"))["active_multipole"]
        o = group[gid].get("active_multipole")
        if list(o) == list(r):
            c = "same"
        elif [_rename(t) for t in r] == list(o):
            c = "same after renaming the orbitals"
        else:
            c = "different"
        total[c] += 1
        if c == "different":
            diff.append((gid, list(o), list(r)))
    save("active", {"total": dict(total), "diff": diff})
    print("active:", dict(total), diff[:5])


# ==================================================
if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) == 1 and args[0] in SECTIONS:
        globals()[args[0]]()
    else:
        for sec in args or SECTIONS:
            t0 = time.time()
            subprocess.run([sys.executable, os.path.abspath(__file__), sec], check=True)
            print(f"  ({sec}: {time.time() - t0:.0f} s)", flush=True)
