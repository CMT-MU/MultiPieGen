"""Site symmetry of Wyckoff positions with vectors.

Assumptions
  * Symmetry operations use the column-vector convention r' = R r + t. A 4x4 matrix is
    [[R, t], [0, 1]]; a 3x3 matrix is a point-group operation (t = 0).
    Pass a list, or an ndarray (N,4,4) / (N,3,3) of sympy entries.
  * Positions are fractional coordinates with free parameters x, y, z. Vector components are
    X, Y, Z, in the same (crystal-axes) basis as the positions.
  * Operations are numbered in the order given, starting at 1 (the same indices as in the data).
  * For a non-primitive lattice (A, B, C, I, F, R) the list must contain the centering translations
    too (e.g. 8 operations for C2/m); use expand_centering() to build it from the point operations,
    and list every conventional-cell copy of a site in the entry. Alternatively pass the primitive
    operations (without centering) together with `centering="C"` etc.: positions that differ by a
    centering vector are then the same site, and the entries list the primitive-cell bonds only.
  * The third element of an entry (the indices) is optional; it is only used for index_matches_given.
  * With 4x4 matrices (space groups) positions are compared modulo integer lattice translations.
    With 3x3 matrices (point groups) they are compared exactly. `lattice_mod` overrides this.
  * Positions must be exact expressions (integers / rationals). Floats are not supported.

Site symmetry = the set of operations that fix the representative position (the first entry
of the list) for all values of x, y, z.
The oriented symbol (ITA style: ..m, m2m, 4/mm.m, ...) is stored in SiteSymmetry.symbol. The
orientation refers to the representative position (the first ITA coordinate triplet is assumed).
Following the symmetry directions of the crystal system ([001], [100]/[010], [1-10]/[110], ...),
the symmetry elements along each direction are written in order, and directions without any
symmetry element are filled with ".". The crystal must be given in the standard axes of the
conventional cell (monoclinic: any unique axis; trigonal: hexagonal or rhombohedral axes).
The site symmetry itself ignores the vectors. When vectors are passed, the bond symmetry is
computed as well (fields bond_*): the operations that fix the representative position AND map
the representative vector V onto +V or -V. Only the rotation part R acts on a vector (partial
translations are excluded), and the comparison holds identically in the vector parameters
X, Y, Z. The direction of the bond is ignored, so an operation that maps V to -V belongs to the
bond symmetry; those operations are listed in bond_flips (the same meaning as the minus sign in
[[1, -2]] in the data). With vectors, an operation maps the representative bond onto entry j only
if both the position and the vector match, so several entries may share one position.
"""

from __future__ import annotations
import os

import re
from dataclasses import dataclass

import sympy as sp

# ==================================================
header = """
Wyckoff bond for all space and point groups (fractional coordinate, conventional cell, no plus set).
- wyckoff_bond_data.
    { group: { Wyckoff: (vector, center, mapping, bond_symmetry) } }
"""

x, y, z, X, Y, Z = sp.symbols("x y z X Y Z")
_NAMES = {"x": x, "y": y, "z": z, "X": X, "Y": Y, "Z": Z}
_POS = (x, y, z)
_ZERO = {s: 0 for s in _POS}

# (det, trace) of the rotation part -> operation type (holds for integer matrices in any lattice basis)
_TYPE = {
    (1, 3): "1",
    (1, -1): "2",
    (1, 0): "3",
    (1, 1): "4",
    (1, 2): "6",
    (-1, -3): "-1",
    (-1, 1): "m",
    (-1, 0): "-3",
    (-1, -1): "-4",
    (-1, -2): "-6",
}
_ORDER = ["1", "2", "3", "4", "6", "-1", "m", "-3", "-4", "-6"]

# The 32 crystallographic point groups: number of operations of each type. This tuple of
# counts determines each of the 32 groups uniquely.
_GROUPS = {
    "1": {"1": 1},
    "-1": {"1": 1, "-1": 1},
    "2": {"1": 1, "2": 1},
    "m": {"1": 1, "m": 1},
    "2/m": {"1": 1, "2": 1, "-1": 1, "m": 1},
    "222": {"1": 1, "2": 3},
    "mm2": {"1": 1, "2": 1, "m": 2},
    "mmm": {"1": 1, "2": 3, "-1": 1, "m": 3},
    "4": {"1": 1, "2": 1, "4": 2},
    "-4": {"1": 1, "2": 1, "-4": 2},
    "4/m": {"1": 1, "2": 1, "4": 2, "-1": 1, "m": 1, "-4": 2},
    "422": {"1": 1, "2": 5, "4": 2},
    "4mm": {"1": 1, "2": 1, "4": 2, "m": 4},
    "-42m": {"1": 1, "2": 3, "m": 2, "-4": 2},
    "4/mmm": {"1": 1, "2": 5, "4": 2, "-1": 1, "m": 5, "-4": 2},
    "3": {"1": 1, "3": 2},
    "-3": {"1": 1, "3": 2, "-1": 1, "-3": 2},
    "32": {"1": 1, "2": 3, "3": 2},
    "3m": {"1": 1, "3": 2, "m": 3},
    "-3m": {"1": 1, "2": 3, "3": 2, "-1": 1, "m": 3, "-3": 2},
    "6": {"1": 1, "2": 1, "3": 2, "6": 2},
    "-6": {"1": 1, "3": 2, "m": 1, "-6": 2},
    "6/m": {"1": 1, "2": 1, "3": 2, "6": 2, "-1": 1, "m": 1, "-3": 2, "-6": 2},
    "622": {"1": 1, "2": 7, "3": 2, "6": 2},
    "6mm": {"1": 1, "2": 1, "3": 2, "6": 2, "m": 6},
    "-6m2": {"1": 1, "2": 3, "3": 2, "m": 4, "-6": 2},
    "6/mmm": {"1": 1, "2": 7, "3": 2, "6": 2, "-1": 1, "m": 7, "-3": 2, "-6": 2},
    "23": {"1": 1, "2": 3, "3": 8},
    "m-3": {"1": 1, "2": 3, "3": 8, "-1": 1, "m": 3, "-3": 8},
    "432": {"1": 1, "2": 9, "3": 8, "4": 6},
    "-43m": {"1": 1, "2": 3, "3": 8, "m": 6, "-4": 6},
    "m-3m": {"1": 1, "2": 9, "3": 8, "4": 6, "-1": 1, "m": 9, "-3": 8, "-4": 6},
}
_SIGNATURE = {tuple(d.get(t, 0) for t in _ORDER): name for name, d in _GROUPS.items()}
assert len(_SIGNATURE) == len(_GROUPS) == 32


# ---------------------------------------------------------------- input handling
def _parse(s):
    """sympify a string, accepting implicit multiplication such as '2X', '-2x+1/2', '1/2x', '2(x+y)'."""
    s = re.sub(r"(?<![A-Za-z_])(\d+)\s*(?=[xyzXYZ(])", r"\1*", s)  # 2X -> 2*X, 2(x+y) -> 2*(x+y)
    s = re.sub(r"\)\s*(?=[xyzXYZ(])", ")*", s)  # (x+y)z -> (x+y)*z
    return sp.sympify(s, locals=_NAMES)


def _vectors(obj):
    """Turn '[[x,y,z],...]' or a nested list into a list of 3-component column vectors."""
    if isinstance(obj, str):
        obj = _parse(obj)
    vs = [sp.Matrix(_parse(v) if isinstance(v, str) else v) for v in obj]
    if any(v.shape != (3, 1) for v in vs):
        raise ValueError("each position / vector must have exactly 3 components")
    return vs


def _normalize_ops(ops):
    """-> ([(R, t), ...], whether any 4x4 matrix was given)"""
    pairs, has_t = [], False
    for op in ops:
        M = sp.Matrix(op)
        if M.shape == (3, 3):
            R, t = M, sp.zeros(3, 1)
        elif M.shape == (4, 4):
            if list(M[3, :]) != [0, 0, 0, 1]:
                raise ValueError("the last row of a 4x4 matrix must be [0, 0, 0, 1]")
            R, t, has_t = M[:3, :3], M[:3, 3], True
        else:
            raise ValueError(f"a symmetry operation must be 3x3 or 4x4, got {M.shape}")
        pairs.append((R, t))
    return pairs, has_t


def _key(v, lattice_mod, shifts=None):
    """Reduce a position to (linear part in x,y,z, constant part); the constant is taken mod 1 if lattice_mod.
    shifts: the translations [0, c1, ...] of a centered lattice; positions that differ by one of them
    (mod 1) get the same key."""
    if shifts and len(shifts) > 1:
        v = sp.Matrix(v)
        return min((_key(v + s, lattice_mod) for s in shifts), key=str)
    key = []
    for e in v:
        e = sp.expand(e)
        if not e.free_symbols <= set(_POS):
            raise ValueError(f"position contains a symbol other than x, y, z: {e}")
        c = e.subs(_ZERO)
        key.append((sp.expand(e - c), c % 1 if lattice_mod else c))
    return tuple(key)


def _centering_shifts(centering):
    """None -> None. A lattice letter or a list of vectors -> [0, c1, c2, ...] (sympy column vectors, mod 1)."""
    if centering is None:
        return None
    vecs = CENTERING[centering] if isinstance(centering, str) else centering
    shifts = [sp.zeros(3, 1)]
    for v in vecs:
        c = sp.Matrix([sp.Rational(e) % 1 for e in v])
        if all(c != s for s in shifts):
            shifts.append(c)
    return shifts


def _use_centering(pairs, centering):
    """-> (shifts or None, whether it is an R lattice in hexagonal axes).
    Refuses an operation list that already contains the centering translations."""
    shifts = _centering_shifts(centering)
    if shifts is None or len(shifts) == 1:
        return None, False
    if any(R == sp.eye(3) and any(c % 1 for c in t) for R, t in pairs):
        raise ValueError(
            "`centering` was given, but the operations already contain centering translations: "
            "pass either the full operations, or the primitive operations together with `centering`"
        )
    third = (sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(2, 3))
    return shifts, any(tuple(sorted(s)) == third for s in shifts)


def _sign(v, w):
    """+1 if v == w, -1 if v == -w, 0 otherwise."""
    if all(sp.expand(a - b) == 0 for a, b in zip(v, w)):
        return 1
    if all(sp.expand(a + b) == 0 for a, b in zip(v, w)):
        return -1
    return 0


def _op_type(R):
    return _TYPE.get((R.det(), R.trace()))


def _point_group(Rs):
    """International symbol of the point group formed by the rotation parts (None if not identified)."""
    uniq = list({tuple(R): R for R in Rs}.values())  # drop duplicates caused by pure lattice translations
    types = [_op_type(R) for R in uniq]
    if None in types:
        return None
    return _SIGNATURE.get(tuple(types.count(t) for t in _ORDER))


# ---------------------------------------------------------------- oriented symbol
from functools import reduce
from math import gcd


def _prim(v):
    """Scale a rational vector to a primitive integer vector with a positive first nonzero entry (+/- identified)."""
    v = [sp.Rational(e) for e in v]
    den = reduce(lambda a, b: a * b // gcd(a, b), [e.q for e in v], 1)
    w = [int(e * den) for e in v]
    g = reduce(gcd, [abs(a) for a in w]) or 1
    w = [a // g for a in w]
    for a in w:
        if a:
            if a < 0:
                w = [-b for b in w]
            break
    return tuple(w)


def _axis(R, typ):
    """Direction of the symmetry element: the axis for rotations / rotoinversions, the normal for mirrors. None for 1 and -1."""
    I3 = sp.eye(3)
    if typ in ("2", "3", "4", "6"):
        A = R - I3
    elif typ == "m":
        A = R + I3
    elif typ in ("-3", "-4", "-6"):
        A = -R - I3
    else:
        return None
    ns = A.nullspace()
    return _prim(ns[0]) if len(ns) == 1 else None


# Symmetry directions of each crystal system (ITA Table 2.1.3.1). Each element is a set of
# equivalent directions, listed in the order in which they appear in the symbol.
_P = _prim
_SETS = {
    # Third set: list [110] before [1-10] to match the ITA tables (m.2m, 2mm, ...)
    "tetragonal": [[_P((0, 0, 1))], [_P((1, 0, 0)), _P((0, 1, 0))], [_P((1, 1, 0)), _P((1, -1, 0))]],
    "hexagonal": [
        [_P((0, 0, 1))],
        [_P((1, 0, 0)), _P((0, 1, 0)), _P((-1, -1, 0))],
        [_P((1, -1, 0)), _P((1, 2, 0)), _P((-2, -1, 0))],
    ],
    "rhombohedral": [[_P((1, 1, 1))], [_P((1, -1, 0)), _P((0, 1, -1)), _P((-1, 0, 1))]],
    "cubic": [
        [_P((1, 0, 0)), _P((0, 1, 0)), _P((0, 0, 1))],
        [_P((1, 1, 1)), _P((1, -1, 1)), _P((-1, 1, 1)), _P((-1, -1, 1))],
        [_P((1, -1, 0)), _P((1, 1, 0)), _P((0, 1, -1)), _P((0, 1, 1)), _P((-1, 0, 1)), _P((1, 0, 1))],
    ],
    "orthorhombic": [[_P((1, 0, 0))], [_P((0, 1, 0))], [_P((0, 0, 1))]],
    "triclinic": [],
}
_SYSTEM = {}
for _s, _names in {
    "triclinic": ["1", "-1"],
    "monoclinic": ["2", "m", "2/m"],
    "orthorhombic": ["222", "mm2", "mmm"],
    "tetragonal": ["4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm"],
    "trigonal": ["3", "-3", "32", "3m", "-3m"],
    "hexagonal": ["6", "-6", "6/m", "622", "6mm", "-6m2", "6/mmm"],
    "cubic": ["23", "m-3", "432", "-43m", "m-3m"],
}.items():
    for _n in _names:
        _SYSTEM[_n] = _s


def _direction_sets(group_Rs, r_centered=False):
    """Determine the crystal system and its symmetry directions from the rotation parts of the whole space group.
    r_centered: R lattice in hexagonal axes (ITA writes these with the two-position symbols of the
    rhombohedral system: 3m, -3m, .2/m, ...)."""
    pg = _point_group(group_Rs)
    if pg is None:
        raise ValueError("cannot identify the point group (the operations do not form a group, or the basis is not integral)")
    system = _SYSTEM[pg]
    info = [(_op_type(R), R) for R in group_Rs]
    if system == "monoclinic":
        for typ, R in info:
            if typ in ("2", "m"):
                return system, [[_axis(R, typ)]]
    if system == "trigonal":
        ax = {_axis(R, t) for t, R in info if t in ("3", "-3")}
        if ax == {_P((0, 0, 1))}:
            return "trigonal(hex)", (_SETS["hexagonal"][:2] if r_centered else _SETS["hexagonal"])
        if ax == {_P((1, 1, 1))}:
            return "trigonal(rhomb)", _SETS["rhombohedral"]
        raise ValueError(f"the 3-fold axis of the trigonal group is neither [001] nor [111]: {ax}")
    if system in ("tetragonal", "hexagonal"):
        ax = {_axis(R, t) for t, R in info if t in ("4", "-4", "6", "-6")}
        if ax != {_P((0, 0, 1))}:
            raise ValueError(
                f"the principal axis of the {system} group is not [001] (non-standard orientations are not supported): {ax}"
            )
    return system, _SETS[system]


def _components(Hinfo, n_H, has_inv, d):
    """Combine the symmetry elements along direction d (rotation, rotoinversion, mirror) into one HM symbol."""
    rot = max([int(t) for t, ax in Hinfo if ax == d and t in ("2", "3", "4", "6")], default=1)
    rotinv = {t for t, ax in Hinfo if ax == d and t in ("-4", "-6")}
    mir = any(t == "m" and ax == d for t, ax in Hinfo)
    if rot == 6:
        return "6/m" if mir else "6"
    if rot == 4:
        return "4/m" if mir else "4"
    if rot == 3:
        return "-6" if mir else ("-3" if has_inv else "3")
    if rot == 2:
        if "-4" in rotinv:
            return "-4"
        if mir:
            return "2/m" if n_H == 4 else "m"
        return "2"
    return "m" if mir else ""


def _oriented_symbol(group_Rs, H_Rs, r_centered=False):
    system, sets = _direction_sets(group_Rs, r_centered)
    n_H = len(H_Rs)
    Hinfo = [(_op_type(R), _axis(R, _op_type(R))) for R in H_Rs]
    has_inv = any(t == "-1" for t, _ in Hinfo)
    orbit = {}

    def orb(d):
        if d not in orbit:
            orbit[d] = {_prim(R * sp.Matrix(d)) for R in H_Rs}
        return orbit[d]

    tokens = []
    for k, S in enumerate(sets):
        carriers = []
        for d in S:
            c = _components(Hinfo, n_H, has_inv, d)
            if not c or any(c == c0 and d in orb(d0) for d0, c0 in carriers):
                continue
            carriers.append((d, c))
        cs = [c for _, c in carriers]
        if system == "tetragonal" and k > 0:
            cs.sort(key=lambda c: c != "2")  # tetragonal 2nd/3rd positions: rotation axis first (2m, never m2)
        if system == "cubic" and k == 0:
            if sorted(cs) == ["2", "m", "m"]:
                cs = ["m", "m", "2"]  # ITA writes mm2 (not 2mm) in the first position of cubic groups
            elif cs == ["4/m"] and len(orb(carriers[0][0])) == 3:
                cs = ["m"]  # if the three <100> are equivalent, write "m" as in m-3m (not 4/m)
        tokens.append("".join(cs) or ".")
    if all(t == "." for t in tokens):
        return "-1" if has_inv else "1"
    return "".join(tokens)


# ---------------------------------------------------------------- result
@dataclass
class SiteSymmetry:
    label: str | None
    representative: sp.Matrix  # representative position (first in the list)
    stabilizer: list  # numbers (1-based) of the operations that fix the representative position
    order: int  # order of the site-symmetry group
    point_group: str | None  # international symbol of the site-symmetry group (orientation not distinguished)
    op_types: list  # type of each stabilizer operation ('1','2','m','-1',...)
    orbit_ops: list  # per position: numbers of the operations mapping the representative onto it
    vector_sign: list | None  # same shape as orbit_ops; +1: +V, -1: -V, 0: neither
    unmatched: list  # operations whose image matched none of the given positions (list not closed)
    multiplicity_ok: bool  # len(positions) * order == total number of operations
    index_matches_given: bool | None = None
    symbol: str | None = None  # oriented site-symmetry symbol (ITA style: ..m, m2m, 4/mm.m, ...)
    # Bond symmetry (only when vectors are given): the subgroup of the site symmetry that also maps the
    # representative vector onto itself, ignoring its direction (R V == +-V). Rotation parts only.
    bond_stabilizer: list | None = None  # numbers (1-based) of these operations
    bond_order: int | None = None
    bond_flips: list | None = None  # those among them that reverse the vector (R V == -V)
    bond_op_types: list | None = None
    bond_point_group: str | None = None  # international symbol, orientation not distinguished
    bond_symbol: str | None = None  # oriented symbol, same conventions as `symbol`
    bond_multiplicity_ok: bool | None = None  # len(positions) * bond_order == total number of operations

    def signed_orbit_ops(self):
        """Same format as the third element of the data (operations that flip the vector get a minus sign)."""
        if self.vector_sign is None:
            raise ValueError("no vectors were given")
        out = []
        for ops, signs in zip(self.orbit_ops, self.vector_sign):
            if 0 in signs:
                raise ValueError("some operation maps the representative vector to neither +V nor -V")
            out.append(sorted((i * s for i, s in zip(ops, signs)), key=abs))
        return out

    def __str__(self):
        head = f"{self.label}: " if self.label else ""
        s = (
            f"{head}site symmetry {self.symbol or self.point_group} (order {self.order}), "
            f"ops {self.stabilizer}, multiplicity {len(self.orbit_ops)}"
        )
        if self.bond_stabilizer is not None:
            s += (
                f", bond symmetry {self.bond_symbol or self.bond_point_group} (order {self.bond_order}), "
                f"ops {self.bond_stabilizer}, flipping {self.bond_flips}"
            )
            if not self.bond_multiplicity_ok:
                s += " [len(positions)*bond_order != |G|]"
        if self.vector_sign is not None and all(0 not in v for v in self.vector_sign):
            s += f", signed ops {self.signed_orbit_ops()}"
        if self.index_matches_given is not None:
            s += f", matches given: {self.index_matches_given}"
        if self.unmatched:
            s += f", NOT CLOSED (unmatched ops {self.unmatched})"
        if not self.multiplicity_ok:
            s += ", multiplicity*order != |G|"
        return s


# ---------------------------------------------------------------- main routine
def site_symmetry(ops, positions, vectors=None, *, lattice_mod=None, label=None, oriented_errors_ok=False, centering=None):
    """Determine the site symmetry of a Wyckoff position (`positions`).

    ops       : symmetry operations (3x3 or 4x4, sympy entries). Their order gives the numbers (1-based).
    positions : '[[x,y,z],...]' or a nested list. The first entry is the representative position.
    vectors   : the same number of vectors '[[X,Y,Z],...]' (optional).
    centering : lattice letter ("A", "B", "C", "I", "F", "R") or a list of centering vectors. If given, `ops`
                are the primitive operations WITHOUT the centering translations, and positions that differ by
                a centering vector are the same site, so `positions` need only list the primitive-cell sites.
    """
    pairs, has_t = _normalize_ops(ops)
    shifts, r_from_shifts = _use_centering(pairs, centering)
    if lattice_mod is None:
        lattice_mod = has_t or shifts is not None
    P = _vectors(positions)
    keys = [_key(p, lattice_mod, shifts) for p in P]
    where = {}
    for j, k in enumerate(keys):
        where.setdefault(k, []).append(j)

    V = None
    if vectors is not None:
        V = _vectors(vectors)
        if len(V) != len(P):
            raise ValueError("the number of positions and vectors differ")

    orbit_ops = [[] for _ in P]
    stabilizer, unmatched = [], []
    for i, (R, t) in enumerate(pairs, start=1):
        k = _key(R * P[0] + t, lattice_mod, shifts)
        if k == keys[0]:
            stabilizer.append(i)  # site symmetry: depends on the position only
        js = where.get(k)
        if js is not None and V is not None:
            # With vectors, an operation maps the representative bond onto bond j only if the
            # position matches AND the vector matches (R V0 == +-V_j); several bonds may share a position.
            js = [j for j in js if _sign(R * V[0], V[j]) != 0]
        if not js:
            unmatched.append(i)
            continue
        for j in js:
            orbit_ops[j].append(i)

    vector_sign = None
    bond_stabilizer = bond_flips = bond_Rs = None
    if V is not None:
        vector_sign = [[_sign(pairs[i - 1][0] * V[0], V[j]) for i in orbit_ops[j]] for j in range(len(P))]
        # Bond symmetry: the operations (rotation part only; the translation part is never applied
        # to a vector) that fix the representative position AND map the representative vector onto
        # itself, ignoring its direction (R V == +V or R V == -V for all values of X, Y, Z).
        sgn = {i: _sign(pairs[i - 1][0] * V[0], V[0]) for i in stabilizer}
        bond_stabilizer = [i for i in stabilizer if sgn[i] != 0]
        bond_flips = [i for i in bond_stabilizer if sgn[i] == -1]  # operations that reverse the vector
        bond_Rs = [pairs[i - 1][0] for i in bond_stabilizer]

    Rs = [pairs[i - 1][0] for i in stabilizer]
    H_Rs = list({tuple(R): R for R in Rs}.values())
    all_Rs = list({tuple(R): R for R, _ in pairs}.values())
    r_centered = any(
        R == sp.eye(3) and tuple(sorted(c % 1 for c in t)) == (sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(2, 3))
        for R, t in pairs
    )

    r_centered = r_centered or r_from_shifts

    def oriented(Hs):
        try:
            return _oriented_symbol(all_Rs, list({tuple(R): R for R in Hs}.values()), r_centered)
        except ValueError:
            if not oriented_errors_ok:
                raise
            return None

    symbol = oriented(Rs)
    bond_kw = {}
    if bond_stabilizer is not None:
        bond_kw = dict(
            bond_stabilizer=bond_stabilizer,
            bond_order=len(bond_stabilizer),
            bond_flips=bond_flips,
            bond_op_types=[_op_type(R) for R in bond_Rs],
            bond_point_group=_point_group(bond_Rs),
            bond_symbol=oriented(bond_Rs),
            bond_multiplicity_ok=len(P) * len(bond_stabilizer) == len(pairs),
        )
    return SiteSymmetry(
        **bond_kw,
        symbol=symbol,
        label=label,
        representative=P[0],
        stabilizer=stabilizer,
        order=len(stabilizer),
        point_group=_point_group(Rs),
        op_types=[_op_type(R) for R in Rs],
        orbit_ops=orbit_ops,
        vector_sign=vector_sign,
        unmatched=unmatched,
        multiplicity_ok=len(set(keys)) * len(stabilizer) == len(pairs),  # distinct positions x site order
    )


# Centering translations of the conventional cell (ITA). R is given in hexagonal axes (obverse setting).
_H = sp.Rational(1, 2)
CENTERING = {
    "P": [],
    "A": [(0, _H, _H)],
    "B": [(_H, 0, _H)],
    "C": [(_H, _H, 0)],
    "I": [(_H, _H, _H)],
    "F": [(0, _H, _H), (_H, 0, _H), (_H, _H, 0)],
    "R": [(sp.Rational(2, 3), sp.Rational(1, 3), sp.Rational(1, 3)), (sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(2, 3))],
}


def expand_centering(ops, centering, order="centering-major"):
    """Expand the operations of a primitive cell with the centering translations of a non-primitive lattice.

    ops       : the operations WITHOUT centering (4x4, or 3x3 for t = 0), in the order of the data.
    centering : a lattice letter "P", "A", "B", "C", "I", "F", "R" (R: hexagonal axes), or a list of
                translation vectors, e.g. [(1/2, 1/2, 0)]. The identity translation (0, 0, 0) is always
                included automatically, so list only the non-trivial centering vectors ("P" = none).
    order     : "centering-major" (default, ITA numbering): all operations with the first translation
                come first, then all operations shifted by the 1st centering vector, then the 2nd, ...
                e.g. C2/m: 1-4 are the point operations, 5-8 are the same four plus (1/2,1/2,0).
                "operation-major": each operation is followed by its centered copies.
    Returns 4x4 sympy matrices whose translation parts are reduced to [0, 1).
    """
    pairs, _ = _normalize_ops(ops)
    if any(R == sp.eye(3) and any(c % 1 for c in t) for R, t in pairs):
        raise ValueError("the operations already contain centering translations (a pure translation was found)")
    vecs = CENTERING[centering] if isinstance(centering, str) else centering
    shifts = [sp.zeros(3, 1)]  # the identity translation is always the first one
    for v in vecs:  # a zero vector (mod 1) or a repeated vector is ignored
        c = sp.Matrix([sp.Rational(e) % 1 for e in v])
        if all(c != s for s in shifts):
            shifts.append(c)
    if order == "centering-major":
        items = [(R, t + c) for c in shifts for R, t in pairs]
    elif order == "operation-major":
        items = [(R, t + c) for R, t in pairs for c in shifts]
    else:
        raise ValueError(f"unknown order: {order!r}")
    out = []
    for R, t in items:
        M = sp.eye(4)
        M[:3, :3] = R
        M[:3, 3] = t.applyfunc(lambda e: e % 1)
        out.append(M)
    return out


def analyze_entry(ops, entry, *, label=None, **kw):
    """Analyze one (vector string, position string[, index]) tuple.
    The third element (the indices) is optional: it is only used to check the computed indices against it,
    the result is stored in index_matches_given."""
    vec, pos, *rest = entry
    given = rest[0] if rest else None
    res = site_symmetry(ops, pos, vec, label=label, **kw)
    if given is not None:
        _check_given(res, given, ops, label)
    return res


def _check_given(res, given, ops, label):
    top = max((abs(v) for g in given for v in g), default=0)
    if top > len(ops):
        warnings.warn(
            f"{label or ''}: the given indices refer to operation numbers up to {top}, but only "
            f"{len(ops)} operations were passed (missing centering translations? or pass "
            f"the primitive operations together with `centering=`)",
            stacklevel=3,
        )
    if res.vector_sign is not None:
        res.index_matches_given = res.signed_orbit_ops() == [sorted(list(g), key=abs) for g in given]


@dataclass
class BondSymmetry:
    """Symmetry of a bond: the vector v and its site c taken together as ONE object.

    The bond is the first entry of the list. `stabilizer` are the operations that fix the site c and map v
    onto +v or -v (only the rotation part acts on v), and `symbol` is the oriented symbol of that group.
    There is no separate site symbol here; the site-only symmetry is in site_symmetry().
    """

    label: str | None
    stabilizer: list  # numbers (1-based) of the operations that keep the bond (v, c) invariant
    flips: list  # those among them that reverse v (R v == -v)
    order: int
    point_group: str | None  # international symbol of the group (orientation not distinguished)
    symbol: str | None  # oriented symbol (ITA style: ..m, m2m, 4/mm.m, ...)
    op_types: list
    multiplicity: int  # number of bonds in the list
    signed_index: list  # per bond: operations that map the representative bond onto it (-: v reversed)
    closed: bool  # every operation maps the bond onto a listed bond and multiplicity * order == |G|
    index_matches_given: bool | None = None

    def __str__(self):
        head = f"{self.label}: " if self.label else ""
        s = (
            f"{head}bond symmetry {self.symbol or self.point_group} (order {self.order}), ops {self.stabilizer}"
            f", flipping {self.flips}, multiplicity {self.multiplicity}, signed ops {self.signed_index}"
        )
        if self.index_matches_given is not None:
            s += f", matches given: {self.index_matches_given}"
        if not self.closed:
            s += ", NOT CLOSED (the bond list is incomplete or inconsistent)"
        return s


def analyze_bond(ops, entry, *, label=None, **kw):
    """Symmetry of the bond (v, c) as one object, for one (vector string, position string[, index]) entry.

    Returns a BondSymmetry: the operations that keep the representative bond (the first entry) invariant,
    with ONE oriented symbol for the group. Keyword arguments are those of site_symmetry (e.g. centering=).
    The third element of `entry` (the indices) may be left out."""
    vec, pos, *rest = entry
    given = rest[0] if rest else None
    r = site_symmetry(ops, pos, vec, label=label, **kw)
    res = BondSymmetry(
        label=label,
        stabilizer=r.bond_stabilizer,
        flips=r.bond_flips,
        order=r.bond_order,
        point_group=r.bond_point_group,
        symbol=r.bond_symbol,
        op_types=r.bond_op_types,
        multiplicity=len(r.orbit_ops),
        signed_index=r.signed_orbit_ops(),
        closed=bool(r.bond_multiplicity_ok and not r.unmatched),
    )
    if given is not None:
        res.index_matches_given = res.signed_index == [sorted(list(g), key=abs) for g in given]
        top = max((abs(v) for g in given for v in g), default=0)
        if top > len(ops):
            warnings.warn(
                f"{label or ''}: the given indices refer to operation numbers up to {top}, but only "
                f"{len(ops)} operations were passed (missing centering translations? or pass "
                f"the primitive operations together with `centering=`)",
                stacklevel=2,
            )
    return res


def analyze_table(ops_by_group, table, *, bond=False, **kw):
    """Analyze a whole {group name: {label: entry}} table.
    ops_by_group is {group name: array of symmetry operations}.
    Returns {group name: {label: SiteSymmetry}}, or BondSymmetry objects if bond=True.
    `centering` may be a single letter or a dict {group name: letter}."""
    func = analyze_bond if bond else analyze_entry
    cen = kw.pop("centering", None)
    return {
        grp: {
            lab: func(ops_by_group[grp], ent, label=lab, centering=cen.get(grp) if isinstance(cen, dict) else cen, **kw)
            for lab, ent in wy.items()
        }
        for grp, wy in table.items()
        if grp in ops_by_group
    }


@dataclass
class SetSymmetry:
    """Symmetry of a set of (v, c) pairs: the operations that map the whole set onto itself."""

    label: str | None
    n_elements: int
    stabilizer: list  # numbers (1-based) of the operations that map the set onto itself
    order: int
    point_group: str | None  # international symbol of the group (orientation not distinguished)
    symbol: str | None  # oriented symbol (same conventions as SiteSymmetry.symbol)
    op_types: list
    permutations: dict  # op number -> [(index of the image element, +1/-1), ...] for every element j
    orbits: list  # orbits of the elements under the symmetry group (lists of element indices)
    closed_under_group: bool  # True if every operation of the group maps the set onto itself
    escapes: dict  # op number -> elements whose image is NOT in the set (only for ops outside the group)

    def __str__(self):
        head = f"{self.label}: " if self.label else ""
        s = (
            f"{head}set symmetry {self.symbol or self.point_group} (order {self.order}) of {self.n_elements} "
            f"elements, ops {self.stabilizer}, orbits {self.orbits}"
        )
        if not self.closed_under_group:
            s += f", NOT invariant under {len(self.escapes)} operation(s) of the group"
        return s


def set_symmetry(
    ops, vectors, centers, *, ignore_direction=True, lattice_mod=None, label=None, oriented_errors_ok=False, centering=None
):
    """Symmetry of the set {(v_j, c_j)} of vectors v_j attached to sites c_j.

    An operation (R, t) is a symmetry of the set if it maps every element (v, c) onto an element of
    the set: the site goes to R c + t (compared modulo lattice translations) and the vector goes to R v
    (only the rotation part acts on a vector; partial translations are excluded). Components are compared
    identically in the parameters x, y, z / X, Y, Z. The set need not be an orbit: for a single element this
    is the bond symmetry, for a complete orbit it is the whole group, and for a cluster of bonds it is the
    subgroup that keeps the cluster invariant.

    ops              : symmetry operations (3x3 or 4x4); include the centering translations (expand_centering).
    vectors, centers : same length; element j is (vectors[j], centers[j]).
    ignore_direction : if True (default) v and -v are the same bond; the sign of each mapping is recorded
                       (permutations[op][j] = (image index, -1) if v is reversed). If False, v must be kept.
    """
    pairs, has_t = _normalize_ops(ops)
    shifts, r_from_shifts = _use_centering(pairs, centering)
    if lattice_mod is None:
        lattice_mod = has_t or shifts is not None
    V, C = _vectors(vectors), _vectors(centers)
    if len(V) != len(C):
        raise ValueError("the numbers of vectors and centers differ")

    def vkey(v):
        a = tuple(sp.expand(e) for e in v)
        return min(a, tuple(sp.expand(-e) for e in v), key=str) if ignore_direction else a

    index = {}
    for j, (v, c) in enumerate(zip(V, C)):
        k = (_key(c, lattice_mod, shifts), vkey(v))
        if k in index:
            raise ValueError(
                f"elements {index[k]} and {j} are identical (modulo lattice translations"
                f"{' and the sign of v' if ignore_direction else ''})"
            )
        index[k] = j

    perms, escapes, stab = {}, {}, []
    for i, (R, t) in enumerate(pairs, start=1):
        img, esc = [], []
        for j, (v, c) in enumerate(zip(V, C)):
            rv = R * v
            m = index.get((_key(R * c + t, lattice_mod, shifts), vkey(rv)))
            if m is None:
                esc.append(j)
            else:
                img.append((m, _sign(rv, V[m]) if ignore_direction else 1))
        if esc:
            escapes[i] = esc
        else:
            perms[i] = img
            stab.append(i)

    parent = list(range(len(V)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for img in perms.values():
        for j, (m, _) in enumerate(img):
            parent[find(j)] = find(m)
    groups = {}
    for j in range(len(V)):
        groups.setdefault(find(j), []).append(j)

    Rs = [pairs[i - 1][0] for i in stab]
    all_Rs = list({tuple(R): R for R, _ in pairs}.values())
    r_centered = any(
        R == sp.eye(3) and tuple(sorted(c % 1 for c in t)) == (sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(2, 3))
        for R, t in pairs
    )
    r_centered = r_centered or r_from_shifts
    try:
        symbol = _oriented_symbol(all_Rs, list({tuple(R): R for R in Rs}.values()), r_centered)
    except ValueError:
        if not oriented_errors_ok:
            raise
        symbol = None
    return SetSymmetry(
        label=label,
        n_elements=len(V),
        stabilizer=stab,
        order=len(stab),
        point_group=_point_group(Rs),
        symbol=symbol,
        op_types=[_op_type(R) for R in Rs],
        permutations=perms,
        orbits=list(groups.values()),
        closed_under_group=len(stab) == len(pairs),
        escapes=escapes,
    )


def analyze_set(ops, entry, *, label=None, **kw):
    """set_symmetry for an entry (vector string, position string, ...) of the data: element j is
    (vectors[j], positions[j]). The third element (the indices) is not used."""
    return set_symmetry(ops, entry[0], entry[1], label=label, **kw)


def diagnose(ops, entry, *, only_differing=False, **kw):
    """Explain, operation by operation, how the computed bond mapping compares with the given indices.

    For every operation i this shows the image of the representative position, the image R*V0 of the
    representative vector, and the entries (with +/- for kept/reversed vector) that i maps the
    representative bond onto, as computed here and as given in the data. Operations where the two
    differ are flagged, which usually shows at once whether the cause is the numbering of the
    operations (e.g. an offset or a different order), the convention of the matrices (transposed /
    inverse / Cartesian basis), a missing centering translation, or the vector basis.
    """
    vec, pos, given = entry
    pairs, _ = _normalize_ops(ops)
    P, V = _vectors(pos), _vectors(vec)
    res = site_symmetry(ops, pos, vec, **kw)
    got, want = {}, {}
    for j, (ops_j, sg) in enumerate(zip(res.orbit_ops, res.vector_sign)):
        for i, s in zip(ops_j, sg):
            got.setdefault(i, []).append((j, s))
    for j, g in enumerate(given):
        for v in g:
            want.setdefault(abs(v), []).append((j, 1 if v > 0 else -1))

    def fmt(lst):
        return ",".join(f"{j}{'+' if s > 0 else '-'}" for j, s in sorted(lst)) or "none"

    lines, n_diff = [], 0
    for i, (R, t) in enumerate(pairs, start=1):
        a, b = sorted(got.get(i, [])), sorted(want.get(i, []))
        n_diff += a != b
        if only_differing and a == b:
            continue
        lines.append(
            f"{i:3d}  image {list(R * P[0] + t)}  R*V0 {list(R * V[0])}  "
            f"computed entry {fmt(a):8s} given entry {fmt(b):8s}{'  <<< differs' if a != b else ''}"
        )
    extra = sorted(set(want) - set(range(1, len(pairs) + 1)))
    head = (
        f"{n_diff} of {len(pairs)} operations differ (entry numbers are 0-based positions in the list; "
        f"+: V kept, -: V reversed)"
    )
    if extra:
        head += (
            f"\nMISMATCH: the data refers to operation numbers {extra}, but only {len(pairs)} operations "
            f"were passed. The list passed is shorter than the one used to make the data "
            f"(missing centering translations? e.g. a C-centered group needs 2x the point operations)."
        )
    return head + "\n" + "\n".join(lines)


# ==================================================
def create_bond_symmetry():
    from data_wyckoff_bond import wyckoff_bond_data  # same directory (run as a script)
    from multipie import Group

    so = {}
    center = {}
    for key in wyckoff_bond_data.keys():
        g = Group(key)
        so[key] = g.symmetry_operation["fractional"]
        if g.info.lattice == "0":
            center[key] = None
        else:
            center[key] = g.info.lattice

    result = analyze_table(so, wyckoff_bond_data, bond=True, centering=center)

    bond = {}
    for k, v in result.items():
        bond[k] = {}
        for wp, ss in v.items():
            if not ss.index_matches_given:
                print("not match", k, wp)
            v, c, op = wyckoff_bond_data[k][wp]
            bond[k][wp] = (v, c, op, ss.symbol)

    return bond


# ==================================================
def create_file():
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "data_wyckoff_bond_sym.py")

    wyckoff_bond = create_bond_symmetry()

    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print(f"wyckoff_bond_data = {wyckoff_bond}", file=f)


# ==================================================
if __name__ == "__main__":
    create_file()
