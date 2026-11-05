# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPieGen developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
For multipolar spherical harmonics (sympy).

The internal basis is no longer fixed to tesseral.
It is defined by an InternalBasis object, which is passed explicitly by the `basis` argument
(no module-level mutable state; the default is "tesseral").
"""

import functools

import numpy as np
import sympy as sp
from sympy.physics.wigner import wigner_3j

TOL = 1e-11

# ==================================================
# max. rank of spherical harmonics.
max_rank = 11

# ==================================================
# max. internal rank of spherical harmonics.
max_s_rank = 3

# ==================================================
# rank name.
rank_name = {
    0: "monopole",
    1: "dipole",
    2: "quadrupole",
    3: "octupole",
    4: "hexadecapole",
    5: "dotriacontapole",
    6: "tetrahexacontapole",
    "monopole": 0,
    "dipole": 1,
    "quadrupole": 2,
    "octupole": 3,
    "hexadecapole": 4,
    "dotriacontapole": 5,
    "tetrahexacontapole": 6,
}

# ==================================================
_x, _y, _z = sp.symbols("x y z", real=True)
_r = sp.sqrt(_x**2 + _y**2 + _z**2)
_rp = sp.sqrt(_x**2 + _y**2)


# ==================================================
def _create_unitary_matrix_tesseral(l):
    """
    Create unitary matrix to tesseral basis.

    Args:
        l (int): rank.

    Returns:
        - (ndarray) -- unitary matrix, [[sympy]] (2l+1,2l+1), <m|gamma> (m=l, l-1, ..., -l).
    """
    U = np.full((2 * l + 1, 2 * l + 1), sp.S(0))
    for m in range(l, 0, -1):
        U[2 * (l - m), l - m] = (-1) ** m / sp.sqrt(2)
        U[2 * (l - m), l + m] = 1 / sp.sqrt(2)
        U[2 * (l - m) + 1, l - m] = -sp.I * ((-1) ** m) / sp.sqrt(2)
        U[2 * (l - m) + 1, l + m] = sp.I / sp.sqrt(2)
    U[2 * l, l] = sp.S(1)

    return U.T


# ==================================================
def _create_unitary_matrix_spherical(l):
    """
    Create unitary matrix to spherical (complex) basis, i.e., identity.

    Args:
        l (int): rank.

    Returns:
        - (ndarray) -- unitary matrix, [[sympy]] (2l+1,2l+1), <m|gamma> (m=l, l-1, ..., -l).
    """
    U = np.full((2 * l + 1, 2 * l + 1), sp.S(0))
    for i in range(2 * l + 1):
        U[i, i] = sp.S(1)

    return U


# ==================================================
# built-in definitions, f(s) -> U[m, alpha] (m=s, s-1, ..., -s).
_BUILTIN_DEFINITIONS = {
    "tesseral": _create_unitary_matrix_tesseral,
    "spherical": _create_unitary_matrix_spherical,
}


# ==================================================
def _as_sympy_array(U, s):
    """
    Convert to (2s+1,2s+1) object ndarray of sympy.
    """
    U = np.asarray(U)
    n = 2 * s + 1
    if U.shape != (n, n):
        raise ValueError(f"internal basis for s={s} must have shape {(n, n)}, got {U.shape}.")
    ret = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            ret[i, j] = sp.sympify(U[i, j])
    return ret


# ==================================================
def _check_unitary(U, s):
    """
    Check unitarity, U^dagger U = 1.
    """
    P = U.T.conjugate() @ U
    for i in range(2 * s + 1):
        for j in range(2 * s + 1):
            d = sp.simplify(P[i, j] - (1 if i == j else 0))
            if d == 0:
                continue
            try:
                ok = abs(complex(d)) < 1e-9
            except TypeError:
                ok = False
            if not ok:
                raise ValueError(f"internal basis for s={s} is not unitary.")


# ==================================================
class InternalBasis:
    """
    Definition of internal basis, with its own caches (no shared state between instances).

    Args:
        definition (str, callable, dict or InternalBasis, optional): definition of the basis.
            - "tesseral" / "spherical": built-in.
            - callable: f(s) -> array-like (2s+1,2s+1), U[m,alpha] (m=s,s-1,...,-s), unitary,
              same convention as _create_unitary_matrix_tesseral.
            - dict: {s: U} for s = 0..max_s_rank.
        name (str, optional): label.

    Note:
        - for s>0, alpha-th component of Olm^{(s)}(k) is sum_m Olm_m * conj(U[m,alpha]).
        - reuse the same instance to benefit from the caches.
    """

    def __init__(self, definition="tesseral", name=None):
        if isinstance(definition, InternalBasis):
            func, label = definition._func, definition.name
        elif isinstance(definition, str):
            if definition not in _BUILTIN_DEFINITIONS:
                raise ValueError(f"unknown internal basis: {definition}")
            func, label = _BUILTIN_DEFINITIONS[definition], definition
        elif isinstance(definition, dict):
            func, label = dict(definition).__getitem__, "custom(dict)"
        elif callable(definition):
            func, label = definition, "custom"
        else:
            raise TypeError("definition must be str, callable, dict or InternalBasis.")

        self._func = func
        self.name = name if name is not None else label
        self._U = {}
        self._Uc = {}
        self._sym = {}  # symbolic harmonics cache.
        self._num = {}  # lambdified harmonics cache.

        for s in range(max_s_rank + 1):  # validate all ranks at construction.
            self.U(s)

    def __repr__(self):
        return f"InternalBasis({self.name!r})"

    def U(self, s):
        """
        Get unitary matrix <m|alpha> (object ndarray of sympy) for internal rank s.
        """
        try:
            return self._U[s]
        except KeyError:
            pass
        if not (0 <= s <= max_s_rank):
            raise ValueError(f"internal rank s must be in 0..{max_s_rank}.")
        U = _as_sympy_array(self._func(s), s)
        _check_unitary(U, s)
        self._U[s] = U
        self._Uc[s] = U.conjugate()
        return U

    def Uc(self, s):
        """
        Get complex conjugate of U.
        """
        self.U(s)
        return self._Uc[s]


# ==================================================
_builtin_bases = {}


def get_internal_basis(basis="tesseral"):
    """
    Get InternalBasis object.

    Args:
        basis (str, callable, dict or InternalBasis, optional): definition of internal basis.

    Returns:
        - (InternalBasis) -- internal basis.

    Note:
        - built-in names are created once and reused.
        - for callable/dict, pass an InternalBasis instance to Olm etc. to keep caches.
    """
    if isinstance(basis, InternalBasis):
        return basis
    if isinstance(basis, str):
        if basis not in _builtin_bases:
            _builtin_bases[basis] = InternalBasis(basis)
        return _builtin_bases[basis]
    return InternalBasis(basis)


# ==================================================
def is_transform_real(U):
    """
    Check if transformed harmonics is real.

    Args:
        U (array-like): coefficient vector, <m|gamma>, for gamma-harmonics (m=l,l-1,...,-l).

    Returns:
        - (bool) -- is real ?
    """
    U = np.asarray(U)
    l = (U.shape[0] - 1) // 2
    if U.dtype == object:
        sgn = np.asarray([sp.S(1) if m % 2 == 0 else sp.S(-1) for m in range(l, -l - 1, -1)])
    else:
        sgn = np.asarray([1.0 if m % 2 == 0 else -1.0 for m in range(l, -l - 1, -1)])
    Ub = sgn * U[::-1].conjugate()

    if U.dtype == object:
        return np.array_equal(U, Ub)
    else:
        return np.allclose(U, Ub, TOL, 0.01 * TOL)


# ==================================================
@functools.lru_cache(maxsize=None)
def _w3j(j1, j2, j3, m1, m2, m3):
    """
    Wigner 3j symbol (cached).
    """
    return wigner_3j(j1, j2, j3, m1, m2, m3)


# ==================================================
@functools.lru_cache(maxsize=None)
def _olm_base(l, m):
    """
    Spherical harmonics (scaled), sqrt(4pi/(2l+1)) r^l Ylm, symbolic in x,y,z (cached).
    """
    if m == 0:
        olm = _r**l * sp.legendre(l, _z / _r)
    else:
        s = m > 0
        m = abs(m)
        olm = sp.sqrt(sp.factorial(l - m) / sp.factorial(l + m)) * _r**l * sp.assoc_legendre(l, m, _z / _r) / _rp ** (m)
        if s:
            olm *= (_x + sp.I * _y) ** m
        else:
            olm *= (-1) ** m * (_x - sp.I * _y) ** m

    return sp.radsimp(sp.factor(olm))


# ==================================================
def _olm_symbolic(l, m, s, k, basis, factor):
    """
    Multipolar spherical harmonics, symbolic in x,y,z (cached in basis).

    Returns:
        - (tuple) -- sympy expressions (length 2s+1; 1 for s=0).
    """
    key = (l, m, s, k, factor)
    res = basis._sym.get(key)
    if res is not None:
        return res

    if factor:
        fin = lambda e: sp.factor(sp.radsimp(e))
    else:
        fin = lambda e: sp.expand(sp.radsimp(e))

    if s == 0:
        res = (fin(_olm_base(l, m)),)
    else:
        n_s = 2 * s + 1
        ev = basis.Uc(s)
        acc = [sp.S(0)] * n_s
        n1 = max(m - l - k, -s)
        n2 = min(m + l + k, s)
        for n in range(n1, n2 + 1):
            w3j = _w3j(l + k, l, s, m - n, -m, n)
            if w3j == 0:
                continue
            wo = w3j * _olm_base(l + k, m - n)
            row = ev[s - n]
            for a in range(n_s):
                if row[a] != 0:
                    acc[a] += wo * row[a]
        pre = sp.I ** (k + s) * (-1) ** (l + m) * sp.sqrt(2 * l + 1)
        res = tuple(fin(pre * e) for e in acc)

    basis._sym[key] = res
    return res


# ==================================================
def _olm_numeric(l, m, s, k, basis, factor, r):
    """
    Numerical evaluation by lambdified function (cached in basis).

    Returns:
        - (list or None) -- complex values, None if not finite (fallback to substitution).
    """
    key = (l, m, s, k, factor)
    f = basis._num.get(key)
    if f is None:
        f = sp.lambdify((_x, _y, _z), list(_olm_symbolic(l, m, s, k, basis, factor)), "numpy")
        basis._num[key] = f
    try:
        with np.errstate(all="ignore"):
            vals = [complex(v) for v in f(*r)]
    except (TypeError, ValueError, ZeroDivisionError):
        return None
    if not all(np.isfinite(v.real) and np.isfinite(v.imag) for v in vals):
        return None
    return vals


# ==================================================
def Olm(l, m, s=0, k=0, rv=None, factor=False, basis="tesseral"):
    """
    Multipolar spherical harmonics (scaled), Olm^{(s)}(k).

    Args:
        l (int): rank.
        m (int): component (-l<= m <= l).
        s (int, optional): internal rank.
        k (int, optional): internal component (-s<= k <=s).
        rv (array-like, optional): (x,y,z) symbol/str/sympy/float values.
        factor (bool, optional): factorize resultant expression ? (sometime slow)
        basis (str, callable, dict or InternalBasis, optional): definition of internal basis.

    Returns:
        - (sympy or complex) -- multipolar spherical harmonics.

    Note:
        - internal basis (s>0) in order of alpha of the given basis,
          for tesseral, [(s,"c"),(s,"s"),(s-1,"c"),...,(0,"c")].
        - for s > 0, return ndarray of [sympy(2s+1)].
    """
    assert 0 <= s <= max_s_rank
    basis = get_internal_basis(basis)

    if rv is None:
        vals = _olm_symbolic(l, m, s, k, basis, factor)
        numeric = False
    else:
        rv = np.asarray(rv)
        vals = None
        if rv.dtype.kind in "iufc":
            vals = _olm_numeric(l, m, s, k, basis, factor, rv[:3].tolist())
        numeric = vals is not None
        if not numeric:
            sub = {_x: rv[0], _y: rv[1], _z: rv[2]}
            vals = [sp.radsimp(e.subs(sub)) for e in _olm_symbolic(l, m, s, k, basis, factor)]
            if any(bool(v.free_symbols) for v in vals) or rv.dtype == object:
                if factor:
                    vals = [sp.factor(sp.radsimp(v)) for v in vals]
                else:
                    vals = [sp.expand(sp.radsimp(v)) for v in vals]
            else:
                vals = [complex(v) for v in vals]
                numeric = True

    if numeric:
        molm = np.asarray(vals, dtype=complex)
    else:
        molm = np.empty(len(vals), dtype=object)
        molm[:] = vals

    if molm.size == 1:
        molm = molm.item()

    return molm


# ==================================================
def create_harmonics_set(l, s=0, k=0, rv=None, sv=None, factor=False, basis="tesseral"):
    """
    Create harmonics set (m=l, l-1, ..., -l).

    Args:
        l (int): rank.
        s (int, optional): internal rank.
        k (int, optional): internal component (-s<= k <=s).
        rv (array-like, optional): (x,y,z) symbol/str/sympy/float values.
        sv (array-like, optional): internal variable.
        factor (bool, optional): factorize resultant expression ? (sometime slow)
        basis (str, callable, dict or InternalBasis, optional): definition of internal basis.

    Returns:
        - (ndarray) -- multipolar spherical harmonics set.
    """
    basis = get_internal_basis(basis)
    olm = np.asarray([Olm(l, m, s, k, rv, factor, basis) for m in range(l, -l - 1, -1)])

    if sv is not None and s > 0:
        b = create_internal_basis(s, sv, factor=factor, basis=basis)
        olm = olm @ b

    return olm


# ==================================================
_factor_expand = np.frompyfunc(lambda i: sp.factor(sp.expand(i)), 1, 1)


def transform_harmonics(Olm_set, U):
    """
    Transform harmonics by O_{lg} = sum_m U_{m,g}O_{lm}.

    Args:
        Olm_set (array-like): Olm [m, point, internal comp.].
        U (array-like): U [m,g].

    Returns:
        - (ndarray) -- transformed harmonics [g, point, internal comp.].
    """
    Olm_set = np.asarray(Olm_set)
    U = np.asarray(U)

    Olg = _factor_expand(U.T @ Olm_set)

    if Olg.size == 1:
        Olg = Olg.item()

    return Olg


# ==================================================
def create_internal_basis(s, sv=None, factor=False, basis="tesseral"):
    """
    Create internal basis (harmonics of rank s transformed by the unitary matrix of the basis).

    Args:
        s (int, optional): internal rank.
        sv (array-like, optional): (x,y,z) symbol/str/sympy/float values.
        factor (bool, optional): factorize resultant expression ? (sometime slow)
        basis (str, callable, dict or InternalBasis, optional): definition of internal basis.

    Returns:
        - (ndarray) -- internal basis.
    """
    basis = get_internal_basis(basis)
    return transform_harmonics(create_harmonics_set(s, rv=sv, factor=factor), basis.U(s))
