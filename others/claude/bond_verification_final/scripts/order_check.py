"""Compare the file's operations (rebuilt from the general-position entry) with pyxtal (ITA order) and spglib."""
import sys; sys.path.insert(0, ".")
import pickle, numpy as np, spglib, sympy as sp
import check_complete as C
r = pickle.load(open("complete_all.pkl", "rb"))
sg = [n for n in r if "^" in n]
print("pyxtal same order:", sum(r[n]["same_as_ita"] is True for n in sg), "/", len(sg), "differ:", [n for n in sg if r[n]["same_as_ita"] is not True])
types = {k: spglib.get_spacegroup_type(k) for k in range(1, 531)}
fr = lambda M: np.array([[float(sp.Rational(e)) for e in row] for row in M.tolist()])
oo = os_ = 0
for n in sg:
    ops = C.rebuild_ops(C.data[n])[1]; m = len(ops)
    Rf = [fr(M[:3, :3]).astype(int) for M in ops]; Tf = [fr(M[:3, 3:4]).ravel() % 1 for M in ops]
    o = s = False
    for h in [k for k, t in types.items() if t.number == r[n]["number"]]:
        db = spglib.get_symmetry_from_database(h); R, T = db["rotations"], db["translations"]
        o |= len(R) >= m and all(np.array_equal(R[i], Rf[i]) and np.allclose((T[i] - Tf[i] + .5) % 1 - .5, 0, atol=1e-6) for i in range(m))
        Sdb = {(tuple(a.ravel()), tuple(np.round(t % 1, 4))) for a, t in zip(R, T)}
        s |= {(tuple(a.ravel()), tuple(np.round(t, 4))) for a, t in zip(Rf, Tf)} <= Sdb
    oo += o; os_ += s
print("spglib same order:", oo, "/", len(sg), " same set:", os_, "/", len(sg))
