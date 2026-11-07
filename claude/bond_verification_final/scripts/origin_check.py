"""Which ITA origin choice do the file's operations use? (groups with two origin choices)
Compare the file's operation set with spglib's database for each Hall number of the group."""
import sys; sys.path.insert(0, ".")
import numpy as np, spglib, sympy as sp
import check_complete as C
types = {k: spglib.get_spacegroup_type(k) for k in range(1, 531)}
names = {C.sg_number(n): n for n in C.data if "^" in n}
fr = lambda M: np.array([[float(sp.Rational(e)) for e in r] for r in M.tolist()])
for num in range(1, 231):
    hs = [k for k, t in types.items() if t.number == num]
    if not any(types[h].choice == "2" for h in hs):
        continue
    ops = C.rebuild_ops(C.data[names[num]])[1]
    S_f = {(tuple(fr(M[:3, :3]).astype(int).ravel()), tuple(np.round(fr(M[:3, 3:4]).ravel() % 1, 4))) for M in ops}
    out = {}
    for h in hs:
        db = spglib.get_symmetry_from_database(h)
        S_db = {(tuple(r.ravel()), tuple(np.round(t % 1, 4))) for r, t in zip(db["rotations"], db["translations"])}
        out[types[h].choice] = S_f <= S_db
    print(num, names[num], "origin 1:", out.get("1"), "origin 2:", out.get("2"))
