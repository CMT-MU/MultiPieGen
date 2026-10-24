import numpy as np
from multipie.util.util import str_to_sympy
from multipie import Group

TOL = 1e-6


# ==================================================
def is_int(v):
    v = np.asarray(v, dtype=float)
    return bool(np.all(np.abs(v - np.round(v)) < TOL))


# ==================================================
def centering(g):
    """Centering translations of group g as 3-vectors; (0,0,0) is always included.

    NOTE: adapt this function if "plus_set" is stored in another form
          (here: a sequence of 3-vectors, or of 4x4 matrices, convertible to float).
    """
    cen = [np.zeros(3)]
    for c in g.symmetry_operation["plus_set"]:
        c = np.asarray(c, dtype=float).ravel()
        cen.append(c.reshape(4, 4)[:3, 3] if c.size == 16 else c[:3])
    return cen


# ==================================================
def same_mod_lattice(a, b, cen):
    """4x4 operations a and b are equal up to a lattice translation (Z^3 + centering vectors)."""
    if np.abs(a[:3, :3] - b[:3, :3]).max() > TOL:
        return False
    dt = a[:3, 3] - b[:3, 3]
    return any(is_int(dt - c) for c in cen)


# ==================================================
def n_classes(so):
    """number of distinct linear parts = order of the point group."""
    return len({tuple(np.round(o[:3, :3]).astype(int).ravel()) for o in so})


# ==================================================
def check(c_group, p_group):
    cg = Group(c_group)
    pg = Group(p_group)
    cosets, R = data[c_group][p_group]  # cosets[0] = parent operation numbers (1-based) of child operations
    R = np.array(str_to_sympy(R).astype(float))  # 4x4, x_parent = M x_child + d
    Ri = np.linalg.inv(R)

    p_so = pg.symmetry_operation["fractional"].astype(float)
    c_so = cg.symmetry_operation["fractional"].astype(float)
    p_cen = centering(pg)
    c_cen = centering(cg)

    err = []
    op = list(cosets[0])
    n_c, n_p = n_classes(c_so), n_classes(p_so)

    def find(a):  # 1-based number of the parent operation equal to a modulo lattice
        for k, b in enumerate(p_so):
            if same_mod_lattice(a, b, p_cen):
                return k + 1
        return None

    # --- 1. number / injectivity -----------------------------------------------------------
    if len(op) != n_c:
        err.append(f"len(op)={len(op)} but the child point group has {n_c} operations.")
    if len(set(op)) != len(op):
        err.append("op is not injective.")

    # --- 2. each child operation, transformed to the parent's coordinates, is the mapped parent operation
    for no, m in enumerate(op):
        t = R @ c_so[no] @ Ri  # child operation in parent coordinates (not reduced)
        q = p_so[m - 1]
        if not same_mod_lattice(t, q, p_cen):
            what = "linear part" if np.abs(t[:3, :3] - q[:3, :3]).max() > TOL else "translation (not a lattice vector)"
            err.append(f"child op {no + 1} -> parent op {m}: {what} differs.")

    # --- 3. the image is closed under products (a subgroup of the parent point group)
    image = set(op)
    for a in op:
        for b in op:
            k = find(p_so[a - 1] @ p_so[b - 1])
            if k not in image:
                err.append(f"image not closed: parent ops {a}, {b} give {k}.")

    # --- 4. M maps the child lattice (Z^3 + centering) into the parent lattice
    for v in list(np.eye(3)) + c_cen[1:]:
        w = R[:3, :3] @ v
        if not any(is_int(w - c) for c in p_cen):
            err.append(f"child lattice vector {np.round(v, 4).tolist()} is not mapped into the parent lattice.")

    # --- 5. coset decomposition: every parent operation exactly once, row k = g_k o (row 0)
    flat = [x for row in cosets for x in row]
    if sorted(flat) != list(range(1, n_p + 1)):
        err.append("cosets do not contain every parent operation exactly once.")
    for row in cosets[1:]:
        g = row[0]
        for a, x in zip(op, row):
            k = find(p_so[g - 1] @ p_so[a - 1])
            if k != x:
                err.append(f"coset row {row[:4]}...: g({g}) o a({a}) = {k}, not {x}.")
                break
    return err


# ==================================================
from data_parent_sg import data_parent_sg as data

err = 0
for c_group, v in data.items():
    for p_group in v.keys():
        print(f"{p_group} => {c_group}")
        lst = check(c_group, p_group)
        if lst:
            print(f"error in {p_group} => {c_group}.")
            for i in lst:
                print("   ", i)
            err += 1
if err == 0:
    print("no error.")
