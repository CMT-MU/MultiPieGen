"""pyxtal 同梱の ITA 表 (一般位置の操作と各 Wyckoff 位置の座標) と spglib の記号を突き合わせるデータ層。

records(num) -> list of dict(letter, mult, coords(list[str]), symbol(spglib), ok)
  coords[0] は ITA の第1座標 (向き付き記号の基準)。
  centering を含む場合 coords には conventional cell 全体の座標が並ぶ。
"""
import ast
import csv
import functools
import sys

import numpy as np
import sympy as sp
import spglib

sys.path.insert(0, "/home/claude")
csv.field_size_limit(10**9)
DB = "/tmp/claude-0/-home-claude/ebc6dc5f-4f46-5d72-ab5d-2c96f63a86f3/scratchpad/pyx/x/pyxtal/database/"
x, y, z = sp.symbols("x y z")
LETTERS = "abcdefghijklmnopqrstuvwxyz" + "A"


@functools.lru_cache(None)
def _wl():
    out = {}
    with open(DB + "wyckoff_list.csv") as f:
        for row in csv.reader(f):
            if row and row[0].isdigit() and len(row) > 1 and row[1].strip():
                out[int(row[0])] = ast.literal_eval(row[1])
    return out


def parse_triplet(s):
    import re
    out = []
    for t in s.split(","):
        t = re.sub(r"(\d)\s*([xyz])", r"\1*\2", t.strip())
        out.append(sp.sympify(t, locals={"x": x, "y": y, "z": z}))
    return out


def op_matrix(s):
    """'-x+y, x, z+1/2' -> 4x4 sympy 行列 (列ベクトル規約 r' = R r + t)"""
    exprs = parse_triplet(s)
    M = sp.eye(4)
    for i, e in enumerate(exprs):
        e = sp.expand(e)
        for j, v in enumerate((x, y, z)):
            M[i, j] = e.coeff(v)
        M[i, 3] = e.subs({x: 0, y: 0, z: 0})
    return M


def ops_for(num):
    """一般位置の座標リストから 4x4 操作 (ITA の並び順)。"""
    gen = _wl()[num][0]
    return [op_matrix(s) for s in gen]


def lattice_for(num):
    if num <= 2:
        return np.array([[3.1, 0.2, 0.1], [0.3, 4.3, 0.2], [0.1, 0.4, 5.7]])
    if num <= 15:
        b = np.deg2rad(100)
        return np.array([[4.0, 0, 0], [0, 5.0, 0], [6.0 * np.cos(b), 0, 6.0 * np.sin(b)]])
    if num <= 74:
        return np.diag([3.1, 4.3, 5.7])
    if num <= 142:
        return np.diag([3.1, 3.1, 5.3])
    if num >= 195:
        return np.diag([4.1, 4.1, 4.1])
    a, c = 3.3, 5.2
    return np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])


VALS = {x: 0.1234, y: 0.2711, z: 0.3619}


def numeric(coord_str):
    return np.array([float(e.subs(VALS)) for e in parse_triplet(coord_str)])


def mat_np(M):
    return np.array(M.tolist(), dtype=float)


def _halls(num):
    """この空間群の hall 番号の候補 (既定 + 原点選択 2 / 六方軸 H があればそれも)。"""
    hs = [h for h in range(1, 531) if spglib.get_spacegroup_type(h).number == num]
    return hs


def oracle_symbols(num, wps):
    """各 Wyckoff 位置を spglib で評価し (letter, symbol, ok, mult) を返す。
    spglib の setting (hall 番号) は、第1座標の Wyckoff 文字と多重度が合うものを探して決める。"""
    ops = [mat_np(M) for M in ops_for(num)]
    lat = lattice_for(num)
    n = len(wps)
    halls = _halls(num)
    marker_pt = np.array([0.0731, 0.1957, 0.0513])
    marker = np.array([(M[:3, :3] @ marker_pt + M[:3, 3]) % 1.0 for M in ops])
    res = []
    for i, coords in enumerate(wps):
        letter = LETTERS[n - 1 - i]
        p = numeric(coords[0])
        orb = []
        for M in ops:
            q = (M[:3, :3] @ p + M[:3, 3]) % 1.0
            if not any(np.allclose(q, o, atol=1e-6) or np.allclose(np.abs(q - o), 1, atol=1e-6) for o in orb):
                orb.append(q)
        orb = np.array(orb)
        cell = (lat, np.vstack([orb, marker]), [1] * len(orb) + [2] * len(marker))
        got = None
        for h in [None] + halls:
            ds = spglib.get_symmetry_dataset(cell, symprec=1e-4, hall_number=h or 0)
            if ds is None or ds.number != num:
                continue
            if (ds.wyckoffs[0] == letter
                    and len(set(ds.wyckoffs[:len(orb)])) == 1
                    and np.allclose(np.array(ds.transformation_matrix), np.eye(3), atol=1e-4)
                    and np.allclose(np.array(ds.origin_shift) - np.round(ds.origin_shift), 0, atol=1e-4)):
                got = (letter, ds.site_symmetry_symbols[0], True, len(orb))
                break
        res.append(got or (letter, None, False, len(orb)))
    return res


@functools.lru_cache(None)
def records(num):
    wps = _wl()[num]
    orc = oracle_symbols(num, wps)
    out = []
    for coords, (letter, sym, ok, mult) in zip(wps, orc):
        out.append(dict(letter=letter, mult=mult, coords=coords, symbol=sym, ok=ok))
    return out


if __name__ == "__main__":
    for num in map(int, sys.argv[1:]):
        recs = records(num)
        print(f"== {num}: {len(recs)} WPs, {sum(r['ok'] for r in recs)} ok")
        for r in reversed(recs):
            print(f"  {r['mult']:>3}{r['letter']}  {r['coords'][0]:<22} {r['symbol']}  {'' if r['ok'] else '(letter/setting mismatch)'}")
