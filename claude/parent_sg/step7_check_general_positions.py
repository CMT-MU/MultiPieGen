"""Check that, for every child group G with a parent S, the i-th general position of G
equals the parent's general position f(i) (+ centering) pulled back to G's coordinates.
Exact rational arithmetic.  Inputs: parent_op_correspondence.csv, parent_transformations.csv, ita_ops.json"""
import csv, json, re, collections
from fractions import Fraction
import sympy as sp

def parse(s):
    R = sp.zeros(3, 3); t = sp.zeros(3, 1)
    for i, comp in enumerate(s.replace(' ', '').split(',')):
        for sign, num, var in re.findall(r'([+-]?)(\d+/\d+|\d+)?([xyz])?', comp):
            if not num and not var: continue
            sg = -1 if sign == '-' else 1
            if var: R[i, 'xyz'.index(var)] += sg
            else:   t[i] += sg * sp.Rational(num)
    return R, t

def fmt(R, t):
    out = []
    for i in range(3):
        s = ''
        for j, v in enumerate('xyz'):
            c = R[i, j]
            if c: s += ('+' if c > 0 and s else '-' if c < 0 else '') + (v if abs(c) == 1 else f'{abs(c)}{v}')
        if t[i] != 0: s += ('+' if t[i] > 0 and s else '') + str(t[i])
        out.append(s or '0')
    return ','.join(out)

rat = lambda x: sp.Rational(x)
isint = lambda v: all(c.is_integer for c in v)
ITA = json.load(open('ita_ops.json'))
corr = collections.defaultdict(list)
for r in csv.DictReader(open('parent_op_correspondence.csv')): corr[int(r['G'])].append(r)
trans = {int(r['G']): r for r in csv.DictReader(open('parent_transformations.csv'))}

res = {}; bad_lin = []; bad_lat = []
for G, rows in sorted(corr.items()):
    T = trans[G]
    M = sp.Matrix(3, 3, [rat(T[f'M{i}{j}']) for i in range(1, 4) for j in range(1, 4)])
    d = sp.Matrix([rat(T[f'd{i}']) for i in (1, 2, 3)]); Mi = M.inv()
    # child's full translation lattice: Z^3 + its centering vectors (ITA order)
    cenG = [sp.Matrix([sp.nsimplify(x, rational=True, tolerance=1e-9) for x in c]) for c in ITA[str(G)]['cen']]
    strict = True; mono = True; prev = (0, -1)
    for r in rows:
        Ri, ti = parse(r['G_op']); Rp, tp = parse(r['parent_op'])
        cen = sp.Matrix([rat(x) for x in r['parent_centering_vector'].split()])
        Rb = Mi * Rp * M
        tb = Mi * (Rp * d + tp + cen - d)          # parent op (with centering) in G's coordinates
        if Rb != Ri: bad_lin.append((G, r['G_op_no'])); continue
        diff = tb - ti
        if not isint(M * diff): bad_lat.append((G, r['G_op_no'])); continue   # must lie in M^-1 Z^3
        if not any(isint(diff - c) for c in cenG): strict = False          # not a translation of G itself
        key = (int(r['parent_op_no']), int(r['parent_centering_index']))
        if key < prev: mono = False
        prev = key
    res[G] = (strict, mono)
n = len(res)
print('groups checked:', n)
print('linear-part mismatches:', len(bad_lin), bad_lin[:5])
print('translation outside parent lattice:', len(bad_lat), bad_lat[:5])
print('all G ops identical to parent ops incl. translation (mod G lattice):', sum(v[0] for v in res.values()))
print('  of which differ only by a parent-lattice translation not in G:', [g for g, v in res.items() if not v[0]][:60])
print('parent op numbers increasing with child op number:', sum(v[1] for v in res.values()), 'of', n)
print('  non-monotone groups:', [g for g, v in res.items() if not v[1]][:40])
