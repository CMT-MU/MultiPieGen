"""Left-coset decomposition of the parent's operations (ITA numbers, modulo the parent's lattice)
with respect to the image A = [f(1),...,f(n)] of the child's operations:
    S/T_S = A + g2*A + g3*A + ...      row_k = [ number of (g_k o a_i) for a_i in A ]   (a_i applied first)
g_k = smallest parent operation number not yet covered, so every row starts with its own minimum
and row 1 is A itself.  Inputs: parent_op_correspondence.csv, ita_ops.json"""
import csv, json, collections
import numpy as np

ITA = json.load(open('ita_ops.json'))
corr = collections.defaultdict(list)
for r in csv.DictReader(open('parent_op_correspondence.csv')): corr[int(r['G'])].append(r)
sym = {int(r['G']): r['G_symbol'] for r in csv.DictReader(open('parent_op_correspondence.csv'))}
psym = {int(r['G']): r['parent_symbol'] for r in csv.DictReader(open('parent_op_correspondence.csv'))}

from fractions import Fraction
trans = {int(r['G']): r for r in csv.DictReader(open('parent_transformations.csv'))}
def to_num(x):                                   # exact: all entries are dyadic rationals, so float() is exact
    f = Fraction(x); v = float(f); assert Fraction(v) == f
    return int(v) if v == int(v) else v
def Md(G):
    T = trans[G]
    return ([[to_num(T[f'M{i}{j}']) for j in range(1, 4)] for i in range(1, 4)], [to_num(T[f'd{i}']) for i in (1, 2, 3)])

key = lambda R: tuple(np.array(R, dtype=int).flatten())
out = {}; stats = collections.Counter()
for G, rows in sorted(corr.items()):
    S = int(rows[0]['parent'])
    Rs = [np.array(R, dtype=int) for R, t in ITA[str(S)]['ops']]          # parent ops, ITA order (linear parts)
    num = {key(R): j + 1 for j, R in enumerate(Rs)}
    assert len(num) == len(Rs)
    A = [int(r['parent_op_no']) for r in sorted(rows, key=lambda r: int(r['G_op_no']))]
    assert len(set(A)) == len(A) and A[0] == 1
    mul = lambda g, a: num[key(Rs[g - 1] @ Rs[a - 1])]                    # apply a first, then g
    assert {mul(a, b) for a in A for b in A} == set(A)                    # A is a subgroup of the parent's point group
    covered = set(A); cosets = [A]
    while len(covered) < len(Rs):
        g = min(set(range(1, len(Rs) + 1)) - covered)
        row = [mul(g, a) for a in A]
        assert not (set(row) & covered) and len(set(row)) == len(A) and row[0] == g
        cosets.append(row); covered |= set(row)
    assert sorted(x for c in cosets for x in c) == list(range(1, len(Rs) + 1))
    assert len(cosets) * len(A) == len(Rs)
    M, d = Md(G)
    Md4 = [M[i] + [d[i]] for i in range(3)] + [[0, 0, 0, 1]]
    out[G] = dict(parent=S, cosets=cosets, M=M, d=d, Md=Md4)
    stats[len(cosets)] += 1

L = ["# Left-coset decomposition of the parent's symmetry operations (ITA numbers, modulo the parent's lattice translations)",
     "# by the image of the child's operations.  cosets[0] = [parent op number of child op 1, 2, ...] (as in parent_op_correspondence);",
     "# cosets[k] = [number of g_k o (cosets[0][i]) for each i], g_k = smallest operation number not yet covered (a[i] applied first).",
     "# Row elements are in the child's operation order.  Parent lattice translations (index m of k-type pairs) are not part of this table.",
     "# Md: the same transformation as a 4x4 augmented matrix [[M, d], [0, 0, 0, 1]] acting on (x, y, z, 1)^T, like the symmetry operations.",
     "#     parent op (R'|t') = Md (R|t) Md^-1  modulo the parent's lattice translations.",
     "# M, d: x_parent = M x_child + d  (ITA standard settings, origin choice 2; M acts on column vectors of fractional coordinates).",
     "", "data_parent_coset_decomposition = {"]
for G, v in out.items():
    L.append(f"    {G}: {{'parent': {v['parent']}, 'cosets': {v['cosets']}, 'M': {v['M']}, 'd': {v['d']}, 'Md': {v['Md']}}},  # {sym[G]} -> {psym[G]}")
L.append("}\n")
open('parent_coset_decomposition.py', 'w').write("\n".join(L))
with open('parent_coset_decomposition.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['G', 'G_symbol', 'parent', 'parent_symbol', 'n_cosets', 'coset_no', 'ops', 'M', 'd', 'Md'])
    for G, v in out.items():
        for k, c in enumerate(v['cosets'], 1): w.writerow([G, sym[G], v['parent'], psym[G], len(v['cosets']), k, ' '.join(map(str, c)), str(v['M']), str(v['d']), str(v['Md'])])
print('groups:', len(out), ' number of cosets -> groups:', dict(sorted(stats.items())))
for G in (3, 14, 62, 123, 166, 215):
    print(G, sym[G], '->', out[G]['parent'], psym[G], out[G]['cosets'] if len(str(out[G]['cosets'])) < 150 else str(out[G]['cosets'])[:150] + '…')
