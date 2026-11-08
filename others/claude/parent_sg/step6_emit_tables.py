import json,csv,sys,numpy as np,spglib,ast
SPD='.'
sys.argv=['x','/dev/null']
exec(open(f'{SPD}/step5_embeddings.py').read().split("if __name__")[0])
emb={int(k):v for k,v in json.load(open(f'{SPD}/emb_all.json')).items()}
rows=list(csv.reader(open('/tmp/pyx/ex/pyxtal/database/wyckoff_list.csv')))
gp={int(r[0]):ast.literal_eval(r[1])[0] for r in rows[1:] if r[0].isdigit() and r[1]}
sym=lambda n: spglib.get_spacegroup_type(first[n]).international_short
def xyz(n,i): return gp[n][i].replace(' ','')
def fr(x):
    from fractions import Fraction
    f=Fraction(float(x)).limit_denominator(48); return str(f)
nbad=0; table=[]; opmap={}; trans=[]
for G,v in sorted(emb.items()):
    S=v['S']; M=np.array(v['M']); d=np.array(v['d'])
    r=verify(G,S,M,d)            # independent re-verification
    assert r is not None and r[0]==v['f'] and r[1]==v['cen'],G
    f,cs=r
    assert len(set(zip(f,cs)))==len(f) and f[0]==1 and cs[0]==0
    sG=np.array(ITA[G]['shift']); sS=np.array(ITA[S]['shift'])
    dita=M@sG+d-sS
    cenS=ITA[S]['cen']
    for i,(j,c) in enumerate(zip(f,cs)):
        table.append([G,sym(G),i+1,xyz(G,i),S,sym(S),j,xyz(S,j-1),c,' '.join(fr(x) for x in cenS[c]),v['m']])
    opmap[G]=dict(parent=S,op=f,centering=cs)
    trans.append([G,sym(G),S,sym(S),v['m']]+[fr(x) for x in M.flatten()]+[fr(x) for x in dita])
out='/home/claude/parent_sg'
with open(f'{out}/parent_op_correspondence.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['G','G_symbol','G_op_no','G_op','parent','parent_symbol','parent_op_no','parent_op','parent_centering_index','parent_centering_vector','lattice_index']); w.writerows(table)
with open(f'{out}/parent_transformations.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['G','G_symbol','parent','parent_symbol','lattice_index']+[f'M{i}{j}' for i in range(1,4) for j in range(1,4)]+['d1','d2','d3']); w.writerows(trans)
L=["# For each space group G (ITA number): parent space group and the correspondence of ITA symmetry-operation numbers.",
"# op[i-1] = number (1-based, ITA order) of the parent operation to which operation i of G corresponds, modulo the parent's lattice translations.",
"# centering[i-1] = index of the parent's centering translation (ITA order, 0 = (0,0,0)) combined with that parent operation.",
"# Coordinates: ITA standard settings, origin choice 2 where two choices exist (inversion at the origin); hexagonal axes for rhombohedral groups.",
"# Embedding: x_parent = M x_G + d (see parent_transformations.csv).",
"","data_parent_op_correspondence = {"]
for G,v in opmap.items():
    L.append(f"    {G}: {{'parent': {v['parent']}, 'op': {v['op']}, 'centering': {v['centering']}}},  # {sym(G)} -> {sym(v['parent'])}")
L.append("}\n")
open(f'{out}/parent_op_correspondence.py','w').write("\n".join(L))
print(len(table),'rows; groups',len(opmap))
import collections
print('centering used:',collections.Counter(c for r in table for c in [r[8]]))
for G in (4,14,123,43):
    print(G,[ (r[2],r[3],r[6],r[7],r[8]) for r in table if r[0]==G][:6])
