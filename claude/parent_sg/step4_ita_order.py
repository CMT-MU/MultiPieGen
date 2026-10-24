import csv,ast,re,json,itertools,numpy as np,spglib
from fractions import Fraction
SP='.'
rows=list(csv.reader(open('/tmp/pyx/ex/pyxtal/database/wyckoff_list.csv')))
gp={int(r[0]):ast.literal_eval(r[1])[0] for r in rows[1:] if r[0].isdigit() and r[1]}
def parse(s):
    R=np.zeros((3,3),dtype=int); t=np.zeros(3)
    for i,comp in enumerate(s.replace(' ','').split(',')):
        for m in re.finditer(r'([+-]?)(\d+/\d+|\d+)?([xyz])?',comp):
            sign,num,var=m.groups()
            if not num and not var: continue
            sg=-1 if sign=='-' else 1
            if var:
                R[i,'xyz'.index(var)]+=sg*(1 if not num else 1)
            else:
                f=Fraction(num); t[i]+=sg*float(f)
    return R,t
first={}
for h in range(1,531): first.setdefault(spglib.get_spacegroup_type(h).number,h)
def modl(v): v=np.asarray(v)%1.0; v=np.round(v,6)%1.0; return tuple(np.round(v*1e5).astype(int)%100000)
ita={}; bad=[]
grid=[k/8 for k in range(8)]
for n in range(1,231):
    sym=spglib.get_symmetry_from_database(first[n])
    sops=list(zip(sym['rotations'],sym['translations']))
    cen=[t for r,t in sops if (r==np.eye(3,dtype=int)).all()]
    nc=len(cen); full=[parse(s) for s in gp[n]]
    npg=len(full)//nc
    reps=full[:npg]
    # check block structure: full[k*npg+i] == reps[i] + cen_k
    cenref=[modl(t) for t in cen]
    ok=len(full)==len(sops)
    sset={ (tuple(r.flatten()),modl(t)) for r,t in sops}
    def match(shift):
        # shifted set: conjugate by origin shift s:  (R,t)->(R,t+(I-R)s)
        S=set()
        for R,t in full:
            S.add((tuple(R.flatten()),modl(t+(np.eye(3)-R)@shift)))
        return S==sset
    found=None
    if match(np.zeros(3)): found=np.zeros(3)
    else:
        for s in itertools.product(grid,repeat=3):
            if match(np.array(s)): found=np.array(s);break
    if found is None: bad.append(n); continue
    cenita=[full[k*npg][1].tolist() for k in range(nc)]
    for k in range(nc):
        for i in range(npg):
            Rk,tk=full[k*npg+i]; assert (Rk==full[i][0]).all() and modl(tk-full[i][1]-np.array(cenita[k]))==modl(np.zeros(3)), (n,k,i)
    ita[n]=dict(shift=[float(x) for x in found],
        ops=[(R.tolist(),(t+(np.eye(3)-R)@found).tolist()) for R,t in reps], ncen=nc, cen=cenita)
print('unmatched',bad)
print('nonzero shifts',[(n,v['shift']) for n,v in ita.items() if any(v['shift'])])
json.dump({str(k):v for k,v in ita.items()},open(f'{SP}/ita_ops.json','w'))
