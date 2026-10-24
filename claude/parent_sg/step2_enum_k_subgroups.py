import spglib, numpy as np, itertools, json, sys, math
from multiprocessing import Pool
SPD='.'
exec(open(f'{SPD}/common.py').read())

def hnf_basis(vecs):
    A=[list(v) for v in vecs]
    basis=[]
    for col in range(3):
        while True:
            nz=[r for r in A if r[col]!=0]
            if len(nz)<=1: break
            nz.sort(key=lambda r:abs(r[col]))
            p=nz[0]
            for r in nz[1:]:
                q=r[col]//p[col]
                for j in range(3): r[j]-=q*p[j]
        nz=[r for r in A if r[col]!=0]
        if nz:
            p=nz[0]; A.remove(p); basis.append(p)
    return basis

def sublattices(m):
    out=[]
    for d1 in range(1,m+1):
        if m%d1: continue
        for d2 in range(1,m//d1+1):
            if (m//d1)%d2: continue
            d3=m//d1//d2
            for a in range(d2):
              for b in range(d3):
                for c in range(d3):
                    out.append(((d1,d2,d3),np.array([[d1,a,b],[0,d2,c],[0,0,d3]])))
    return out

def gens_of(Q):
    mats=[np.array(k).reshape(3,3) for k in Q]
    gens=[];cur=closure([])
    while cur!=Q:
        best=None
        for g in mats:
            c=closure(gens+[g])
            if best is None or len(c)>len(best[1]): best=(g,c)
        gens.append(best[0]); cur=best[1]
    return gens

def run(S, mmax, sizes_allowed):
    h=first[S]; ty=spglib.get_spacegroup_type(h)
    sym=spglib.get_symmetry_from_database(h)
    R=sym['rotations']; T=sym['translations']
    ops=list(zip(R,T)); P={key(r):r for r in R}
    L=lattice(S)
    cvecs=[np.round(t*12).astype(int) for r,t in ops if (r==np.eye(3,dtype=int)).all()]
    vecs=[np.array(v)*12 for v in np.eye(3,dtype=int)]+cvecs
    B=np.array(hnf_basis([list(map(int,v)) for v in vecs]),dtype=float)/12.0
    rng=np.random.default_rng(777+S)
    xs=[rng.random(3)*0.8+0.07 for _ in range(3)]
    subs=[q for q in subgroups(P) if len(q) in sizes_allowed]
    ropd={}
    for r,t in ops: ropd.setdefault(key(r),t)
    res=[]
    for m in range(1,mmax+1):
        for (dd,Hm) in sublattices(m):
            TB=Hm@B; iTB=np.linalg.inv(TB)
            reps=[np.array(n)@B for n in itertools.product(range(dd[0]),range(dd[1]),range(dd[2]))]
            def canon(t):
                c=(t@iTB)%1.0
                c=np.round(c,6)%1.0
                return tuple(np.round(c*1e6).astype(np.int64)%1000000)
            for q in subs:
                ok=True
                for k in q:
                    Rm=np.array(k).reshape(3,3)
                    c=(TB@Rm.T)@iTB
                    if np.abs(c-np.round(c)).max()>1e-7: ok=False;break
                if not ok: continue
                gs=gens_of(q)
                lifts=[[(np.array(g),ropd[key(g)]+w) for w in reps] for g in gs]
                seen=set()
                for combo in itertools.product(*lifts):
                    # closure mod T'
                    elems={ (key(np.eye(3,dtype=int))): np.zeros(3) }
                    if canon(np.zeros(3))!=canon(elems[key(np.eye(3,dtype=int))]): continue
                    frontier=[(np.eye(3,dtype=int),np.zeros(3))]
                    bad=False
                    while frontier and not bad:
                        new=[]
                        for (a,ta) in frontier:
                            for (g,tg) in combo:
                                Rn=g@a; tn=g@ta+tg
                                kk=key(Rn)
                                if kk in elems:
                                    if canon(elems[kk])!=canon(tn): bad=True;break
                                else:
                                    elems[kk]=tn; new.append((Rn,tn))
                            if bad: break
                        frontier=new
                    if bad or len(elems)!=len(q): continue
                    sig=tuple(sorted((k,canon(t)) for k,t in elems.items()))
                    if sig in seen: continue
                    seen.add(sig)
                    pos=[];num=[]
                    for sp,x in enumerate(xs,1):
                        for k,t in elems.items():
                            p=(np.array(k).reshape(3,3)@x+t)@iTB
                            pos.append(list(p%1.0)); num.append(sp)
                    ds=spglib.get_symmetry_dataset((TB@L,np.array(pos),num),symprec=1e-4)
                    res.append((ds.number,m,len(q)))
    return S,sorted(set(res))

def work(a): return run(*a)
if __name__=='__main__':
    mmax=int(sys.argv[1]); outp=sys.argv[2]
    d={int(k):v for k,v in json.load(open(f'{SPD}/subs.json')).items()}
    pgmap=PGMAP
    order={}
    for S,v in d.items(): order[v['pg']]=max(k for k,_ in v['subs'])
    jobs=[]
    for S in range(1,231):
        par=d[S]['pg']
        sizes={order[g] for g,p in pgmap.items() if p==par}
        if sizes: jobs.append((S,mmax,sizes))
    with Pool(2) as pool:
        out=dict(pool.imap_unordered(work,jobs))
    json.dump(out,open(outp,'w'))
