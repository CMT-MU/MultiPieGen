import sys; sys.path.insert(0,".")
import pickle, itertools, collections, multiprocessing as mp
import sympy as sp
import check_complete as C
S=C.S
res=pickle.load(open("complete_all.pkl","rb"))
todo=[(n,W) for n,v in res.items() for W,p in v["per_W"].items() if p["missing"]]
def work(n):
    d=C.data[n]; v=res[n]
    ops=C.rebuild_ops(d)[1]; space="^" in n
    cen=v.get("centering"); shifts=S._centering_shifts(cen)
    out=[]
    for (nn,W) in todo:
        if nn!=n: continue
        ents=[(l,e) for l,e in d.items() if l.split("@")[1]==W]
        c=S._parse(ents[0][1][1])[0]; cM=sp.Matrix(c)
        params=[p for p in (C.x,C.y,C.z) if cM.has(p)]
        r=S.site_symmetry(ops,[list(c)],lattice_mod=space,centering=cen,oriented_errors_ok=True)
        H=list({tuple(ops[i-1][:3,:3]):ops[i-1][:3,:3] for i in r.stabilizer}.values())
        extra=[]
        A=sp.Matrix([[sp.diff(e,p) for p in params] for e in c]); bb=sp.Matrix([e.subs({p:0 for p in params}) for e in c])
        cands=sorted({sp.Rational(k,m) for m in (1,2,3,4,6,8,12) for k in range(m)})
        for M in ops:
            R,t=M[:3,:3],M[:3,3]
            RA=R*A
            for sg in itertools.product([1,-1],repeat=len(params)):
                if RA!=A*sp.diag(*sg): continue
                r=R*bb+t-bb
                ok=any(all(((A*sp.Matrix(sv)-r)[i]).is_integer for i in range(3)) for sv in itertools.product(cands,repeat=len(params)))
                if ok: extra.append(R); break
        Rs=list({tuple(R):R for R in H+extra}.values())
        # closure under products
        changed=True
        while changed:
            changed=False
            for a in list(Rs):
                for b in list(Rs):
                    m=a*b
                    if tuple(m) not in {tuple(q) for q in Rs}: Rs.append(m); changed=True
        classes=C.vector_strata(H)
        # merge classes under Rs
        listed=set()
        for l,(vs,cs,m) in ents:
            ce=S._parse(cs)[0]
            g=None
            for M in ops:
                if S._key(M[:3,:3]*cM+M[:3,3],space,shifts)==S._key(sp.Matrix(ce),space,shifts): g=M[:3,:3];break
            listed.add(C.image(g.inv(),C.entry_subspace(S._parse(vs)[0])))
        reach={C.image(R,s) for s in listed for R in Rs}
        miss=[C.form(min(cl,key=lambda s:(len(s),str(s)))) for cl in classes if not (cl & reach)]
        out.append((n,W,len(extra),miss))
    return out
if __name__=="__main__":
    names=sorted({n for n,_ in todo})
    with mp.Pool(2) as p: rr=p.map(work,names)
    for o in rr:
        for x in o: print(x)
