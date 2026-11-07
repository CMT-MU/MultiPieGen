import sys; sys.path.insert(0,".")
import pickle, itertools, sympy as sp
import check_complete as C
S=C.S
res=pickle.load(open("complete_all.pkl","rb"))
B=[(n,W) for n,v in res.items() for W,p in v["per_W"].items() if p["missing"]]
B.sort()
import csv
rows=[]
for n,W in B:
    d=C.data[n]; v=res[n]; ops=C.rebuild_ops(d)[1]; cen=v["centering"]; shifts=S._centering_shifts(cen) or [sp.zeros(3,1)]
    ents=[(l,e) for l,e in d.items() if l.split("@")[1]==W]
    c=S._parse(ents[0][1][1])[0]; cM=sp.Matrix(c)
    params=[p for p in (C.x,C.y,C.z) if cM.has(p)]
    r=S.site_symmetry(ops,[list(c)],lattice_mod=True,centering=cen,oriented_errors_ok=True)
    H=list({tuple(ops[i-1][:3,:3]):ops[i-1][:3,:3] for i in r.stabilizer}.values())
    classes=C.vector_strata(H)
    A=sp.Matrix([[sp.diff(e,p) for p in params] for e in c]); bb=sp.Matrix([e.subs({p:0 for p in params}) for e in c])
    k=len(params)
    perms=[]
    for pm in itertools.permutations(range(k)):
        for sg in itertools.product([1,-1],repeat=k):
            Sg=sp.zeros(k,k)
            for i in range(k): Sg[pm[i],i]=sg[i]
            perms.append(Sg)
    cands=sorted({sp.Rational(a,m) for m in (1,2,3,4,6,8,12) for a in range(m)})
    ext=[]
    for M in ops:
        for sh in shifts:
            R=M[:3,:3]; t=M[:3,3]+sh
            for Sg in perms:
                if R*A!=A*Sg: continue
                rr=R*bb+t-A*0-bb
                if any(all(((A*sp.Matrix(sv)-rr)[i]).is_integer for i in range(3)) for sv in itertools.product(cands,repeat=k)):
                    ext.append(R); break
    ext=list({tuple(R):R for R in ext}.values())
    # entry subspaces at the representative site
    lab_of={}
    for l,(vs,cs,m) in ents:
        ce=S._parse(cs)[0]
        g=next(M[:3,:3] for M in ops if S._key(M[:3,:3]*cM+M[:3,3],True,S._centering_shifts(cen))==S._key(sp.Matrix(ce),True,S._centering_shifts(cen)))
        sub=C.image(g.inv(),C.entry_subspace(S._parse(vs)[0]))
        for ci,cl in enumerate(classes):
            if sub in cl: lab_of.setdefault(ci,[]).append(l)
    listed={ci for ci in lab_of}
    print(n,W,"extended ops",len(ext))
    for ci,cl in enumerate(classes):
        if ci in listed: continue
        s=min(cl,key=lambda s:(len(s),str(s)))
        hits=set()
        for R in ext:
            im=C.image(R,s)
            for cj,cl2 in enumerate(classes):
                if im in cl2 and cj in listed: hits.add((cj,tuple(lab_of[cj])))
        hl=sorted({l for cj,labs in hits for l in labs})
        print("   missing",C.form(s),"->",hl)
        rows.append((n,W,len(H),C.form(s),";".join(hl),"yes" if hl else "NO"))

w=csv.writer(open("equivalence_41.csv","w",newline=""))
w.writerow(["group","wyckoff","H_order","H-class_not_listed","equivalent_listed_entries","resolved"]); w.writerows(rows)
print(len(rows),"missing H-classes in",len(B),"positions; unresolved:",sum(r[-1]=="NO" for r in rows))
