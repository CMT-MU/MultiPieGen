import spglib, numpy as np, itertools, json, sys
from multiprocessing import Pool
SPD='.'
exec(open(f'{SPD}/common.py').read())
src=open(f'{SPD}/step2_enum_k_subgroups.py').read()
exec(src[src.index('def hnf_basis'):src.index('def run(')])   # hnf_basis, sublattices, gens_of

ITA={int(k):v for k,v in json.load(open(f'{SPD}/ita_ops.json')).items()}
CH={int(k):v for k,v in json.load(open(f'{SPD}/chosen.json')).items() if v}
D={int(k):v for k,v in json.load(open(f'{SPD}/subs.json')).items()}
I3=np.eye(3)

def prim_basis(n):
    cv=[np.round(np.array(c)*24).astype(int) for c in ITA[n]['cen']]
    vecs=[v*24 for v in np.eye(3,dtype=int)]+cv
    return np.array(hnf_basis([list(map(int,v)) for v in vecs]),dtype=float)/24.0

def ita_ops(n): return [(np.array(R),np.array(t)) for R,t in ITA[n]['ops']]

def embeddings(S,m,G):
    """yield (M0,d0): x_S = M0 x_G + d0 for each subgroup H of S (index m) of type G"""
    h=first[S]; sym=spglib.get_symmetry_from_database(h)
    R=sym['rotations']; T=sym['translations']; ops=list(zip(R,T))
    P={key(r):r for r in R}; L=lattice(S)
    BS=prim_basis(S)
    rng=np.random.default_rng(777+S); xs=[rng.random(3)*0.8+0.07 for _ in range(3)]
    qsize=len(ITA[G]['ops'])
    ropd={}
    for r,t in ops: ropd.setdefault(key(r),t)
    out=[]
    subs=[q for q in subgroups(P) if len(q)==qsize]
    for (dd,Hm) in sublattices(m):
        TB=Hm@BS; iTB=np.linalg.inv(TB)
        reps=[np.array(nn)@BS for nn in itertools.product(range(dd[0]),range(dd[1]),range(dd[2]))]
        def canon(t):
            c=np.round((t@iTB)%1.0,6)%1.0
            return tuple(np.round(c*1e6).astype(np.int64)%1000000)
        for q in subs:
            ok=True
            for k in q:
                Rm=np.array(k).reshape(3,3); c=(TB@Rm.T)@iTB
                if np.abs(c-np.round(c)).max()>1e-7: ok=False;break
            if not ok: continue
            gs=gens_of(q)
            lifts=[[(np.array(g),ropd[key(g)]+w) for w in reps] for g in gs]
            seen=set()
            for combo in itertools.product(*lifts):
                elems={key(np.eye(3,dtype=int)):np.zeros(3)}
                frontier=[(np.eye(3,dtype=int),np.zeros(3))]; bad=False
                while frontier and not bad:
                    new=[]
                    for (a,ta) in frontier:
                        for (g,tg) in combo:
                            Rn=g@a; tn=g@ta+tg; kk=key(Rn)
                            if kk in elems:
                                if canon(elems[kk])!=canon(tn): bad=True;break
                            else: elems[kk]=tn; new.append((Rn,tn))
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
                cellH=(TB@L,np.array(pos),num)
                ds0=spglib.get_symmetry_dataset(cellH,symprec=1e-4)
                if ds0 is None or ds0.number!=G: continue
                ds=spglib.get_symmetry_dataset(cellH,symprec=1e-4,hall_number=first[G])
                if ds is None: continue
                Pm=np.array(ds.transformation_matrix)@np.linalg.inv(TB).T
                M0=np.linalg.inv(Pm); d0=-M0@np.array(ds.origin_shift)
                out.append((M0,d0))
    # dedupe
    uniq={}
    for M0,d0 in out:
        uniq[(tuple(np.round(M0.flatten(),6)),tuple(np.round(d0%1.0,6)))]=(M0,d0)
    return list(uniq.values())

def verify(G,S,M,d,tol=1e-6):
    """apply x_S=M x_G+d to ITA ops of G; return (f, cen) or None"""
    opsG=ita_ops(G); opsS=ita_ops(S); BS=prim_basis(S); iBS=np.linalg.inv(BS)
    cenS=[np.array(c) for c in ITA[S]['cen']]
    Rs={tuple(R.flatten()):j for j,(R,t) in enumerate(opsS)}
    Mi=np.linalg.inv(M); f=[];cs=[]
    for R,t in opsG:
        Rp=M@R@Mi
        if np.abs(Rp-np.round(Rp)).max()>tol: return None
        kk=tuple(np.round(Rp).astype(int).flatten())
        if kk not in Rs: return None
        j=Rs[kk]; tp=M@t+d-Rp@d; diff=tp-opsS[j][1]
        c=diff@iBS
        if np.abs(c-np.round(c)).max()>tol: return None
        kc=None
        for k,cv in enumerate(cenS):
            e=diff-cv
            if np.abs(e-np.round(e)).max()<tol: kc=k;break
        if kc is None: return None
        f.append(j+1);cs.append(kc)
    BG=prim_basis(G)
    for v in BG:
        w=M@v; c=w@iBS
        if np.abs(c-np.round(c)).max()>tol: return None
    return f,cs

def normalizer_Ls(G):
    """candidate linear parts (conv coords) of affine normalizer elements of G: list of (L, perm)"""
    opsG=ita_ops(G); n=len(opsG); BG=prim_basis(G)
    Rs={tuple(R.flatten()):i for i,(R,t) in enumerate(opsG)}
    vals=np.array(list(itertools.product((-1,0,1),repeat=9)),dtype=float).reshape(-1,3,3)
    det=np.linalg.det(vals); U=vals[np.abs(np.abs(det)-1)<1e-9]
    Lm=np.einsum('ji,kjl,lm->kim',BG,U,np.linalg.inv(BG.T).T) if False else None
    BT=BG.T; BTi=np.linalg.inv(BT)
    Lall=BT@U@BTi   # (K,3,3), column convention
    Linv=np.linalg.inv(Lall)
    out=[]
    Rarr=np.array([R for R,t in opsG],dtype=float)
    conj=np.einsum('kab,nbc,kcd->knad',Lall,Rarr,Linv)
    rnd=np.round(conj)
    okint=(np.abs(conj-rnd).max(axis=(2,3))<1e-7).all(axis=1)
    for k in np.where(okint)[0]:
        perm=[]
        good=True
        for i in range(n):
            kk=tuple(rnd[k,i].astype(int).flatten())
            if kk not in Rs: good=False;break
            perm.append(Rs[kk])
        if good and len(set(perm))==n: out.append((Lall[k],perm))
    return out

_sgrid=None
def sol_s(G,L,perm,cache):
    key_=tuple(np.round(L.flatten(),6))
    if key_ in cache: return cache[key_]
    global _sgrid
    if _sgrid is None:
        g=np.array(list(itertools.product(range(24),repeat=3)),dtype=float)/24.0; _sgrid=g
    opsG=ita_ops(G); BG=prim_basis(G); iBG=np.linalg.inv(BG)
    mask=np.ones(len(_sgrid),bool)
    for i,(R,t) in enumerate(opsG):
        Rp=opsG[perm[i]][0]; tg=opsG[perm[i]][1]
        v=(L@t-tg)[None,:]+_sgrid@(I3-Rp).T
        c=v@iBG
        mask&=(np.abs(c-np.round(c)).max(axis=1)<1e-7)
        if not mask.any(): break
    res=_sgrid[mask]
    cache[key_]=res
    return res

def best_embedding(G):
    S,m=CH[G]
    A0s=embeddings(S,m,G)
    if not A0s: return G,None
    Ls=normalizer_Ls(G)
    cache={}
    BSS=prim_basis(S); cenS=[np.array(c) for c in ITA[S]['cen']]
    def redd(d):
        best=None
        for cv in cenS:
            e=d-cv; e=e-np.round(e)   # reduce mod Z^3
            c=np.abs(e).sum()
            if best is None or c<best[0]-1e-12: best=(c,e)
        return best
    cand=[]
    for (M0,d0) in A0s:
        Ms=[M0@L for L,_ in Ls]
        c1=np.array([np.abs(M-I3).sum() for M in Ms])
        order=np.argsort(c1,kind='stable')
        bestc1=None
        for idx in order:
            if bestc1 is not None and c1[idx]>bestc1+1e-9: break
            L,perm=Ls[idx]
            sols=sol_s(G,L,perm,cache)
            if len(sols)==0: continue
            bestc1=c1[idx]
            M=M0@L
            dall=sols@M0.T+d0
            best=None
            E=dall[:,None,:]-np.array(cenS)[None,:,:]; E=E-np.round(E)
            cdall=np.abs(E).sum(axis=2); kmin=cdall.argmin(axis=1)
            cdm=cdall[np.arange(len(sols)),kmin]
            order2=np.argsort(cdm,kind='stable')
            tried=0
            for ii in order2:
                if tried>=12 and cdm[ii]>cdm[order2[0]]+1e-9: break
                dred=E[ii,kmin[ii]]
                r=verify(G,S,M,dred)
                tried+=1
                if r is None: continue
                f,cs=r
                cand.append(((round(c1[idx],6),round(float(cdm[ii]),6),f),M,dred,f,cs))
    if not cand: return G,None
    cand.sort(key=lambda x:x[0])
    k,M,d,f,cs=cand[0]
    return G,dict(S=S,m=m,M=M.tolist(),d=d.tolist(),f=f,cen=cs,nA0=len(A0s),nL=len(Ls),ncand=len(cand))

if __name__=='__main__':
    Gs=[int(x) for x in sys.argv[2:]] if len(sys.argv)>2 else sorted(CH)
    with Pool(2) as pool:
        res=dict(pool.imap_unordered(best_embedding,Gs))
    json.dump({str(k):v for k,v in res.items()},open(sys.argv[1],'w'))
    print('done',sum(v is not None for v in res.values()),'/',len(res))
