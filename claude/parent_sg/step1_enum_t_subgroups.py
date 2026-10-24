import spglib, numpy as np, itertools, json, sys
rng=np.random.default_rng(12345)

first={}
for h in range(1,531):
    t=spglib.get_spacegroup_type(h)
    first.setdefault(t.number,h)

def lattice(n):
    if n<=2:  return np.array([[5.0,0,0],[1.1,6.3,0],[0.7,1.9,7.9]])
    if n<=15:
        b=np.radians(100.3); return np.array([[5.1,0,0],[0,6.2,0],[7.3*np.cos(b),0,7.3*np.sin(b)]])
    if n<=74: return np.diag([5.1,6.2,7.3])
    if n<=142: return np.diag([5.1,5.1,7.3])
    if n<=194:
        a=5.1; return np.array([[a,0,0],[-a/2,a*np.sqrt(3)/2,0],[0,0,7.3]])
    return np.diag([5.4,5.4,5.4])

def key(R): return tuple(int(x) for x in R.flatten())
def mul(A,B): return A@B

def closure(gens):
    I=np.eye(3,dtype=int)
    S={key(I):I}
    frontier=[I]
    while frontier:
        new=[]
        for a in frontier:
            for g in gens:
                c=g@a
                k=key(c)
                if k not in S: S[k]=c; new.append(c)
        frontier=new
    return frozenset(S.keys())

def subgroups(P):  # P: dict key->matrix
    mats=list(P.values())
    cyc={closure([m]) for m in mats}
    allsg=set(cyc)
    changed=True
    while changed:
        changed=False
        cur=list(allsg)
        for a in cur:
            for c in cyc:
                if c<=a: continue
                gens=[np.array(k).reshape(3,3) for k in a|c]
                j=closure(gens)
                if j not in allsg: allsg.add(j); changed=True
    return allsg

LABEL_FIX={}
def pg_label(n,sch):
    if sch=='C3i': sch='C3i'
    if sch=='S6': sch='C3i'
    m1 = {149,151,153}   # D3-1 (312 type)
    if sch=='D3'  and n in m1: return 'D3-1'
    if sch=='C3v' and n in {157,159}: return 'C3v-1'
    if sch=='D3d' and n in {162,163}: return 'D3d-1'
    if sch=='D3h' and n in {189,190}: return 'D3h-1'   # label unused for parent lookup
    if sch=='D2d' and n in range(115,121): return 'D2d-1'  # -4m2 type (115-120)
    return sch

def species_positions(ops,xs):
    pos=[];num=[]
    for sp,x in enumerate(xs,1):
        seen=[]
        for R,t in ops:
            p=(R@x+t)%1.0
            seen.append(p)
        pos+= [list(p) for p in seen]; num+=[sp]*len(seen)
    return np.array(pos),num

out={}
for n in range(1,231):
    h=first[n]; ty=spglib.get_spacegroup_type(h)
    sym=spglib.get_symmetry_from_database(h)
    R=sym['rotations']; T=sym['translations']
    ops=list(zip(R,T))
    P={key(r):r for r in R}
    lab=pg_label(n,ty.pointgroup_schoenflies)
    L=lattice(n)
    xs=[rng.random(3)*0.8+0.07 for _ in range(3)]
    subs=subgroups(P)
    res={}
    for q in subs:
        hops=[(r,t) for r,t in ops if key(r) in q]
        pos,num=species_positions(hops,xs)
        ds=spglib.get_symmetry_dataset((L,pos,num),symprec=1e-4)
        res[tuple(sorted(q))]=ds.number
    out[n]=dict(pg=lab,ordpg=len(P),H=ty.international_short,
                subs=[(len(k),v) for k,v in res.items()])
    print(n,ty.international_short,lab,len(subs),file=sys.stderr)
json.dump({n:dict(pg=v['pg'],H=v['H'],subs=v['subs']) for n,v in out.items()},open(sys.argv[1],'w'))
