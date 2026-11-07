import sys; sys.path.insert(0, ".")
import sys, pickle, multiprocessing as mp
import check_complete as C
def f(n):
    try: return n, C.check_group(n)
    except Exception as e: return n, dict(name=n, problems=[f"EXC {type(e).__name__}: {e}"], per_W={})
if __name__=="__main__":
    names=list(C.data)
    out={}
    with mp.Pool(2) as p:
        for i,(n,r) in enumerate(p.imap_unordered(f,names),1):
            out[n]=r
            if i%20==0: print(i,flush=True)
    pickle.dump(out,open("complete_all.pkl","wb"))
