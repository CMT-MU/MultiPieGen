import json,collections,csv,spglib
SP='.'
d={int(k):v for k,v in json.load(open(f'{SP}/subs.json')).items()}
res={int(k):v for k,v in json.load(open(f'{SP}/res_m4.json')).items()}
PG=json.load(open(f'{SP}/pgmap.json'))
c=collections.defaultdict(set)
for S,lst in res.items():
    for g,m,q in lst:
        if (g!=S or m>1) and PG.get(d[g]['pg'])==d[S]['pg']: c[g].add((m,S))
first={}
for h in range(1,531): first.setdefault(spglib.get_spacegroup_type(h).number,h)
sym=lambda n: spglib.get_spacegroup_type(first[n]).international_short
rows=[];chosen={}
for G in range(1,231):
    pg=d[G]['pg']
    if pg not in PG: chosen[G]=None; rows.append([G,sym(G),pg,'-','-','-','','']); continue
    L=sorted(c[G]); mmin=L[0][0]
    at=[S for m,S in L if m==mmin]
    chosen[G]=(at[0],mmin)
    rows.append([G,sym(G),pg,PG[pg],at[0],mmin,len(at),' '.join(f"{S}({m})" for m,S in L)])
json.dump({str(k):v for k,v in chosen.items()},open(f'{SP}/chosen.json','w'))
with open(f'{SP}/candidates.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(['G','G_symbol','G_pg','parent_pg','chosen_parent','lattice_index','n_tied_at_that_index','all_candidates(S(index))']);w.writerows(rows)
cnt=collections.Counter(v[1] for v in chosen.values() if v)
print('index distribution',dict(cnt))
print('ties:',sum(1 for r in rows if r[6]!='' and r[6]>1),'unique:',sum(1 for r in rows if r[6]==1))
for G in (1,2,3,4,5,7,14,15,19,43,51,62,76,80,88,144,151,166,169,205,212,213,215):
    r=rows[G-1]; print(r[:7], sym(r[4]) if r[4] not in ('-',) else '')
