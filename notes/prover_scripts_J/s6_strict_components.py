# On the strict survivors (all a_i > 0 at some legal beta; closed under in-twin moves since a is preserved),
# compare D-classes with connected components of in-twin switching (single moves).
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import collections, pickle
n=6
K=keys(); pos={k:i for i,k in enumerate(K)}
surv=[int(l) for l in open('/home/user/CLTN/certificates/n6_survivors.idx') if l.strip()]
data={}
for g in surv:
    A=graph_of_idx(g)
    bd=fastregion.beta_data(A,list(range(n)))
    r=fastregion.region(A,list(range(n)),None,bd=bd,positivity_only=True)
    if r['verdict']!='positivity_possible': continue
    data[g]=(tuple(bd['d']), tuple(tuple(p) for p in bd['phi']), tuple(tuple(p) for p in bd['omega']), tuple(sorted(tuple(p) for p in bd['a'])))
strict=sorted(data)
print('strict survivors', len(strict))
dcls=collections.defaultdict(list)
for g in strict: dcls[data[g]].append(g)
print('D-classes on strict set', len(dcls), 'max size', max(len(c) for c in dcls.values()))
parent={g:g for g in strict}
def find(a):
    while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
    return a
def in_twins(A):
    return [(u,v) for u in range(n) for v in range(u+1,n) if A[u][v]==A[v][u] and all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v))]
leave=0; cross=0; moves=0
for g in strict:
    A=graph_of_idx(g)
    for (u,v) in in_twins(A):
        for w in range(n):
            if w in (u,v) or A[u][w]==A[v][w]: continue
            B=[r[:] for r in A]; B[u][w],B[v][w]=A[v][w],A[u][w]
            h=pos[enum_cache.key_to_int(graphs.canonical_form(B)[0])]; moves+=1
            if h not in data: leave+=1; continue
            if data[h]!=data[g]: cross+=1
            a,b=find(g),find(h)
            if a!=b: parent[a]=b
comp=collections.defaultdict(list)
for g in strict: comp[find(g)].append(g)
print('moves',moves,'leaving strict set',leave,'crossing D-classes',cross,'components',len(comp))
split=[]
for key,c in dcls.items():
    roots={find(g) for g in c}
    if len(roots)>1: split.append((len(c),len(roots),c[0]))
split.sort(reverse=True)
print('D-classes with >1 twin component:', len(split)); print(split[:25])
print('size distribution of D-classes on strict set', sorted(collections.Counter(len(c) for c in dcls.values()).items()))
pickle.dump({'data':data,'comp':dict(comp),'dcls':{i:c for i,c in enumerate(dcls.values())},'split':split}, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','wb'))
