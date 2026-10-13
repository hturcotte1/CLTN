# Union-find over the 34,071 survivors using single in-twin switching moves.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import time, collections
n=6
K = keys()
pos = {k:i for i,k in enumerate(K)}
surv = [int(l) for l in open('/home/user/CLTN/certificates/n6_survivors.idx') if l.strip()]
survset = set(surv)
cls = classes()
cls_of = {}
for ci,c in enumerate(cls):
    for g in c: cls_of[g]=ci
assert len(cls_of)==len(surv)==34071
parent = {g:g for g in surv}
def find(a):
    while parent[a]!=a:
        parent[a]=parent[parent[a]]; a=parent[a]
    return a
def union(a,b):
    a,b=find(a),find(b)
    if a!=b: parent[a]=b
def in_twins(A):
    out=[]
    for u in range(n):
        for v in range(u+1,n):
            if A[u][v]!=A[v][u]: continue
            if all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v)):
                out.append((u,v))
    return out
t0=time.time()
nmoves=0; outside=0; crossclass=0
for cnt,g in enumerate(surv):
    A = graph_of_idx(g)
    for (u,v) in in_twins(A):
        for w in range(n):
            if w in (u,v) or A[u][w]==A[v][w]: continue
            B=[r[:] for r in A]; B[u][w],B[v][w]=A[v][w],A[u][w]
            h = pos[enum_cache.key_to_int(graphs.canonical_form(B)[0])]
            nmoves+=1
            if h not in survset: outside+=1; continue
            if cls_of[h]!=cls_of[g]: crossclass+=1
            union(g,h)
    if cnt%5000==0: print(cnt, time.time()-t0, flush=True)
comp=collections.defaultdict(list)
for g in surv: comp[find(g)].append(g)
comps=list(comp.values())
print('moves',nmoves,'leaving survivor set',outside,'crossing D-classes',crossclass)
print('components',len(comps),'D-classes',len(cls))
# compare partitions
same = all(len({cls_of[g] for g in c})==1 for c in comps)
print('each component inside one D-class:', same)
compsize_by_class = collections.Counter()
mism=[]
for ci,c in enumerate(cls):
    roots={find(g) for g in c}
    if len(roots)>1: mism.append((len(c),len(roots),ci))
print('D-classes split into >1 component:', len(mism))
mism.sort(reverse=True); print(mism[:30])
json.dump({'components':comps,'mismatch':mism}, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s3_out.json','w'))
