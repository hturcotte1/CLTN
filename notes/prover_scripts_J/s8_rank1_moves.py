# General rank-one moves N' = N + y x^T with x in {0,+-1}^n, y in {0,+-1}^n, N' a valid nonedge matrix,
# x ⊥ N^k 1 and x ⊥ N^k y for all k (Lemma: then D is preserved, a as a vector).
# Exact integer arithmetic (numpy int64 on values <= 5^6; no floats anywhere).
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import numpy as np, pickle, collections, time
n=6
S=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','rb'))
data=S['data']; strict=sorted(data)
K=keys(); pos={k:i for i,k in enumerate(K)}
cand=np.array(list(itertools.product([-1,0,1],repeat=n)),dtype=np.int64)
cand=cand[np.any(cand!=0,axis=1)]
# x up to sign: first nonzero entry positive
first=np.array([r[np.nonzero(r)[0][0]] for r in cand]); candx=cand[first>0]
def krylov(M,v):
    out=[v]
    for _ in range(n-1): out.append(M@out[-1])
    return np.array(out)
parent={g:g for g in strict}
def find(a):
    while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
    return a
stats=collections.Counter(); t0=time.time()
twinlike=collections.Counter()
for cnt,g in enumerate(strict):
    A=np.array(graph_of_idx(g),dtype=np.int64)
    N=np.ones((n,n),dtype=np.int64)-np.eye(n,dtype=np.int64)-A.T
    V=krylov(N,np.ones(n,dtype=np.int64))
    xs=candx[np.all(candx@V.T==0,axis=1)]
    for x in xs:
        # allowed y_w per row
        allowed=[]
        for w in range(n):
            al=[0]
            for s in (1,-1):
                if x[w]!=0: continue
                row=N[w]+s*x
                if np.all((row==0)|(row==1)): al.append(s)
            allowed.append(al)
        if all(len(a)==1 for a in allowed): continue
        Xk=krylov(N.T,x)          # rows: (N^T)^k x
        for y in itertools.product(*allowed):
            y=np.array(y,dtype=np.int64)
            if not y.any(): continue
            if np.any(Xk@y!=0): continue
            N2=N+np.outer(y,x)
            A2=(np.ones((n,n),dtype=np.int64)-np.eye(n,dtype=np.int64)-N2).T
            assert np.all((A2==0)|(A2==1)) and np.all(np.diag(A2)==0)
            h=pos[enum_cache.key_to_int(graphs.canonical_form(A2.tolist())[0])]
            stats['moves']+=1
            supp=int(np.count_nonzero(x)); twinlike[(supp,int(np.count_nonzero(y)))]+=1
            if h not in data: stats['leave']+=1; continue
            if data[h]!=data[g]: stats['cross']+=1; continue
            a,b=find(g),find(h)
            if a!=b: parent[a]=b
    if cnt%5000==0: print(cnt, round(time.time()-t0), dict(stats), flush=True)
comp=collections.defaultdict(list)
for g in strict: comp[find(g)].append(g)
dcls=collections.defaultdict(list)
for g in strict: dcls[data[g]].append(g)
split=[(len(c),len({find(g) for g in c}),c[0]) for c in dcls.values() if len({find(g) for g in c})>1]
split.sort(reverse=True)
print('stats',dict(stats)); print('|supp x|,|supp y| histogram', sorted(twinlike.items()))
print('components',len(comp),'D-classes',len(dcls),'split D-classes',len(split)); print(split[:20])
pickle.dump({'comp':dict(comp),'split':split}, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s8_out.pkl','wb'))
