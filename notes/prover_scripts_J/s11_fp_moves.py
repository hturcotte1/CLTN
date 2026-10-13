# Task 3 in general: does each kind of D-preserving move preserve FP(G) and the stable supports?
# Exact: ctln.analyze_support on every support, at legal (eps, delta). Sampled over n=5 graphs (all 9,608 scanned
# for moves; up to CAP moves tested per kind) -- deterministic order.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import numpy as np, collections, random, sys
from itertools import combinations
n=int(sys.argv[1]) if len(sys.argv)>1 else 5
CAP=int(sys.argv[2]) if len(sys.argv)>2 else 150
keys_n=enum_cache.load_keys(n)
params=[(F(1,4),F(1,2)),(F(1,10),F(1,8))]
cand=np.array(list(itertools.product([-1,0,1],repeat=n)),dtype=np.int64); cand=cand[np.any(cand!=0,axis=1)]
first=np.array([r[np.nonzero(r)[0][0]] for r in cand]); candx=cand[first>0]
def krylov(M,v):
    out=[v]
    for _ in range(n-1): out.append(M@out[-1])
    return np.array(out)
def moves(A):
    A=np.array(A,dtype=np.int64); N=np.ones((n,n),dtype=np.int64)-np.eye(n,dtype=np.int64)-A.T
    V=krylov(N,np.ones(n,dtype=np.int64))
    for x in candx[np.all(candx@V.T==0,axis=1)]:
        allowed=[]
        for w in range(n):
            al=[0]
            for s in (1,-1):
                if x[w]==0 and np.all(((N[w]+s*x)==0)|((N[w]+s*x)==1)): al.append(s)
            allowed.append(al)
        Xk=krylov(N.T,x)
        for y in itertools.product(*allowed):
            y=np.array(y,dtype=np.int64)
            if not y.any() or np.any(Xk@y!=0): continue
            B=(np.ones((n,n),dtype=np.int64)-np.eye(n,dtype=np.int64)-(N+np.outer(y,x))).T
            supp=list(np.nonzero(x)[0])
            if len(supp)==2 and x[supp[0]]==-x[supp[1]]:
                u,v=supp
                twin = all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v)) and A[u][v]==A[v][u]
                kind = ('twin_bidirected' if A[u][v]==1 else 'twin_nonadjacent') if twin else 'pair_nontwin'
            else:
                kind='support%d'%len(supp)
            yield kind, B.tolist(), x.tolist(), y.tolist()
def fpdata(A,eps,delta):
    fp={};st=set()
    for m in range(1,n+1):
        for s in combinations(range(n),m):
            r=ctln.analyze_support(A,list(s),eps,delta)
            if r['stable'] is None or r.get('BOUNDARY'):
                fp[s]='UNDECIDED'
            elif r['is_fixed_point']:
                fp[s]=(tuple(r['x']),tuple(r['charpoly']))
                if r['stable_fixed_point']: st.add(s)
    return fp,st
tested=collections.Counter(); fp_differ=collections.Counter(); st_differ=collections.Counter(); undec=collections.Counter()
examples={}
for idx,kv in enumerate(keys_n):
    A=graphs.from_key(enum_cache.key_from_int(kv,n),n)
    for kind,B,x,y in moves(A):
        if tested[kind]>=CAP: continue
        tested[kind]+=1
        for eps,delta in params:
            f1,s1=fpdata(A,eps,delta); f2,s2=fpdata(B,eps,delta)
            if 'UNDECIDED' in f1.values() or 'UNDECIDED' in f2.values(): undec[kind]+=1
            if set(f1)!=set(f2):
                fp_differ[kind]+=1
                if kind not in examples: examples[kind]=(idx,A,B,x,y,str(eps),str(delta),sorted(set(f1)^set(f2)),sorted(s1),sorted(s2))
                break
            if s1!=s2: st_differ[kind]+=1
            if any(f1[s]!=f2[s] for s in f1): fp_differ[kind+'_values']+=1
    if all(tested[k]>=CAP for k in ('twin_bidirected','twin_nonadjacent','pair_nontwin','support4')): break
print('n',n,'tested',dict(tested)); print('FP-set differs',dict(fp_differ)); print('stable set differs',dict(st_differ)); print('undecided encountered',dict(undec))
for k,e in examples.items():
    idx,A,B,x,y,eps,delta,sd,s1,s2=e
    print(f'EXAMPLE {k}: graph idx {idx} (n={n}) A={show(A)}  B={show(B)}  x={x} y={y}  at eps={eps}, delta={delta}: FP symmetric difference {sd}; stable A {s1} B {s2}')
