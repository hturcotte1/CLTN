# n=6 STRICT survivor pairs in one D-class (non-adjacent in-twin move) with different FP(G) / stable sets. Exact.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import pickle
from itertools import combinations
n=6
data=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','rb'))['data']
K=keys(); pos={k:i for i,k in enumerate(K)}
def fps(A,eps,delta):
    out={}
    for m in range(1,n+1):
        for s in combinations(range(n),m):
            r=ctln.analyze_support(A,list(s),eps,delta)
            assert r['stable'] is not None and not r.get('BOUNDARY')
            if r['is_fixed_point']: out[s]=r['stable_fixed_point']
    return out
found=0
for g in sorted(data):
    A=graph_of_idx(g)
    for u in range(n):
        for v in range(u+1,n):
            if not (A[u][v]==A[v][u]==0 and all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v))): continue
            for w in range(n):
                if w in (u,v) or A[u][w]==A[v][w]: continue
                B=[r[:] for r in A]; B[u][w],B[v][w]=A[v][w],A[u][w]
                h=pos[enum_cache.key_to_int(graphs.canonical_form(B)[0])]
                assert data[h]==data[g]
                for eps,delta in ((F(1,4),F(1,2)),(F(1,10),F(1,8))):
                    f1,f2=fps(A,eps,delta),fps(B,eps,delta)
                    s1=sorted(s for s in f1 if f1[s]); s2=sorted(s for s in f2 if f2[s])
                    if h!=g and (sorted(map(len,s1))!=sorted(map(len,s2)) or len(f1)!=len(f2)):
                        print(f'STRICT survivors idx {g}, {h}: same D, non-adjacent in-twins u={u}, v={v}, move at w={w}; eps={eps}, delta={delta}, q={ctln.q_of(eps,delta)}')
                        print('   A ',show(A),' stable supports',s1,' |FP|',len(f1),' FP sizes',sorted(map(len,f1)))
                        print('   B ',show(B),' stable supports',s2,' |FP|',len(f2),' FP sizes',sorted(map(len,f2)))
                        found+=1; break
                if found>=4: raise SystemExit
