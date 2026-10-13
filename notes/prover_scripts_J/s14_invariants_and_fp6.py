# (1) Within every D-class of the 34,071 survivors: which graph invariants are constant? (exact integers)
# (2) n=6 survivor pairs in one D-class with different FP(G)/stable sets at a legal parameter (exact).
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import collections
from itertools import combinations
n=6
cls=classes()
def inv(A):
    cp=tuple(fastregion.icharpoly([r[:] for r in A]))
    cpJ=tuple(fastregion.icharpoly([[A[i][j]+1 for j in range(n)] for i in range(n)]))   # chi_(A+J)
    perp=tuple(x-y for x,y in zip(cp,cpJ))                                                 # = n chi_(A|perp) (shifted index)
    two=sum(A[i][j]*A[j][i] for i in range(n) for j in range(i+1,n))
    edges=sum(map(sum,A))
    import numpy as np
    M=np.array(A,dtype=object); P=np.identity(n,dtype=object); tw=[]
    for k in range(n+1): tw.append(int(np.ones(n,dtype=object)@P@np.ones(n,dtype=object))); P=P@M
    return {'indeg':tuple(sorted(indeg(A))),'outdeg':tuple(sorted(outdeg(A))),'edges':edges,'two_cycles':two,'chiA':cp,'chiAperp':perp,'total_walks':tuple(tw),
            'chiAT_eq_chiA':True,'domfree':graphs.is_domination_free(A)}
const=collections.Counter(); nontriv=0
for c in cls:
    if len(c)<2: continue
    nontriv+=1
    I=[inv(graph_of_idx(g)) for g in c]
    for k in I[0]:
        if all(x[k]==I[0][k] for x in I): const[k]+=1
print('non-singleton D-classes:',nontriv); print('number of those on which the invariant is constant:',dict(const))
# (2) FP(G) differences inside survivor D-classes, via non-adjacent in-twin moves (strict survivors)
K=keys(); pos={k:i for i,k in enumerate(K)}
surv=set(int(l) for l in open('/home/user/CLTN/certificates/n6_survivors.idx') if l.strip())
eps,delta=F(1,4),F(1,2)
def fps(A):
    out={}
    for m in range(1,n+1):
        for s in combinations(range(n),m):
            r=ctln.analyze_support(A,list(s),eps,delta)
            assert r['stable'] is not None and not r.get('BOUNDARY')
            if r['is_fixed_point']: out[s]=r['stable_fixed_point']
    return out
found=0
for g in sorted(surv):
    A=graph_of_idx(g)
    for u in range(n):
        for v in range(u+1,n):
            if not (A[u][v]==A[v][u]==0 and all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v))): continue
            for w in range(n):
                if w in (u,v) or A[u][w]==A[v][w]: continue
                B=[r[:] for r in A]; B[u][w],B[v][w]=A[v][w],A[u][w]
                h=pos[enum_cache.key_to_int(graphs.canonical_form(B)[0])]
                if h not in surv: continue
                f1,f2=fps(A),fps(B)
                if set(f1)!=set(f2) or {s for s in f1 if f1[s]}!={s for s in f2 if f2[s]}:
                    print(f'n=6 survivors idx {g} and {h} (same D-class; non-adjacent in-twins u={u}, v={v}, w={w}):')
                    print('   A',show(A),' FP',sorted(f1),' stable',sorted(s for s in f1 if f1[s]))
                    print('   B',show(B),' FP',sorted(f2),' stable',sorted(s for s in f2 if f2[s]))
                    found+=1
                if found>=2: raise SystemExit
