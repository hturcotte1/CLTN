exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
n=6
K=keys(); pos={k:i for i,k in enumerate(K)}
surv=[int(l) for l in open('/home/user/CLTN/certificates/n6_survivors.idx') if l.strip()]; survset=set(surv)
def in_twins(A):
    return [(u,v) for u in range(n) for v in range(u+1,n) if A[u][v]==A[v][u] and all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v))]
shown=0
for g in surv:
    A=graph_of_idx(g)
    for (u,v) in in_twins(A):
        for w in range(n):
            if w in (u,v) or A[u][w]==A[v][w]: continue
            B=[r[:] for r in A]; B[u][w],B[v][w]=A[v][w],A[u][w]
            h=pos[enum_cache.key_to_int(graphs.canonical_form(B)[0])]
            if h in survset: continue
            ba=fastregion.beta_data(A,list(range(n))); bb=fastregion.beta_data(B,list(range(n)))
            print('g',g,'h',h,(u,v,w),'a equal',ba['a']==bb['a'],'d eq',ba['d']==bb['d'])
            ra=fastregion.region(A,list(range(n)),None,bd=ba,positivity_only=True); rb=fastregion.region(B,list(range(n)),None,bd=bb,positivity_only=True)
            print('  verdicts', ra['verdict'], rb['verdict'], 'A',show(A),'B',show(B))
            Bc=graphs.from_key(enum_cache.key_from_int(K[h],6),6)
            rc=fastregion.region(Bc,list(range(n)),None,positivity_only=True); print('  canonical B verdict', rc['verdict'])
            shown+=1
            if shown>=4: raise SystemExit
