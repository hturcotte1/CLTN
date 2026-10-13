# Components of the strict survivors under BIDIRECTED in-twin switching only (the FP(G)-preserving move, J.4),
# and the structure of the largest D-class.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import pickle, collections
n=6
S6=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','rb')); data=S6['data']; strict=sorted(data)
K=keys(); pos={k:i for i,k in enumerate(K)}
parent={g:g for g in strict}
def find(a):
    while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
    return a
def twins(A, bidir):
    return [(u,v) for u in range(n) for v in range(u+1,n) if A[u][v]==A[v][u]==(1 if bidir else 0) and all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v))]
has_bt=0; has_nt=0
for g in strict:
    A=graph_of_idx(g)
    bt=twins(A,True); nt=twins(A,False); has_bt+=bool(bt); has_nt+=bool(nt)
    for (u,v) in bt:
        for w in range(n):
            if w in (u,v) or A[u][w]==A[v][w]: continue
            B=[r[:] for r in A]; B[u][w],B[v][w]=A[v][w],A[u][w]
            h=pos[enum_cache.key_to_int(graphs.canonical_form(B)[0])]
            assert data[h]==data[g]
            a,b=find(g),find(h)
            if a!=b: parent[a]=b
comp=collections.defaultdict(list)
for g in strict: comp[find(g)].append(g)
print('strict survivors with a bidirected in-twin pair:',has_bt,' with a non-adjacent in-twin pair:',has_nt)
print('components under bidirected in-twin switching:',len(comp),' (D-classes: 17307)')
print('nontrivial components:',sum(1 for c in comp.values() if len(c)>1),' graphs in them:',sum(len(c) for c in comp.values() if len(c)>1), ' max size', max(len(c) for c in comp.values()))
# largest D-class
dcls=collections.defaultdict(list)
for g in strict: dcls[data[g]].append(g)
big=max(dcls.values(),key=len)
print('largest D-class size',len(big),' in-degree multiset', sorted(indeg(graph_of_idx(big[0]))))
outs=collections.Counter(tuple(sorted(outdeg(graph_of_idx(g)))) for g in big); print('  out-degree sequences in it:',len(outs), outs.most_common(5))
A=graph_of_idx(big[0]); print('  example member', show(A), ' bidirected twins', twins(A,True), ' nonadjacent twins', twins(A,False))
# in-twin classes (equivalence classes of the in-twin relation) in that example
