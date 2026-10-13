# For D-classes still split after all rank-one moves: take one graph from each of two components,
# find all relabellings matching in-walk vectors (A^T)^k 1, and report the minimal rank of N_B' - N_A
# and the structure of the Krylov space.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import pickle, collections
n=6
S6=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','rb')); data=S6['data']
S8=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s8_out.pkl','rb'))
root={g:r for r,c in S8['comp'].items() for g in c}
dcls=collections.defaultdict(list)
for g in data: dcls[data[g]].append(g)
split=[c for c in dcls.values() if len({root[g] for g in c})>1]
split.sort(key=len)
def walkvecs(A):
    At=sp.Matrix(A).T; v=sp.ones(n,1); out=[v]
    for k in range(n): v=At*v; out.append(v)
    return out
rankhist=collections.Counter(); krydim=collections.Counter()
for ci,c in enumerate(split):
    comps=collections.defaultdict(list)
    for g in c: comps[root[g]].append(g)
    reps=[v[0] for v in comps.values()]
    A=graph_of_idx(reps[0]); B=graph_of_idx(reps[1])
    WA=sp.Matrix.hstack(*walkvecs(A)); NA=sp.Matrix(Nmat(A))
    best=None
    for p in itertools.permutations(range(n)):
        Bp=[list(r) for r in graphs.permute(B,p)]
        if sp.Matrix.hstack(*walkvecs(Bp))!=WA: continue
        r=(sp.Matrix(Nmat(Bp))-NA).rank()
        if best is None or r<best[0]: best=(r,p,Bp)
    assert best is not None, 'no walk-matching relabelling -- contradicts Theorem J.2'
    kd=WA.rank(); rankhist[best[0]]+=1; krydim[kd]+=1
    if ci<6:
        print('class size',len(c),'components',len(comps),'Krylov dim',kd,'min rank of difference',best[0])
        print('   A ',show(A)); print('   Bp',show(best[2])); print('   N_B-N_A', (sp.Matrix(Nmat(best[2]))-NA).tolist())
print('min-rank histogram over the',len(split),'split classes:',dict(rankhist)); print('Krylov dim histogram',dict(krydim))
