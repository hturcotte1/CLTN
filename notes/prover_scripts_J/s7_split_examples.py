exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import pickle, collections
S=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','rb'))
data=S['data']; comp=S['comp']
root={}
for r,c in comp.items():
    for g in c: root[g]=r
dcls=collections.defaultdict(list)
for g in data: dcls[data[g]].append(g)
small=[c for c in dcls.values() if len({root[g] for g in c})>1]
small.sort(key=len)
print('sizes of split classes (smallest 10):',[len(c) for c in small[:10]])
n=6
def best_relabel(A,B):
    best=None
    for p in itertools.permutations(range(n)):
        Bp=graphs.permute(B,p)
        diff=[(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]!=Bp[i][j]]
        if best is None or len(diff)<len(best[0]): best=(diff,p,Bp)
    return best
for c in small[:6]:
    reps={}
    for g in c: reps.setdefault(root[g],g)
    gs=list(reps.values())
    A=graph_of_idx(gs[0]); B=graph_of_idx(gs[1])
    diff,p,Bp=best_relabel(A,B)
    print('class',c,'comps',len(reps))
    print('  A ',show(A)); print('  Bp',show(Bp),' diff',diff)
    NA=sp.Matrix(Nmat(A)); NB=sp.Matrix(Nmat([list(r) for r in Bp]))
    print('  N_B - N_A rank', (NB-NA).rank(), (NB-NA).tolist())
    # left null/eigen vectors of N_A orthogonal to 1
    for lam in set(NA.eigenvals().keys()):
        if not lam.is_rational: continue
        ns=(NA.T-lam*sp.eye(n)).nullspace()
        print('   eigenvalue',lam,'left eigvecs',[list(v) for v in ns])
