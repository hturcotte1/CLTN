# Theorem J.2 check: partition of the survivors by D(A) == partition by
#   I(A) := (charpoly of A, multiset over vertices i of the sequence (((A^T)^k 1)_i)_{k=0..n}).
# Exact integers throughout (sympy charpoly over ZZ; integer matrix powers).
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import collections
n=6
surv=[int(l) for l in open('/home/user/CLTN/certificates/n6_survivors.idx') if l.strip()]
cls=classes(); cls_of={g:ci for ci,c in enumerate(cls) for g in c}
def inv(A):
    At=[[A[j][i] for j in range(n)] for i in range(n)]
    v=[1]*n; seqs=[[1] for _ in range(n)]
    for k in range(n):
        v=[sum(At[i][j]*v[j] for j in range(n)) for i in range(n)]
        for i in range(n): seqs[i].append(v[i])
    cp=tuple(fastregion.icharpoly([row[:] for row in A]))
    return (cp, tuple(sorted(tuple(s) for s in seqs)))
I={g:inv(graph_of_idx(g)) for g in surv}
icls=collections.defaultdict(list)
for g in surv: icls[I[g]].append(g)
P1=sorted(sorted(c) for c in cls); P2=sorted(sorted(c) for c in icls.values())
print('D-classes',len(P1),'I-classes',len(P2),'partitions identical:',P1==P2)
