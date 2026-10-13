# Is D(A) already determined by weaker data? Compare the D-partition of the 34,071 survivors with partitions by
#   W1 = (chi_A, chi_(A+J), in-degree multiset)          [generalized spectrum + in-degrees]
#   Wk = (chi_A, multiset of per-vertex in-walk sequences truncated at length k), k = 1..6
# All exact integers.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import collections
n=6
surv=[int(l) for l in open('/home/user/CLTN/certificates/n6_survivors.idx') if l.strip()]
cls=classes(); D={g:ci for ci,c in enumerate(cls) for g in c}
rows={}
for g in surv:
    A=graph_of_idx(g)
    cp=tuple(fastregion.icharpoly([r[:] for r in A])); cpJ=tuple(fastregion.icharpoly([[A[i][j]+1 for j in range(n)] for i in range(n)]))
    At=[[A[j][i] for j in range(n)] for i in range(n)]
    v=[1]*n; seqs=[[1] for _ in range(n)]
    for k in range(n):
        v=[sum(At[i][j]*v[j] for j in range(n)) for i in range(n)]
        for i in range(n): seqs[i].append(v[i])
    rows[g]=(cp,cpJ,seqs)
def npart(f):
    P=collections.defaultdict(set)
    for g in surv: P[f(g)].add(D[g])
    merged=sum(1 for s in P.values() if len(s)>1)
    return len(P), merged
print('D-classes',len(cls))
print('W1 (chi_A, chi_(A+J), in-degrees): classes, #W1-classes merging >1 D-class =', npart(lambda g:(rows[g][0],rows[g][1],tuple(sorted(s[1] for s in rows[g][2])))))
for k in range(1,n+1):
    print(f'W{k} (chi_A, in-walk sequences up to length {k}):', npart(lambda g:(rows[g][0],tuple(sorted(tuple(s[:k+1]) for s in rows[g][2])))))
