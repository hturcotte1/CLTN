# Exact check of Remark J.4': for random supports containing in-twins relative to sigma, -I+W_sigma has the
# exact eigenvalue -eps (bidirected) or +delta (non-adjacent): xi^T(-I+W_sigma) = lam xi^T. Exact Fractions.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
from itertools import combinations
n=6; K=keys(); cnt={'bi':0,'non':0}; bad=0; unstable_ok=0
for idx in range(0,len(K),997):
    A=graph_of_idx(idx)
    for eps,delta in ((F(1,4),F(1,2)),(F(1,10),F(1,8))):
        W=ctln.weight_matrix(A,eps,delta)
        for m in range(2,n+1):
            for s in combinations(range(n),m):
                for u,v in combinations(s,2):
                    if A[u][v]!=A[v][u] or any(A[k][u]!=A[k][v] for k in s if k not in (u,v)): continue
                    lam = -eps if A[u][v]==1 else delta
                    M=[[(F(-1) if i==j else F(0))+W[s[i]][s[j]] for j in range(m)] for i in range(m)]
                    iu,iv=s.index(u),s.index(v)
                    row=[M[iu][j]-M[iv][j] for j in range(m)]
                    xi=[F(1) if j==iu else F(-1) if j==iv else F(0) for j in range(m)]
                    if row!=[lam*x for x in xi]: bad+=1
                    cnt['bi' if A[u][v] else 'non']+=1
                    if A[u][v]==0 and m<=5:
                        r=ctln.analyze_support(A,list(s),eps,delta)
                        if r['stable'] is False: unstable_ok+=1
                        elif r['stable'] is True: bad+=1
print('checked pairs',cnt,' eigen-identity failures / stable non-adjacent-twin supports:',bad,' non-adjacent supports confirmed unstable by Hurwitz:',unstable_ok)
