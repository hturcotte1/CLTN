"""W1 and X1 on larger bases (m = 7, 8), exact integer polynomials. Uses determinant lemma for u^T adj v."""
import sys, random
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/referee_theory_1')
from algebra_lib import berkowitz, polysub
random.seed(7)
def padd(a, b):
    L = max(len(a), len(b)); a = [0]*(L-len(a))+a; b=[0]*(L-len(b))+b; return [x+y for x,y in zip(a,b)]
def pmul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def strip(a):
    while a and a[0]==0: a=a[1:]
    return a
def vadju(X, u, v):  # v^T adj(zI - X) u
    n=len(X); Y=[[X[i][j]-u[i]*v[j] for j in range(n)] for i in range(n)]
    return polysub(berkowitz(Y), berkowitz(X))
def nonedge(A): n=len(A); return [[0 if i==j else 1-A[i][j] for j in range(n)] for i in range(n)]
bad=0; cnt=0
for trial in range(12):
    m=random.randint(7,8); H=[[0 if i==j else random.randint(0,1) for j in range(m)] for i in range(m)]
    mu=[random.randint(1,3) for _ in range(m)]
    idx=[a for a in range(m) for _ in range(mu[a])]; n=len(idx)
    B=[[0 if i==j else (1 if idx[i]==idx[j] else H[idx[i]][idx[j]]) for j in range(n)] for i in range(n)]
    lhs=vadju(nonedge(B),[1]*n,[1]*n)
    D=nonedge(H); DM=[[D[i][j]*mu[j] for j in range(m)] for i in range(m)]
    rhs=pmul([1]+[0]*(n-m), vadju(DM,[1]*m,mu))
    cnt+=1; bad+= strip(lhs)!=strip(rhs)
    # X1
    S=[random.randint(0,1) for _ in range(m)]; T=[random.randint(0,1) for _ in range(m)]; k=random.randint(1,6)
    N=m+k; G=[[0]*N for _ in range(N)]
    for i in range(m):
        for j in range(m): G[i][j]=H[i][j]
        for v in range(m,N): G[i][v]=S[i]; G[v][i]=T[i]
    for a in range(m,N):
        for b2 in range(m,N):
            if a!=b2: G[a][b2]=1
    lhs=vadju(nonedge(G),[1]*N,[1]*N)
    L=[[D[i][j]-(1-S[i])-(1-T[j]) for j in range(m)] for i in range(m)]
    zpsi=pmul([1,0], vadju(D,[1]*m,[1]*m)); kchi=[k*c for c in berkowitz(L)]
    rhs=pmul([1]+[0]*(k-1), padd(zpsi,kchi))
    cnt+=1; bad+= strip(lhs)!=strip(rhs)
print('checks',cnt,'mismatches',bad)
