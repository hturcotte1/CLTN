# (1) n=6 sanity check of Theorem J.4: bidirected in-twin moves on the first 120 strict survivors that have one,
#     all 63 supports, two legal parameters; FP set, stable set, x vectors and charpolys must agree exactly.
# (2) G1: D(G1^T) != D(G1); the move G1 -> G2p as a similarity N' = S^{-1} N S with S 1 = 1.
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import pickle
from itertools import combinations
n=6
data=pickle.load(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/s6_out.pkl','rb'))['data']
def rec(A,eps,delta):
    out={}
    for m in range(1,n+1):
        for s in combinations(range(n),m):
            r=ctln.analyze_support(A,list(s),eps,delta)
            out[s]=(r['is_fixed_point'],r['stable_fixed_point'],tuple(r['x']) if r['is_fixed_point'] else None,tuple(r['charpoly']) if r['is_fixed_point'] else None, r['stable'] is None)
    return out
tested=0; bad=0; und=0
for g in sorted(data):
    A=graph_of_idx(g)
    mv=[(u,v,w) for u in range(n) for v in range(u+1,n) if A[u][v]==A[v][u]==1 and all(A[k][u]==A[k][v] for k in range(n) if k not in (u,v)) for w in range(n) if w not in (u,v) and A[u][w]!=A[v][w]]
    if not mv: continue
    u,v,w=mv[0]
    B=[r[:] for r in A]; B[u][w],B[v][w]=A[v][w],A[u][w]
    for eps,delta in ((F(1,4),F(1,2)),(F(1,10),F(1,8))):
        r1,r2=rec(A,eps,delta),rec(B,eps,delta)
        und+=sum(1 for s in r1 if r1[s][4])+sum(1 for s in r2 if r2[s][4])
        if r1!=r2:
            # supports with exactly one twin are allowed to differ only if non-FP in both
            for s in r1:
                if r1[s]!=r2[s] and not (len(set(s)&{u,v})==1 and not r1[s][0] and not r2[s][0]): bad+=1
    tested+=1
    if tested>=120: break
print('bidirected moves tested (n=6 strict survivors):',tested,' violations of J.4:',bad,' undecided supports seen:',und)
# (2)
G1T=[[G1[j][i] for j in range(n)] for i in range(n)]
print('D(G1^T)==D(G1):',Dsig(G1T)==Dsig(G1),' in-degree multisets',sorted(indeg(G1)),sorted(indeg(G1T)))
N=sp.Matrix(Nmat(G1)); G2p=[list(r) for r in graphs.permute(G2,(4,3,5,2,1,0))]; N2=sp.Matrix(Nmat(G2p))
x=sp.Matrix([0,1,-1,0,0,0]); y=N2-N; y=y[:,1]   # N2 - N = e_4 x^T
assert N2-N==y*x.T
zs=sp.symbols('z0:6'); z=sp.Matrix(zs)
sol=sp.solve(list(N*z-y)+[ (x.T*z)[0] ], zs, dict=True)
print('solutions z of N z = y, x^T z = 0:', sol)
if sol:
    zz=z.subs(sol[0]).subs({s:0 for s in zs})
    S=sp.eye(6)+zz*x.T
    print('S =',S.tolist()); print('S*1 == 1:', S*sp.ones(6,1)==sp.ones(6,1), '  S^{-1} N S == N2:', S.inv()*N*S==N2)
