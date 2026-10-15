import sys
sys.path.insert(0,'/home/user/CLTN')
from fractions import Fraction as F
import ctln
H=[[0,0,1,1,1,1,0],[0,0,1,0,1,0,0],[1,1,0,1,1,1,1],[1,1,1,0,1,1,1],[0,0,1,1,0,1,1],[0,1,0,0,1,0,1],[1,1,0,1,0,0,0]]
def blow(H,mu):
    cls=[a for a,m in enumerate(mu) for _ in range(m)]
    n=len(cls)
    return [[0 if i==j else (1 if cls[i]==cls[j] else H[cls[i]][cls[j]]) for j in range(n)] for i in range(n)]
beta=F(53,26)
for a in [0,1,2,3]:
    A=blow(H,[1,1,2,a,1,1,1]); n=len(A)
    res=[]
    for t in [F(27),F(40),F(53),F(107,2),F(54),F(60)]:
        eps=1/(t+1); delta=(beta-1)*eps
        r=ctln.analyze_support(A,list(range(n)),eps,delta)
        res.append((str(t), r['positive'], r['stable']))
    print('a=',a,'n=',n,res)
print('--- detail at t=107/2')
t=F(107,2); eps=1/(t+1); delta=(beta-1)*eps
for a in [0,1]:
    A=blow(H,[1,1,2,a,1,1,1]); n=len(A)
    r=ctln.analyze_support(A,list(range(n)),eps,delta)
    print('a=',a,'eps',eps,'delta',delta,{k:r[k] for k in r if k in ('positive','stable','stability_detail','BOUNDARY','legal','stable_fixed_point','off_support_ok')})
