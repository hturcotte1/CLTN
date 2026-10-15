import sys, time
sys.path.insert(0,'/home/user/CLTN')
from fractions import Fraction as F
import ctln
n=int(sys.argv[1]); k=n-6
G1=[[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
S=[1,1,1,1,1,0]
B=[[0]*n for _ in range(n)]
for i in range(6):
    for j in range(6): B[i][j]=G1[i][j]
    for v in range(6,n): B[i][v]=S[i]; B[v][i]=1
for a in range(6,n):
    for b in range(6,n):
        if a!=b: B[a][b]=1
t=F(int(sys.argv[2])); beta=F(sys.argv[3])
eps=1/(1+t); delta=(beta-1)/(1+t)
t0=time.time(); r=ctln.analyze_support(B,list(range(n)),eps,delta); print(time.time()-t0, r['stable_fixed_point'], r['stable'], r['stability_detail'], r['positive'])
