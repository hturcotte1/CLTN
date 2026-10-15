from fractions import Fraction as F
import sys
G1=[[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
def build(k):
    n=6+k; A=[[0]*n for _ in range(n)]
    for i in range(6):
        for j in range(6): A[i][j]=G1[i][j]
    for a in range(6,n):
        for g in range(6): A[a][g]=1
        for b in range(6,n):
            if a!=b: A[a][b]=1
    for g in range(5):
        for a in range(6,n): A[g][a]=1
    return A
def charpoly(M):  # Faddeev-LeVerrier, exact; det(zI-M)
    n=len(M); c=[F(1)]; Mk=[[F(0)]*n for _ in range(n)]
    I=[[F(int(i==j)) for j in range(n)] for i in range(n)]
    for k in range(1,n+1):
        # Mk = M*Mk_prev + c_{k-1} I
        P=[[sum(M[i][l]*Mk[l][j] for l in range(n) if Mk[l][j]!=0) for j in range(n)] for i in range(n)]
        Mk=[[P[i][j]+c[-1]*I[i][j] for j in range(n)] for i in range(n)]
        AM=[[sum(M[i][l]*Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c.append(-sum(AM[i][i] for i in range(n))/k)
    return c
def routh(c):  # strict Hurwitz via Routh first column; returns True/False/None(undecided)
    r0=c[0::2]; r1=c[1::2]
    if any(x<=0 for x in c): return False
    rows=[r0,r1]
    n=len(c)-1
    for _ in range(n-1):
        a,b=rows[-2],rows[-1]
        if b[0]==0: return None
        new=[]
        for i in range(len(a)-1):
            bi=b[i+1] if i+1<len(b) else F(0)
            new.append((b[0]*a[i+1]-a[0]*bi)/b[0])
        if not new: new=[F(0)]
        rows.append(new)
    col=[r[0] for r in rows if r]
    if any(x==0 for x in col): return None
    return all(x>0 for x in col)
def solve(K):
    n=len(K); M=[row[:]+[F(1)] for row in K]
    for c in range(n):
        p=next(r for r in range(c,n) if M[r][c]!=0); M[c],M[p]=M[p],M[c]
        pv=M[c][c]; M[c]=[x/pv for x in M[c]]
        for r in range(n):
            if r!=c and M[r][c]!=0:
                f=M[r][c]; M[r]=[x-f*y for x,y in zip(M[r],M[c])]
    return [M[i][n] for i in range(n)]
def test(k,beta,t):
    A=build(k); n=len(A)
    K=[[F(int(i==j))+(beta if (i!=j and A[j][i]==0) else 0)+t for j in range(n)] for i in range(n)]
    x=solve(K); pos=all(v>0 for v in x)
    cp=charpoly([[-v for v in r] for r in K])  # det(zI+K)=det(zI-(-K))
    return pos, routh(cp)
if __name__=="__main__":
    beta=F(sys.argv[1]); t=F(sys.argv[2]); ks=[int(a) for a in sys.argv[3:]]
    eps=1/(1+t); delta=(beta-1)*eps; print("eps",eps,"delta",delta,"q",delta*(1-eps)/eps)
    for k in ks: print(k,6+k,test(k,beta,t),flush=True)
