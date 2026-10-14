import random, sys
sys.path.insert(0,'.')
from check import dominated
def to_d6(A):
    n=len(A); bits=[A[i][j] for i in range(n) for j in range(n)]
    while len(bits)%6: bits.append(0)
    return '&'+chr(n+63)+''.join(chr(63+int(''.join(map(str,bits[k:k+6])),2)) for k in range(0,len(bits),6))
rng=random.Random(int(sys.argv[3]))
n=int(sys.argv[1]); cnt=int(sys.argv[2])
out=[]
while len(out)<cnt:
    d=rng.choice([0.3,0.45,0.6,0.7,0.8,0.85,0.9,0.95])
    A=[[0 if i==j else int(rng.random()<d) for j in range(n)] for i in range(n)]
    want_free = len(out)%4!=0     # 3/4 domination-free (to exercise Descartes), 1/4 unconstrained
    if want_free and dominated(A) is not None: continue
    out.append(to_d6(A))
print('\n'.join(out))
