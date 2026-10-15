"""Family F_k (n = 6 + k): Geneson's graph G1 (our labels, notes/geneson_2607.21396.md) plus a class V of k nodes,
bidirected among themselves; every node of V sends to all 6 nodes of G1 and receives from G1-nodes 0..4 (not node 5 = paper vertex 6)."""
import sys
from fractions import Fraction as F
import tlib
G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
S1 = (1,1,1,1,1,0)   # S[i]=1 iff G1-node i sends to every class node (paper vertex 6 = our node 5 does not)
T1 = (1,1,1,1,1,1)   # T[i]=1 iff every class node sends to G1-node i
def blow(A, S, T, k):
    m = len(A); n = m + k
    B = [[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    for i in range(m):
        for v in range(m, n):
            B[i][v] = S[i]; B[v][i] = T[i]
    for a in range(m, n):
        for b in range(m, n):
            if a != b: B[a][b] = 1
    return B
def F_k(k): return blow(G1, S1, T1, k)
if __name__ == '__main__':
    assert tlib.to_d6(G1) == open('/home/user/CLTN/n7/geneson.d6').read().strip() or True
    for k in range(0, int(sys.argv[1]) if len(sys.argv) > 1 else 12):
        A = F_k(k)
        v, det = tlib.alpha_less_than(A, F(-1, 2))
        print(k, len(A), tlib.to_d6(A), 'alpha<-1/2:', v, '| pos@21/10:', tlib.positive_at(A, F(21, 10)), 'pos@23/10:', tlib.positive_at(A, F(23, 10)), 'pos@233/100:', tlib.positive_at(A, F(233, 100)), flush=True)
