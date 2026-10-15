"""Referee helper library: own digraph6 decoder, own family builders. Exact only."""
import sys
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as Fr

def d6_decode(s):
    assert s[0] == '&'
    body = [ord(c) - 63 for c in s[1:]]
    n = body[0]; assert n < 63
    bits = []
    for v in body[1:]:
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    A = [[bits[i*n+j] for j in range(n)] for i in range(n)]
    for i in range(n): assert A[i][i] == 0
    return A

R_D6 = '&F^U|m^QTq?'
S_R = [1,1,0,1,1,1,1]

def build_Rk(k, A0=None, S=S_R):
    """core nodes 0..m-1; class nodes m..m+k-1 bidirected clique, class -> all core,
    core i -> class iff S[i]."""
    A0 = A0 or d6_decode(R_D6)
    m = len(A0); n = m + k
    A = [[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): A[i][j] = A0[i][j]
    for c in range(m, n):
        for j in range(n):
            if j != c and j >= m: A[c][j] = 1
        for i in range(m):
            A[c][i] = 1
            if S[i]: A[i][c] = 1
    return A

def selfjoin(A0, k):
    m = len(A0); n = m*k
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if i // m == j // m: A[i][j] = A0[i % m][j % m]
            else: A[i][j] = 1
    return A
