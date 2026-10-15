"""Referee misc: (a) F_1 isomorphic to the n=7 candidate &FNFL|\\u|o? and present in the repo candidate list;
(b) Lemma I.4 spot check: sign of x(t) for the full support equals sign pattern of w, for several t (exact)."""
import sys, json, gzip, subprocess
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/n7')
from fractions import Fraction as Fr
import ctln
G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
def Fk(kk):
    n = 6 + kk; A = [[0]*n for _ in range(n)]
    for i in range(6):
        for j in range(6): A[i][j] = G1[i][j]
    for v in range(6, n):
        for i in range(6):
            A[v][i] = 1; A[i][v] = 1 if i <= 4 else 0
        for u in range(6, n):
            if u != v: A[v][u] = 1
    return A
def to_d6(A):
    n = len(A); bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(n + 63)
    for i in range(0, len(bits), 6):
        v = 0
        for bb in bits[i:i+6]: v = 2*v + bb
        s += chr(v + 63)
    return s
def canon(s):
    out = subprocess.run(['/home/user/CLTN/n7/labelg', '-q'], input=s + '\n', capture_output=True, text=True).stdout.strip()
    return out
c1 = canon(to_d6(Fk(1))); c2 = canon('&FNFL|\\u|o?')
print('canon F_1 =', c1, ' canon &FNFL|\\u|o? =', c2, ' same:', c1 == c2)
cands = [json.loads(l) for l in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
keys = set()
for c in cands:
    keys.add(canon(c['d6']))
print('candidate records:', len(cands), 'F_1 in candidate list:', c1 in keys, ' (fields:', list(cands[0].keys())[:8], ')')
# (b) Lemma I.4 spot check
for kk in [0, 2, 4]:
    A = Fk(kk); n = len(A)
    for beta in [Fr(21, 10), Fr(23, 10), Fr(2325, 1000), Fr(5, 2)]:
        g = 1 - 1/beta
        w = ctln.solve_exact([[Fr(A[j][i]) + (g if i == j else 0) for j in range(n)] for i in range(n)], [Fr(1)]*n)
        wpos = all(v > 0 for v in w)
        for t in [Fr(1, 10), Fr(1), Fr(7), Fr(50), Fr(10**4)]:
            eps = 1/(1+t); delta = (beta-1)/(1+t)
            rec = ctln.analyze_support(A, list(range(n)), eps, delta)
            assert rec['positive'] == wpos, (kk, beta, t)
print('Lemma I.4 spot check (positivity t-free and equal to w > 0): OK on 60 points')
