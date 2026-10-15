"""
exact_core.py -- generic EXACT primitives (integers / Fractions; no floating point anywhere).
Conventions CLAUDE.md section 2: A[i][j] = 1 iff i -> j.  D := J - I - A (non-edge matrix).  beta = 1 + delta/eps,
r = 1/beta, gamma = 1 - r.
  psi_D(z) = 1^T adj(zI - D) 1 = chi_D(z) - chi_{D+J}(z)   (= m chi_{D|perp}, Lemma I.1; matrix determinant lemma)
  Lemma D1: alpha(A|perp) < -gamma  iff  every root of psi_D has Re z > -r.
  Lemma D2: w = (A^T + gamma I)^{-1} 1 > 0  iff  u = (D^T + rI)^{-1} 1 > 0 and 1^T u > 1.
"""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/notes/research_wf3/theory')
import ctln
from fastregion import icharpoly
from tlib import taylor_shift

def D_of_A(A):
    m = len(A); return [[0 if i == j else 1 - A[i][j] for j in range(m)] for i in range(m)]

def psi(D):
    """hi->lo integer coefficients of psi_D (degree m-1, leading coefficient m)."""
    m = len(D)
    a = icharpoly(D); b = icharpoly([[D[i][j] + 1 for j in range(m)] for i in range(m)])
    out = [x - y for x, y in zip(a, b)][1:]
    assert out[0] == m
    return out

def roots_right_of(p, s):
    """exact three-valued: every root of p (hi->lo, p[0] > 0) has Re z > -s.  Returns (True|False|None, detail).
    Q(w) := (-1)^deg p(-w - s) has roots w = -z - s; Q Hurwitz  <=>  all Re z > -s.  A vanishing minor -> None."""
    p = [F(x) for x in p]
    while p and p[0] == 0: p = p[1:]
    d = len(p) - 1
    ptil = [c * (-1) ** j for j, c in enumerate(p)]          # (-1)^d p(-y), hi->lo in y
    Q = taylor_shift(ptil, F(s))                               # (-1)^d p(-(w + s))
    if Q[0] < 0: Q = [-c for c in Q]
    return ctln.hurwitz_stable(Q)

def strip_zero_roots(p):
    p = list(p)
    k = 0
    while p and p[-1] == 0: p.pop(); k += 1
    return p, k

def limit_stable(A, beta):
    """exact: alpha(A|perp) < -(1 - 1/beta), via psi_D (zero roots stripped: they have Re 0 > -r)."""
    p, _ = strip_zero_roots(psi(D_of_A(A)))
    return roots_right_of(p, 1 / F(beta))

def positivity(A, beta):
    """exact u, U and the verdict w > 0 (Lemma D2); also returns w."""
    D = D_of_A(A); m = len(D); r = 1 / F(beta)
    M = [[F(D[j][i]) + (r if i == j else 0) for j in range(m)] for i in range(m)]
    u = ctln.solve_exact(M, [F(1)] * m)
    if u is None: return {'ok': None, 'reason': 'singular'}
    U = sum(u)
    ok = min(u) > 0 and U > 1
    w = [x / (U - 1) for x in u] if U != 1 else None
    return {'ok': ok, 'u': u, 'U': U, 'w': w}

def sink_class_extension(A, S, k):
    """G[S, all; k]: k new nodes, bidirected among themselves, each sends to every old node, receives from i iff S[i]."""
    m = len(A); n = m + k
    B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    for v in range(m, n):
        for i in range(m):
            B[v][i] = 1
            B[i][v] = 1 if S[i] else 0
        for w in range(m, n):
            if w != v: B[v][w] = 1
    return B

def is_complete(A):
    m = len(A); return all(A[i][j] == 1 for i in range(m) for j in range(m) if i != j)

def ctln_witness(A, beta, t):
    """explicit exact SFP check by ctln.analyze_support at eps = 1/(1+t), delta = (beta-1)/(1+t)."""
    beta = F(beta); t = F(t)
    eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
    n = len(A)
    rec = ctln.analyze_support(A, list(range(n)), eps, delta)
    return eps, delta, rec
