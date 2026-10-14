"""
EXACT. Certifies that the PURE-STABILITY form of Lemma M (T(beta) upward closed; and the weak form)
is FALSE, on a 3-node graph. Positivity fails for this graph (node 1 is dominated by node 2), so
this is NOT a CTLN counterexample of any kind; it only shows that Lemma M needs the positivity
hypothesis a(beta) > 0 (or domination-freeness).
"""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln, graphs

A = [[0, 1, 1],
     [0, 0, 1],
     [1, 0, 0]]
sigma = [0, 1, 2]
print("adjacency (A[i][j]=1 iff i->j):")
for r in A: print("   ", r)
beta = F(22, 5)
for t in [F(11, 25), F(1), F(100)]:
    eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
    rec = ctln.analyze_support(A, sigma, eps, delta)
    print(f"beta={beta} t={t}: eps={eps} delta={delta} q={ctln.q_of(eps, delta)} legal={ctln.legal(eps, delta)}"
          f"  stable={rec['stable']} positive={rec['positive']} x={rec['x']}  charpoly={rec['charpoly']}")
print("dominated nodes:", graphs.dominated_nodes(A))
# exact limit spectrum of A|perp: charpoly of compression with rational basis
import sympy as sp
z = sp.symbols('z')
m = 3
At = sp.Matrix(A)
V = sp.Matrix(m, m - 1, lambda i, j: 1 if i == j else (-1 if i == m - 1 else 0))
Gp = (V.T * V).inv() * V.T * At * V
print("charpoly of A|perp:", sp.factor(Gp.charpoly(z).as_expr()), " roots:", sp.solve(Gp.charpoly(z).as_expr(), z))
print("gamma = 1-1/beta =", 1 - sp.Rational(5, 22))
