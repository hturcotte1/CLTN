"""EXPLORATION, exact: for the top non-monotone examples, is (3) non-monotone in t along a LEGAL vertical line?
And does (1) hold there (would be a counterexample -> must not happen at n=5, which is certified)."""
import sys
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import ctln
exs = [[[0,1,0,0,0],[1,0,1,0,0],[1,0,0,1,1],[1,1,1,0,1],[1,1,1,1,0]],
       [[0,0,0],[1,0,1],[1,1,0]],
       [[0,1,0,0],[1,0,1,0],[1,0,0,1],[1,1,1,0]]]
for A in exs:
    n = len(A)
    for beta in [F(10), F(5), F(3)]:
        row = []
        for t in [F(1)/(beta-2) + F(1, 1000), F(3, 20), F(1, 4), F(1, 2), F(1), F(5), F(100), F(10**6)]:
            if t*(beta-2) <= 1: continue
            eps = 1/(1+t); delta = (beta-1)/(1+t)
            assert ctln.legal(eps, delta)
            r = ctln.analyze_support(A, list(range(n)), eps, delta)
            row.append(f"t={float(t):.4g}:{'S' if r['stable'] else ('U' if r['stable'] is False else '?')}{'+' if r['positive'] else '-'}")
        print(n, f"beta={beta}", " ".join(row))
print("legend: S/U/? = (3) stable/unstable/undecided; +/- = (1) positive / not")
