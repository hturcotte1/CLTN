"""Exact sampling around q_6* with ctln.analyze_support (decision path exact)."""
import sys, json
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
import sympy as sp
import ctln
from graphs import from_key
from enum_cache import key_from_int
x = sp.Symbol('x')
G = {1000255: 1002304259, 1387528: 1066791683, 1424629: 1069347459}
As = {i: from_key(key_from_int(k, 6), 6) for i, k in G.items()}
print({i: [''.join(map(str, r)) for r in A] for i, A in As.items()})
bpoly = sp.Poly([797, -6468, 22544, -44440, 53280, -38016, 13824], x)
iv = bpoly.intervals(eps=sp.Rational(1, 10**40), inf=sp.Rational(1029, 475), sup=sp.Rational(13, 6))
assert len(iv) == 1
bl, bh = F(str(iv[0][0][0])), F(str(iv[0][0][1]))
def pt(beta, t):
    eps = 1 / (t + 1); delta = (beta - 1) * eps
    return eps, delta, (beta - 1) * t / (t + 1)
def run(i, beta, t):
    eps, delta, q = pt(beta, t)
    assert t * (beta - 2) > 1
    r = ctln.analyze_support(As[i], list(range(6)), eps, delta)
    return q, r['positive'], r['off_support_ok'], r['stable'], r.get('stability_detail')
qstar_lo = bl - 1
qstar_hi = bh - 1
# BELOW: beta = beta_c - h, t large -> q within ~h + 1/t below q*
for h, t in [(F(1, 10**4), F(10**5)), (F(1, 10**6), F(10**7)), (F(1, 10**8), F(10**9))]:
    for i in (1000255, 1424629):
        q, *rest = run(i, bl - h, t)
        print('BELOW', i, 'q*-q <', float(qstar_hi - q), rest)
# ABOVE: beta = beta_c + h  (any point with q > q* has beta > q*+1 = beta_c); scan t
for h in [F(1, 10**8), F(1, 10**5), F(1, 10**3), F(1, 50)]:
    beta = bh + h
    tmin = 1 / (beta - 2)
    for i in (1000255, 1424629, 1387528):
        verdicts = []
        for t in [tmin * F(101, 100), tmin * 2, F(30), F(60), F(100), F(10**3), F(10**5), F(10**9)]:
            if t * (beta - 2) <= 1: continue
            q, pos, off, st, det = run(i, beta, t)
            verdicts.append((pos, st))
        print('ABOVE h=', h, i, 'any stable+positive:', any(p and s is True for p, s in verdicts), 'any undecided:', any(s is None for p, s in verdicts))
