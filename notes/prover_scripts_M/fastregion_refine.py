"""EXACT: refine the right endpoints of the stable beta-intervals (fastregion, repo engine) to rational isolating
intervals of width < 1e-12 and compare them: finite-t endpoint vs limit endpoint beta_c."""
import sys
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')
import fastregion
from certify_M_violation import from_d6
A = from_d6(sys.argv[1]); sig = list(range(len(A)))
x = sp.symbols('x')
res = {}
TS = [fastregion.INF] + [F(x) for x in (sys.argv[2:] or ['10', '30', '100', '1000', '10000'])]
for t in TS:
    cert = fastregion.region(A, sig, t)
    if 'stable_intervals' not in cert:
        print('t =', t, 'no stable beta-interval'); continue
    c = cert['stable_intervals'][-1]['hi']
    lo, hi = F(c['beta_interval'][0]), F(c['beta_interval'][1])
    P = sp.Poly(c['beta_minpoly_hi_to_lo'], x)
    ivs = [iv for iv, mult in P.intervals(eps=sp.Rational(1, 10**13)) if iv[0] >= sp.Rational(lo.numerator, lo.denominator) and iv[1] <= sp.Rational(hi.numerator, hi.denominator)]
    assert len(ivs) == 1, ivs
    res['INF' if t is None else str(t)] = ivs[0]
    print('t =', 'INF' if t is None else t, ' right endpoint in', ivs[0], ' ~ %.13f' % float(ivs[0][0]))
bc = res['INF']
for k, iv in res.items():
    if k != 'INF':
        print(f't={k}: finite right endpoint > beta_c exactly: {iv[0] > bc[1]}')
