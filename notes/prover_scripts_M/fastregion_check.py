"""EXACT: fastregion.region (repo engine) at t = 10, 30, 100, 1000 and t = INF for one graph; prints the stable
beta-intervals (isolating intervals of the algebraic endpoints) to show the finite-t stable interval extends
beyond the limit critical value beta_c (the right end of the t = INF interval)."""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import fastregion
from certify_M_violation import from_d6
A = from_d6(sys.argv[1]); sig = list(range(len(A)))
for t in [fastregion.INF, F(10), F(30), F(100), F(1000)]:
    cert = fastregion.region(A, sig, t, crosscheck=(t is not None))
    print('t =', 'INF' if t is None else t, 'verdict', cert['verdict'])
    for c in cert.get('stable_intervals', []):
        def show(e):
            if e['kind'] == 'beta_min': return 'beta_min=' + e['beta']
            if e['kind'] == 'inf': return 'inf'
            return 'root in [%s, %s] (=%.7f) deg %d %s' % (e['beta_interval'][0], e['beta_interval'][1],
                    float(F(e['beta_interval'][0])), len(e['beta_minpoly_hi_to_lo']) - 1, e['vanishing'])
        print('    (', show(c['lo']), ',', show(c['hi']), ')')
