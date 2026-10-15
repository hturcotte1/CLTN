"""Exact instance certificates for family E (family_E.py): for a = 0..AMAX the explicit graph E_a is checked
(1) never limit-stable: l^2+l+1 divides n*chi_{A|perp} exactly (so alpha >= -1/2; Lemma I.3(d)), other roots Re < -1/2 (exact Hurwitz);
(2) a legal stable positive full-support fixed point at an explicit rational (eps, delta) via ctln.analyze_support.
Writes family_E_witnesses.json."""
import json, sys, time
from fractions import Fraction as F
import tlib, ctln, boundary_tools as bt, rstar
from family_E import E
AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
out = []
for av in range(0, AMAX + 1):
    t0 = time.time()
    A = E(av); n = len(A)
    mult, others = bt.boundary_split(A)
    R = rstar.R_at(A)[0]
    w = bt.witness_search(A, ts=(20, 30, 48, 60, 100, 200, 500, 1000), xs=(F(21, 20), F(11, 10), F(5, 4), F(3, 2), F(7, 4)))
    rec = {'a': av, 'n': n, 'd6': tlib.to_d6(A), 'A_rows': [''.join(map(str, r)) for r in A],
           'phi3_divides_n_chi_perp_exactly_once': mult == 1, 'other_roots_Re_lt_-1/2_exact': others, 'R_star_exact': str(R),
           'clique': ctln.is_clique(A, list(range(n))), 'witness': w}
    out.append(rec)
    print(av, n, rec['d6'], mult, others, R, None if not w else (w['t'], w['beta'], w['q'], w['detail']), '%.1fs' % (time.time() - t0), flush=True)
json.dump(out, open('family_E_witnesses.json', 'w'), indent=1)
