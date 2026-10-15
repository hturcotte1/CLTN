import json, sys
from fractions import Fraction as F
import tlib, boundary_tools as bt, ctln
inp, outp = sys.argv[1], sys.argv[2]
keep = json.load(open(inp))
out = []
for c, d6, R in keep:
    A = tlib.from_d6(d6)
    Rf = F(R)
    # asymptotic window: beta = 2 + x/t with 1 < x < 8R (Theorem AS, c ~ t/2)
    xs = [F(21, 20), F(11, 10), F(5, 4), F(3, 2)] + [x for x in (F(7, 4), F(2), F(5, 2), F(3)) if x < 8 * Rf]
    w = bt.witness_search(A, xs=tuple(xs), ts=(20, 30, 48, 60, 100, 200, 500, 1000, 3000, 10**4))
    mult, others = bt.boundary_split(A)
    rec = {'canon': c, 'd6': d6, 'n': len(A), 'R_star': R, 'phi3_multiplicity': mult, 'other_roots_strictly_below_-1/2': others,
           'never_limit_stable_reason': 'l^2+l+1 divides n*chi_perp exactly, so alpha(A|perp) >= -1/2 (Lemma I.3(d): never limit-stable at legal beta)',
           'witness': w}
    out.append(rec)
    print(d6, 'R*', R, 'witness:', None if not w else (w['t'], w['beta'], w['q'], w['detail']), flush=True)
json.dump(out, open(outp, 'w'), indent=1)
