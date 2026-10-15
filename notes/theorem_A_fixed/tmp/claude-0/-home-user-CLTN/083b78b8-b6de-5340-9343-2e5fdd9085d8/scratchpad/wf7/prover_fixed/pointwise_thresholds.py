"""
pointwise_thresholds.py -- EXACT: at single fixed points, the smallest integer k0 such that every Hurwitz minor of
(-1)^(m+1) Psi_k(-1/beta - x) is > 0 for all real k >= k0 (positive leading coefficient, no real root >= k0: Sturm), together
with exact positivity (Lemma X2').  Conclusion per row: the family member is an SFP at that (eps, delta) for every integer k >= k0.
Writes pointwise_thresholds.json.
"""
import json
from fractions import Fraction as F
import fixlib as fl
from theorem_B import FAMILIES

POINTS = [('R', F(21, 10), F(20)), ('R', F(43, 20), F(7)), ('Q', F(21, 10), F(20)),
          ('Q2', F(13, 4), F(1)), ('Q2', F(7, 2), F(1)), ('Y', F(13, 4), F(1)), ('Y', F(4), F(1))]
out = []
for nm, beta, t in POINTS:
    base, classes, rhos = FAMILIES[nm]; A = fl.from_d6(base)
    eps, delta = fl.eps_delta(beta, t)
    assert fl.ctln.legal(eps, delta)
    pos, w, sig = fl.positivity_core(A, classes, beta)
    k0, mins = fl.k_threshold(A, classes, rhos, beta, t, kmin=0)
    n0 = len(A) + sum(rhos) * k0 if k0 is not None else None
    rec = {'family': nm, 'beta': str(beta), 't': str(t), 'eps': str(eps), 'delta': str(delta), 'q': str(fl.ctln.q_of(eps, delta)),
           'positive_all_members': pos, 'k0': k0, 'n_at_k0': n0,
           'statement': ('%s_k is an SFP at (eps, delta) = (%s, %s) for every integer k >= %s' % (nm, eps, delta, k0)) if (pos and k0 is not None) else 'no conclusion'}
    out.append(rec)
    print(rec['statement'], '(n >= %s)' % n0, flush=True)
json.dump(out, open('pointwise_thresholds.json', 'w'), indent=1)
