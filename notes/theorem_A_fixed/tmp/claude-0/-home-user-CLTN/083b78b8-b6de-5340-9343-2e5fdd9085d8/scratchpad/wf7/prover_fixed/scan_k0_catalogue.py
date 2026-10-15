"""EXACT: for each catalogue family with (float-preselected) limit min Re > -1/2, compute at fixed points (beta, t)
 - positivity of every member (Lemma X2: w_G(beta) > 0 and sigma_S(beta) > 1, exact rationals)
 - k0 = smallest integer >= 0 with every Hurwitz minor (polynomial in k) > 0 on [k0, inf)  (exact Sturm counts)
Output: scan_k0_catalogue.json.  The float preselection only chooses which families to try."""
import json, sys, time, sympy as sp
from fractions import Fraction as F
import fixlib as fl
POINTS = [(F(21,10), F(20)), (F(11,5), F(20)), (F(23,10), F(20)), (F(12,5), F(100))]
good = json.load(open('heur_catalogue_good.json'))
res = []
t0 = time.time()
for i, o in enumerate(good):
    A = fl.from_d6(o['base']); S = tuple(o['S'])
    Psi = fl.Psi_expr(A, [S], [1])
    row = {'base': o['base'], 'S': list(S), 'chi_L': o['chi_L'], 'heur_minRe': o['minRe'], 'points': []}
    for beta, t in POINTS:
        pos, w, sig = fl.positivity_core(A, [S], beta)
        co = fl.shifted_coeffs(Psi.subs(fl.a_s, sp.Rational(t) / sp.Rational(beta)), 1 / sp.Rational(beta))
        mins = fl.hurwitz_minors(co)
        k0 = 0
        for mn in mins:
            kk = fl.smallest_k0(mn, 0)
            if kk is None: k0 = None; break
            k0 = max(k0, kk)
        row['points'].append({'beta': str(beta), 't': str(t), 'positive': pos, 'k0': k0})
    res.append(row)
    print(i, o['base'], o['S'], [(p['positive'], p['k0']) for p in row['points']], round(time.time() - t0), flush=True)
json.dump(res, open('scan_k0_catalogue.json', 'w'), indent=0)
