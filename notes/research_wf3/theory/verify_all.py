"""
verify_all.py -- re-runs every EXACT instance claim of notes_theory.md from the certificate JSON files.
No floating point in any decision.  Usage:  python3 verify_all.py
For every witness: rebuild the graph from its digraph6 string, re-run ctln.analyze_support at the stored rational
(eps, delta), and assert legal, x > 0, stable (all Hurwitz minors positive), not a clique.
For every limit claim: tlib.limit_ok_at at the stored beta.  For every boundary (finite-eps-only) claim: l^2+l+1 divides
n*chi_perp exactly (alpha >= -1/2, never limit-stable) and the remaining factor is Hurwitz at shift -1/2.
"""
import json, sys
from fractions import Fraction as F
import tlib, ctln, boundary_tools as bt
def check_point(d6, eps, delta, label):
    A = tlib.from_d6(d6); n = len(A); eps = F(eps); delta = F(delta)
    assert ctln.legal(eps, delta), label
    r = ctln.analyze_support(A, list(range(n)), eps, delta)
    assert r['positive'] and r['stable'] is True and r['stable_fixed_point'] and not r['clique'] and not r['target_free'], (label, r['stability_detail'])
    return ctln.q_of(eps, delta)
n_ok = 0
# family F witnesses
cf = json.load(open('family_F_certificate.json'))
for w in cf['F6_witnesses']:
    check_point(w['d6'], w['eps'], w['delta'], 'F n=%d' % w['n'])
    assert tlib.limit_ok_at(tlib.from_d6(w['d6']), F(w['beta'])); n_ok += 1
print('family F witnesses OK:', len(cf['F6_witnesses']))
# family E witnesses + boundary
ce = json.load(open('family_E_witnesses.json'))
for r in ce:
    w = r['witness']; check_point(r['d6'], w['eps'], w['delta'], 'E a=%d' % r['a'])
    mult, others = bt.boundary_split(tlib.from_d6(r['d6'])); assert mult == 1 and others is True; n_ok += 1
print('family E witnesses + exact alpha = -1/2 OK:', len(ce))
# n=8 boundary candidates
cb = json.load(open('cube_n8_witnesses.json'))
for r in cb:
    w = r['witness']; check_point(r['d6'], w['eps'], w['delta'], 'cube ' + r['d6'])
    mult, others = bt.boundary_split(tlib.from_d6(r['d6'])); assert mult == 1 and others is True; n_ok += 1
print('n=8 finite-eps-only boundary candidates OK:', len(cb))
# high-q instances
best = {}
for f in ['top_n8.json', 'cert_q_910.json', 'cert_q_11.json'] + (['cert_q_12.json'] if __import__('os').path.exists('cert_q_12.json') else []):
    for r in json.load(open(f)):
        if 'eps' not in r: continue
        q = check_point(r['d6'], r['eps'], r['delta'], 'q ' + r['d6'])
        B = tlib.from_d6(r['d6']); assert tlib.limit_ok_at(B, F(r.get('beta_exact', r.get('beta'))))
        n = len(B); best[n] = max(best.get(n, F(0)), q); n_ok += 1
print('high-q instances OK; best certified q per n:', {n: (str(q), float(q)) for n, q in sorted(best.items())})
# family Q / Q2 / Y witnesses (large n: these take a few minutes)
import os
for f, key in [('family_Q_certificate.json', 'Q5_witness'), ('family_Q2_certificate.json', 'W5_witnesses'), ('family_Y_certificate.json', 'Y5_witness')]:
    if not os.path.exists(f): print('missing', f); continue
    c = json.load(open(f)); ws = c.get(key)
    if ws is None: print('no instance in', f); continue
    if isinstance(ws, dict): ws = [ws]
    for w in ws:
        q = check_point(w['d6'], w['eps'], w['delta'], f + ' n=%d' % w['n'])
        assert tlib.limit_ok_at(tlib.from_d6(w['d6']), F(w['beta'])); n_ok += 1
        best[w['n']] = max(best.get(w['n'], F(0)), q)
    print(f, 'instances OK:', len(ws), flush=True)
print('best certified q per n (all records):', {n: (str(q), float(q)) for n, q in sorted(best.items())})
print('ALL CHECKS PASSED (%d records)' % n_ok)
