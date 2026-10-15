"""Exact cross-check: where the certified limit interval ends on the Hurwitz curve (vanishing Delta), its right end
must be beta_c = 1/(1+alpha) computed here: gcd(record minpoly, our beta_c minpoly) has a root in the record's
isolating interval."""
import json, sympy as sp
x = sp.symbols('x')
D = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf5/bulge/'
L = [json.loads(l) for l in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
T = {r['d6']: r for r in json.load(open(D + 'rstar_table.json'))}
ok = bad = 0
for r in L:
    for li in (r['limit_intervals'] or []):
        hi = li['hi']
        if hi.get('vanishing') == ['Delta']:
            g = sp.gcd(sp.Poly(hi['beta_minpoly_hi_to_lo'], x), sp.Poly(T[r['d6']]['beta_c_minpoly_hi_to_lo'], x))
            a, b = map(sp.Rational, hi['beta_interval'])
            if g.degree() > 0 and g.count_roots(a, b) == 1: ok += 1
            else: bad += 1; print('MISMATCH', r['d6'])
print('Delta-ending limit intervals whose right end equals our beta_c exactly:', ok, 'mismatches:', bad)
