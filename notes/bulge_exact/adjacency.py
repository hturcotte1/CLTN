"""Exact check: for every graph with a bounded stable cell, the left end of each bounded cell is the critical
point beta_c = 1/(1+alpha) (gcd of the region2d critical minpoly with the minpoly of beta_c has a root inside the
critical isolating interval), and no stable cell lies left of beta_c with bounded t.  Also, display-only, the
first-order prediction t_hi ~ beta_c^3 R*/(beta_s - beta_c) at the cell sample."""
import json
from fractions import Fraction as F
import sympy as sp
x = sp.symbols('x')
D = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf5/bulge/'
L = [json.loads(l) for l in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
T = {r['d6']: r for r in json.load(open(D + 'rstar_table.json'))}
def P(c): return sp.Poly(c, x)
out = []
for r in L:
    if not r['bounded_in_t_somewhere']: continue
    t = T[r['d6']]
    Mb = P(t['beta_c_minpoly_hi_to_lo'])
    crit = r['region2d'].get('beta_critical') or r.get('beta_critical') or []
    cells = r['region2d']['cells']
    res = []
    for c in cells:
        bad = [tc for tc in c['t_cells'] if tc['stable_fp'] and tc['t_hi'] is not None]
        if not bad: continue
        lo = F(c['beta_cell'][0])
        # critical point immediately left: isolating interval with right end <= lo, maximal
        left = [k for k in crit if F(k['interval'][1]) <= lo]
        if not left:   # cell starts at the domain boundary beta = 2
            ok = (lo == 2 and t['beta_c_minpoly_hi_to_lo'] == [1, -2])
        else:
            k = max(left, key=lambda k: F(k['interval'][1]))
            a, b = sp.Rational(k['interval'][0]), sp.Rational(k['interval'][1])
            g = sp.gcd(P(k['minpoly_hi_to_lo']), Mb)
            ok = g.degree() > 0 and sp.Poly(g, x).count_roots(a, b) == 1
        a0, b0 = F(t['beta_c_interval'][0]), F(t['beta_c_interval'][1])
        while b0 - a0 > F(1, 10**18):     # exact bisection (display value only)
            mid = (a0 + b0) / 2
            if Mb.count_roots(sp.Rational(a0.numerator, a0.denominator), sp.Rational(mid.numerator, mid.denominator)) >= 1: b0 = mid
            else: a0 = mid
        bc = (a0 + b0) / 2
        bs = F(c['beta_sample'])
        pred = t['beta_c3_Rstar_display_only'] / float(bs - bc) if (t['beta_c3_Rstar_display_only'] and bs != bc) else None
        res.append({'beta_cell': c['beta_cell'], 'left_critical_is_beta_c': ok, 't_hi_at_sample': [tc['t_hi'] for tc in bad],
                    'display_only_first_order_t_hi': pred, 'beta_c_display': float(bc), 'beta_sample_display': float(bs)})
    out.append({'d6': r['d6'], 'bounded_cells': res})
json.dump(out, open(D + 'bulge_adjacency.json', 'w'), indent=1)
allok = all(c['left_critical_is_beta_c'] for o in out for c in o['bounded_cells'])
print(len(out), 'graphs; bounded cells', sum(len(o['bounded_cells']) for o in out), 'all left ends = beta_c exactly:', allok)
for o in out[:40]:
    for c in o['bounded_cells']:
        print(o['d6'], c['left_critical_is_beta_c'], c['t_hi_at_sample'], '%.1f' % c['display_only_first_order_t_hi'] if c['display_only_first_order_t_hi'] else None)
