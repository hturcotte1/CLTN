"""Exact certification of high-q instances.  For each graph: (1) exact limit check at a rational beta just below the
heuristic top (alpha < -gamma by exact Hurwitz, w > 0 by exact solve) and (2) an exact legal stable positive full-support
fixed point via ctln.analyze_support at explicit rational (eps, delta) with large t.  q = (beta-1) t/(t+1)."""
import json, sys
import numpy as np
from fractions import Fraction as F
import tlib, ctln
def heur_top(B):
    a = tlib.heur_alpha(B); bc = 1 / (1 + a)
    grid = np.linspace(2.0001, bc - 1e-9, 6000)
    ok = [b for b in grid if tlib.heur_w(B, b).min() > 0]
    return a, bc, (max(ok) if ok else None)
def certify(B, ts=(10**3, 10**4, 10**5, 10**6)):
    a, bc, top = heur_top(B)
    if top is None: return None
    for den in (10000, 2000, 500, 100):
        b = F(int(top * den), den)
        for _ in range(40):
            if b <= 2: break
            if tlib.limit_ok_at(B, b):
                w = tlib.sfp_witness(B, b, ts=ts)
                if w:
                    t, eps, delta, r = w
                    return {'d6': tlib.to_d6(B), 'n': len(B), 'heur_alpha': a, 'heur_beta_c': bc, 'heur_pos_top': top,
                            'beta': str(b), 'limit_ok_exact': True, 't': str(t), 'eps': str(eps), 'delta': str(delta),
                            'q': str(ctln.q_of(eps, delta)), 'q_float_display_only': float(ctln.q_of(eps, delta)),
                            'legal': ctln.legal(eps, delta), 'positive': r['positive'], 'stable': r['stable'],
                            'detail': r['stability_detail'], 'clique': r['clique'], 'x': r['x']}
            b -= F(1, den)
    return None
if __name__ == '__main__':
    d6s = json.load(open(sys.argv[1])); out = []
    for s in d6s:
        r = certify(tlib.from_d6(s))
        if r: out.append(r); print(r['n'], r['d6'], 'beta', r['beta'], 't', r['t'], 'q=%.5f' % r['q_float_display_only'], r['detail'], flush=True)
        else: print('no certificate for', s, flush=True)
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
