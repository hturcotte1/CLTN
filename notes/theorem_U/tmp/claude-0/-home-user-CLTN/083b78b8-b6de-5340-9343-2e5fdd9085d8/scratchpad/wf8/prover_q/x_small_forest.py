"""HEURISTIC search (floats only to rank candidates) for small arms-forest instances; every reported line is then
checked exactly by certify_forest.certify (C1-C3)."""
import sys, itertools, numpy as np
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q')
from certify_forest import certify, closed_form_positivity
res = []
for d in (4, 5, 6):
    for tail in itertools.product(range(0, 6), repeat=d - 1):
        c = [1, 1] + list(tail)
        if c[2] < 1: continue
        # forest constraint: c_j >= 1 needed if c_{j+2} >= 1 for j >= 3
        if any(c[j + 2] >= 1 and c[j] < 1 for j in range(3, d - 1)): continue
        if c[-1] == 0: continue
        m = 3 + 2 * (sum(c[2:]) - 1)
        if m > 40: continue
        mr = np.roots(c).real.min()
        if mr >= 0: continue
        ceil = -1 / mr
        if ceil < 2.6: continue
        res.append((ceil, m, c))
res.sort(key=lambda x: -x[0])
print(len(res), 'float candidates; top:', res[:5])
found = []
for ceil, m, c in res[:400]:
    for beta in (F(13, 4), F(3), F(14, 5), F(5, 2)):
        if beta >= ceil: continue
        if not closed_form_positivity(c, beta)['ok']: continue
        rec = certify(c, beta)
        if rec['certified'] and rec['n0'] <= 70:
            found.append((float(beta), rec['n0'], m, c, rec['C3_k0']))
            print('beta*', beta, 'n0', rec['n0'], 'm', m, 'k0', rec['C3_k0'], 'c', c, flush=True)
            break
found.sort(key=lambda x: (-x[0], x[1]))
print('BEST', found[:8])
