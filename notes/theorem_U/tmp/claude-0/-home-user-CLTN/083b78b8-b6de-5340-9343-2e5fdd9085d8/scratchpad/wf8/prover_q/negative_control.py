"""negative_control.py -- exact: explicit limit test (E1) and positivity (E2) for G_k with k = k0-1 .. k0+1 on the two
small instances of certify_explicit.py, showing the reduced threshold k0 is attained and not vacuous.  No floats.
Run: python3 negative_control.py  (writes negative_control.json)"""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q')
from certify_explicit import build_Gk
from exact_core import psi, roots_right_of, positivity, strip_zero_roots
out = []
for c, k0, beta in (([1, 1, 5, 3, 1], 9, F(3)), ([1, 1, 4, 2, 1], 26, F(13, 4))):
    for k in (k0 - 1, k0, k0 + 1):
        A, Dk, names, m = build_Gk(c, k)
        p, z = strip_zero_roots(psi(Dk))
        lt = roots_right_of(p, 1 / beta)
        rec = {'c': c, 'beta': str(beta), 'k': k, 'n': len(A), 'limit_alpha_lt_-gamma': str(lt[0]), 'detail': lt[1],
               'positive': positivity(A, beta)['ok']}
        out.append(rec); print(rec, flush=True)
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/negative_control.json', 'w'), indent=1)
