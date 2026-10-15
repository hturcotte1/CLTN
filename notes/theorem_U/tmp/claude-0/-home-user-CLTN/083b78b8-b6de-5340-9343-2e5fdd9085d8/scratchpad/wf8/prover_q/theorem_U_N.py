"""theorem_U_N.py N [Wmax] -- same certification as theorem_U.py for a single N (writes theorem_U_instance_N<N>.json).
No floating point."""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q')
from theorem_U import rounded
from certify_forest import certify, closed_form_positivity
from forest_exact import roots_right_of
N = int(sys.argv[1]); Wmax = int(sys.argv[2]) if len(sys.argv) > 2 else 50
t0 = time.time()
for W in range(1, Wmax + 1):
    c = rounded(N, W)
    if roots_right_of(c, F(1, N))[0] is not True: print('W=%d: C1 fails' % W, flush=True); continue
    if not closed_form_positivity(c, F(N))['ok']: print('W=%d: C2 fails' % W, flush=True); continue
    rec = certify(c, F(N)); rec['N'] = N; rec['W'] = W
    print('W=%d certified=%s m=%s k0=%s n0=%s %.1fs' % (W, rec['certified'], rec['m_core'], rec['C3_k0'], rec['n0'], time.time() - t0), flush=True)
    if rec['certified']:
        json.dump(rec, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/theorem_U_instance_N%d.json' % N, 'w'), indent=1)
        break
