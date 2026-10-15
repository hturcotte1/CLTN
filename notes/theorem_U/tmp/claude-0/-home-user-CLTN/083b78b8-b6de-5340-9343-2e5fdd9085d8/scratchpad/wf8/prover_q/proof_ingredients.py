"""
proof_ingredients.py -- EXACT spot-checks of the quantitative steps in the proof of Theorem U (notes_prover_q.md, sec. 4),
at the explicit Rouche bound W = W0(N) = 2N 4^N (N+1)^(2N-2), for N = 3..Nmax.  No floating point.
 (a) c~_1 = 1 and c~_j >= 1 for 2 <= j <= 2N (so the arms forest exists and satisfies the forest constraint);
 (b) Lemma S (strip): every root of p~ has  -3s/2 < Re z < -s/2,  s = 1/(2N)   [two exact Hurwitz tests];
 (c) Lemma P positivity of the arms-forest core at beta = N  (closed form; asymptotic claim of step 3).
These are consistency checks of the proof, not part of it (the proof is valid for every N with no computation).
Run: python3 proof_ingredients.py [Nmax]   (writes proof_ingredients.json)
"""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q')
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/notes/research_wf3/theory')
import ctln
from tlib import taylor_shift
from theorem_U import rounded
from forest_exact import roots_right_of
from certify_forest import closed_form_positivity
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
out = []
for N in range(3, Nmax + 1):
    t0 = time.time()
    W = 2 * N * 4 ** N * (N + 1) ** (2 * N - 2); s = F(1, 2 * N)
    c = rounded(N, W)
    a = c[0] == 1 and c[1] == 1 and all(x >= 1 for x in c[2:])
    left = roots_right_of(c, 3 * s / 2)                                   # all Re z > -3s/2
    right = ctln.hurwitz_stable(taylor_shift([F(x) for x in c], -s / 2))   # p~(w - s/2) Hurwitz <=> all Re z < -s/2
    pos = closed_form_positivity(c, F(N))
    rec = {'N': N, 'W0': W, 'digits_c_d': len(str(c[-1])), 'a_coeffs_ok': a, 'b_Re_gt_-3s/2': str(left[0]),
           'b_Re_lt_-s/2': str(right[0]), 'c_positive_at_beta_N': pos['ok'], 'lambda_fmax': str(pos.get('lambda_fmax')), 'lambda_fmax_lt_1': pos['lambda_fmax'] < 1,
           'u_c2': str(pos.get('u_c2')), 'u_c2_in_(0,1)': 0 < pos['u_c2'] < 1, 'seconds': round(time.time() - t0, 1)}
    rec['pass'] = a and left[0] is True and right[0] is True and pos['ok']
    out.append(rec); print({k: rec[k] for k in rec if k not in ('lambda_fmax', 'u_c2')}, flush=True)
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf8/prover_q/proof_ingredients.json', 'w'), indent=1)
