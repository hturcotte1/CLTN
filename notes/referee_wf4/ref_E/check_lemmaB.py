"""Referee: Lemma B conclusion tested exactly on E_a for a in {0, 2, 5}.
For c in CS: gamma_mid = 1/2 + (3/16)/c must be legal + positive + stable (exact ctln.analyze_support);
bracketing at the largest c (consistency with R* = 1/4 and the legality drift 1/8, exact verdicts):
  gamma = 1/2 + (1/4 + 1/50)/c  -> unstable (Hurwitz False), gamma = 1/2 + (1/4 - 1/50)/c -> stable,
  gamma = 1/2 + (1/8 - 1/50)/c  -> illegal,                    gamma = 1/2 + (1/8 + 1/50)/c -> legal.
Parameters are rational: beta = 1/(1-gamma), t = beta(c-1), eps = 1/(1+t), delta = eps(beta-1)."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E')
sys.path.insert(0, '/home/user/CLTN')
from fractions import Fraction as F
from reflib import d6_decode, blowup
import ctln
H = d6_decode('&FNDL~[uds?')
CS = [F(20), F(100), F(1000), F(10000)]
def params(c, g):
    beta = 1/(1-g); t = beta*(c-1); eps = 1/(1+t); delta = eps*(beta-1); return beta, t, eps, delta
out = []
for av in (0, 2, 5):
    A, _ = blowup(H, [1, 1, 2, av, 1, 1, 1]); n = len(A)
    for c in CS:
        g = F(1, 2) + F(3, 16)/c
        beta, t, eps, delta = params(c, g)
        r = ctln.analyze_support(A, list(range(n)), eps, delta)
        rec = {'a': av, 'n': n, 'c': str(c), 'gamma': str(g), 'beta': str(beta), 'eps': str(eps), 'delta': str(delta),
               'q': str(ctln.q_of(eps, delta)), 'legal': ctln.legal(eps, delta), 'positive': r['positive'],
               'stable': r['stable'], 'detail': r['stability_detail'], 'sfp': r['stable_fixed_point'], 'target_free': r['target_free']}
        out.append(rec); print(json.dumps(rec), flush=True)
    c = CS[-1]
    br = {}
    for name, g in (('R_hi_unstable', F(1, 2) + (F(1, 4) + F(1, 50))/c), ('R_lo_stable', F(1, 2) + (F(1, 4) - F(1, 50))/c),
                    ('L_lo_illegal', F(1, 2) + (F(1, 8) - F(1, 50))/c), ('L_hi_legal', F(1, 2) + (F(1, 8) + F(1, 50))/c)):
        beta, t, eps, delta = params(c, g)
        r = ctln.analyze_support(A, list(range(n)), eps, delta)
        br[name] = {'legal': ctln.legal(eps, delta), 'stable': r['stable'], 'detail': r['stability_detail'], 'positive': r['positive']}
    rec = {'a': av, 'bracket_c': str(c), 'bracket': br}
    out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_E/check_lemmaB.json', 'w'), indent=1)
