"""HEURISTIC (float) PREDICTION TABLE, not a certificate.  For each n=7 strict survivor whose limit verdict is
'stable_somewhere', the first-order coefficient R = Re(1/s'(z0)) at the critical moving zero z0 (Re z0 = -alpha).
Prediction (Lemma AS): the finite-t SFP beta-interval's right end is beta_c + beta_c^3 R / t + O(1/t^2);
R > 0  => T(beta) bounded above for beta slightly > beta_c (Lemma M fails), R < 0 => M holds near (beta_c, inf)."""
import json, numpy as np
from zeros_R import zeros_R
from scan_n7pos import from_d6
rows = []
for line in open('/home/user/CLTN/n7/n7_limit_results.jsonl'):
    r = json.loads(line)
    if r['verdict'] != 'stable_somewhere': continue
    A = from_d6(r['d6']); m = len(A)
    zr = zeros_R(A)
    bounded = [z for z, R, k in zr if k != 'pole']
    alpha = -min(z.real for z, R, k in zr)
    beta_c = 1 / (1 + alpha) if alpha > -1 + 1e-12 else float('inf')
    crit = [(z, R, k) for z, R, k in zr if abs(z.real + alpha) < 1e-7]
    kinds = sorted(set(k for z, R, k in crit))
    Rc = max((R for z, R, k in crit if k == 'moving' and R is not None and np.isfinite(R)), default=None)
    N = np.ones((m, m)) - np.eye(m) - A.T
    C = np.eye(m) + beta_c * N if np.isfinite(beta_c) else None
    pos = None
    if C is not None:
        a = np.linalg.det(C) * np.linalg.solve(C, np.ones(m)); pos = bool((a > 1e-9).all())
    rows.append({"d6": r['d6'], "complete": r.get('complete'), "alpha": alpha, "beta_c": beta_c, "critical_kinds": kinds,
                 "R": Rc, "beta_c3_R": (beta_c ** 3 * Rc) if (Rc is not None and np.isfinite(beta_c)) else None,
                 "positive_at_beta_c": pos,
                 "prediction": (None if Rc is None else ("T(beta) bounded above just right of beta_c (M fails)" if Rc > 0 else "M holds near (beta_c, inf)"))})
json.dump(rows, open('n7_first_order_table.json', 'w'), indent=1)
from collections import Counter
print(len(rows), Counter((r['positive_at_beta_c'], None if r['R'] is None else r['R'] > 0, tuple(r['critical_kinds'])) for r in rows))
