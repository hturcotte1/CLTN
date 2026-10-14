"""Cross-validation: exact Lemma-R curve test (curve1d) vs the float graph-level scan, all 9,608 digraphs on 5 nodes,
plus exact pointwise witnesses for every 'stable_somewhere' (ctln.analyze_support at a legal point near the curve)."""
import sys, numpy as np, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from curve1d import curve_test
from run_scan_small import iso
from scan_graph_level import run
As = np.array(list(iso(5)))
best, argc, ni = run(As)
float_st = best > 1e-9
t0 = time.time(); ex = []; nwit = 0; tally = {}
for b, A in enumerate(As.astype(int).tolist()):
    r = curve_test(A)
    tally[r['verdict']] = tally.get(r['verdict'], 0) + 1
    exact_st = r['verdict'] == 'stable_somewhere'
    if exact_st != float_st[b]: ex.append((b, r['verdict'], best[b]))
    if exact_st:
        # exact pointwise witness: stable sample t on the curve, then step slightly into the legal region
        c = next(c for c in r['cells'] if c['stable'])
        t = F(c['sample_t']); beta = 2 + 1 / t + F(1, 10**6)
        eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
        rec = ctln.analyze_support(A, list(range(5)), eps, delta)
        assert ctln.legal(eps, delta) and rec['stable'] is True, (A, t)
        nwit += 1
print('exact tally', tally, ' float legal-stable', int(float_st.sum()), ' disagreements', len(ex), ex[:10])
print('exact legal stable witnesses verified with ctln.analyze_support:', nwit, ' sec/graph', (time.time() - t0) / len(As))
