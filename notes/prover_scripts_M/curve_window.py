"""EXACT curve test (Lemma R) on the n=7 positivity survivors with edges >= 22 and float alpha in (-1/2, -0.3)
(float used only to SELECT; every verdict is exact). Writes every non-'never' with an exact alpha test."""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from curve1d import curve_test
from certify_M_violation import from_d6
from certify_graph_level_pure import omega_at
lo, hi, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rows = [l.split()[0] for l in open('alpha_window.txt')][lo:hi]
tally = {}; t0 = time.time()
with open(outp, 'w') as fh:
    for i, s in enumerate(rows):
        A = from_d6(s); r = curve_test(A); v = r['verdict']
        if v == 'stable_somewhere':
            hs = ctln.hurwitz_stable(omega_at(A, F(2)))[0]
            v += '|alpha<-1/2' if hs is True else ('|alpha>=-1/2' if hs is False else '|alpha boundary')
        if v != 'never':
            fh.write(json.dumps({"d6": s, "verdict": v, "stable_cells_t": [(c['lo'], c['hi'], c['sample_t']) for c in r.get('cells', []) if c['stable']]}) + "\n"); fh.flush()
        tally[v] = tally.get(v, 0) + 1
        if i % 2000 == 0: print(lo + i, tally, round(time.time() - t0), flush=True)
print('DONE', lo, hi, tally, round(time.time() - t0), flush=True)
