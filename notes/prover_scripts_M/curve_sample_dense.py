"""EXACT curve test (Lemma R) on a random sample of n=7 positivity survivors with edges >= 22 (necessary: B > Nn).
For each legal-stable graph also the exact alpha test (omega(.;2) Hurwitz). Output JSON lines for non-'never'."""
import sys, random, json, time, glob
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from curve1d import curve_test
from certify_M_violation import from_d6
from certify_graph_level_pure import omega_at
k, seed, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
P = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/n7pos/'
lines = []
for f in ['n7_pos_part0.txt', 'n7_pos_part1.txt', 'n7_pos_part2.txt', 'n7_pos_part3.txt']:
    for l in open(P + f):
        sp = l.split()
        if len(sp) >= 2 and int(sp[1]) >= 22: lines.append(sp[0])
random.seed(seed)
sample = random.sample(lines, k)
print('pool', len(lines), 'sample', k, flush=True)
tally = {}; t0 = time.time()
with open(outp, 'w') as fh:
    for i, s in enumerate(sample):
        A = from_d6(s)
        r = curve_test(A); v = r['verdict']
        if v == 'stable_somewhere':
            hs = ctln.hurwitz_stable(omega_at(A, F(2)))[0]
            v += '|alpha<-1/2' if hs is True else ('|alpha>=-1/2 !!' if hs is False else '|alpha boundary')
        if v != 'never':
            fh.write(json.dumps({"d6": s, "verdict": v, "stable_cells_t": [(c['lo'], c['hi'], c['sample_t']) for c in r.get('cells', []) if c['stable']]}) + "\n"); fh.flush()
        tally[v] = tally.get(v, 0) + 1
        if i % 1000 == 0: print(i, tally, round(time.time() - t0), flush=True)
print('DONE', tally, round(time.time() - t0), flush=True)
