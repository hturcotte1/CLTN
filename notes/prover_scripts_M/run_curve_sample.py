import sys, time, random
from curve1d import curve_test
from certify_M_violation import from_d6
path, k, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
lines = []
with open(path) as fh:
    for i, l in enumerate(fh):
        if i >= 400000: break
        lines.append(l.split()[0])
random.seed(seed)
sample = random.sample(lines, k)
t0 = time.time(); tally = {}
for s in sample:
    r = curve_test(from_d6(s))
    tally[r['verdict']] = tally.get(r['verdict'], 0) + 1
    if r['verdict'] != 'never': print(s, r['verdict'], [(c['lo'][:8], (c['hi'] or 'inf')[:8], c['stable']) for c in r['cells']], flush=True)
dt = time.time() - t0
print('tally', tally, 'sec/graph', dt / k)
