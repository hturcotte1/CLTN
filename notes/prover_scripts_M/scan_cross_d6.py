"""HEURISTIC ONLY. Crossing directions for positive (graph, beta) pairs from a digraph6 stream."""
import sys, numpy as np
from crossings import crossings, a_vec
from scan_n7pos import from_d6
path = sys.argv[1]; start = int(sys.argv[2]); stop = int(sys.argv[3])
betas = np.concatenate([2 + np.logspace(-3, 1.3, 40), 1 + np.linspace(0.05, 0.95, 10)])
stats = {}; ex = []
npos = 0
with open(path) as fh:
    for i, line in enumerate(fh):
        if i < start: continue
        if i >= stop: break
        s = line.split()[0]
        A = from_d6(s); m = len(A)
        N = np.ones((m, m)) - np.eye(m) - A.T
        for beta in betas:
            a = a_vec(N, beta)
            if not (a > 1e-9).all(): continue
            npos += 1
            for (t, y, dl) in crossings(N, beta):
                leg = t * (beta - 2) > 1
                key = (beta > 2, leg, dl > 0)
                stats[key] = stats.get(key, 0) + 1
                if dl <= 0 and len(ex) < 20: ex.append((s, beta, t, y, dl, leg))
print('pos pairs', npos)
for k in sorted(stats): print('beta>2=%s legal=%s good=%s' % k, stats[k])
for e in ex: print(e)
