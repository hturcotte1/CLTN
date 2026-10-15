"""HEURISTIC: grow cores by one node (random in/out patterns), keep the best many-sink-limit score (core_design.score)."""
import numpy as np, sys, time
import tlib, core_design as cd
cd.BGRID = [2.5, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0]
rng = np.random.default_rng(int(sys.argv[2]))
cores = [np.array(tlib.from_d6(s)) for s in sys.argv[1].split(',')]
samples = int(sys.argv[3]); levels = int(sys.argv[4])
for L in range(levels):
    cand = []
    t0 = time.time()
    for A in cores:
        m = len(A)
        for s in range(samples):
            B = np.zeros((m + 1, m + 1), int); B[:m, :m] = A
            B[:m, m] = rng.integers(0, 2, m); B[m, :m] = rng.integers(0, 2, m)
            sc, info = cd.score(B)
            if sc > 0: cand.append((sc, tlib.to_d6(B.tolist()), info))
    cand.sort(key=lambda x: -x[0])
    print('level m=%d: best scores' % (len(cores[0]) + 1), [c[0] for c in cand[:8]], cand[:2], '%.0fs' % (time.time() - t0), flush=True)
    cores = [np.array(tlib.from_d6(c[1])) for c in cand[:6]]
