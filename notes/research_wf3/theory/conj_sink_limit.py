"""HEURISTIC test of: for a sink-class extension G[S, all; k] that is positive at beta (w_G > 0 and sigma_S > 1),
the k -> infinity limit matrix L = D - b 1^T has an eigenvalue with Re <= -1/2.
Random digraphs (m = 4..8) and all S; floats only."""
import numpy as np, itertools, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
viol = []; tested = 0
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 4000):
    m = rng.integers(4, 9)
    p = rng.uniform(0.5, 0.9)
    A = (rng.random((m, m)) < p).astype(int); np.fill_diagonal(A, 0)
    D = 1 - np.eye(m, dtype=int) - A
    for beta in (2.05, 2.3, 2.7, 3.2):
        g = 1 - 1 / beta
        try: w = np.linalg.solve(A.T + g * np.eye(m), np.ones(m))
        except np.linalg.LinAlgError: continue
        if w.min() <= 1e-9: continue
        for bits in range(1, 1 << m):
            b = np.array([(bits >> i) & 1 for i in range(m)])
            S = 1 - b
            if (S * w).sum() <= 1 + 1e-9: continue
            tested += 1
            L = D - np.outer(b, np.ones(m))
            mr = np.linalg.eigvals(L).real.min()
            if mr > -0.5 + 1e-7:
                viol.append((A.tolist(), b.tolist(), beta, mr))
print('tested (G,S,beta) positive:', tested, ' violations (min Re spec L > -1/2):', len(viol))
for v in viol[:5]: print(v)
