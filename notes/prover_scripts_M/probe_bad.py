import numpy as np, sys
from scan_n7pos import from_d6
from first_order import analyse
s = sys.argv[1]
A = from_d6(s); m = len(A)
print(A.astype(int))
r = analyse(A); print({k: r[k] for k in ['alpha', 'beta_c', 'k0', 'ReSp', 'pos']})
N = np.ones((m, m)) - np.eye(m) - A.T
bc = r['beta_c']
for db in [1e-4, 1e-3, 3e-3, 1e-2, 3e-2]:
    beta = bc + db
    C = np.eye(m) + beta * N
    a = np.linalg.det(C) * np.linalg.solve(C, np.ones(m))
    row = []
    for t in [1/(beta-2)*1.01, 10, 30, 100, 300, 1e3, 3e3, 1e4, 1e5, 1e6]:
        K = C + t * np.ones((m, m))
        row.append((round(t, 2), round(np.linalg.eigvals(K).real.min(), 7)))
    print(round(beta, 6), 'pos', bool((a > 0).all()), row)
