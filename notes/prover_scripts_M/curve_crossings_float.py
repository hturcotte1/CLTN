"""HEURISTIC (float). Along the legal boundary curve K_b(t) = I + (2+1/t) N + t J, count n_-(t) = #eigenvalues with
Re < 0 on a log grid. 'bad' = n_- increases somewhere (an eigenvalue leaves the stable half-plane as t grows);
'curve-M violation' = n_- goes from 0 to > 0. Compare dominated vs domination-free random graphs."""
import sys, numpy as np
sys.path.insert(0, '/home/user/CLTN')
from graphs import is_domination_free
rng = np.random.default_rng(int(sys.argv[1]))
ts = np.logspace(-2, 4, 1500)
st = {}; ex = []
for trial in range(int(sys.argv[2])):
    m = int(rng.integers(4, 10)); p = rng.uniform(0.5, 0.9)
    A = (rng.random((m, m)) < p).astype(int); np.fill_diagonal(A, 0)
    df = is_domination_free(A.tolist())
    N = np.ones((m, m)) - np.eye(m) - A.T
    neg = np.array([(np.linalg.eigvals(np.eye(m) + (2 + 1 / t) * N + t * np.ones((m, m))).real < -1e-10).sum() for t in ts])
    if (neg == 0).any():
        inc = np.diff(neg) > 0
        lost = ((neg[:-1] == 0) & inc).any()
        key = ('df' if df else 'dom', 'stable-somewhere', 'LOST' if lost else 'kept')
        if lost and df and len(ex) < 5: ex.append(A.tolist())
    else:
        key = ('df' if df else 'dom', 'never')
    st[key] = st.get(key, 0) + 1
for k in sorted(st): print(k, st[k])
for e in ex: print(e)
