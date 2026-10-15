"""EXPLORATION (floats): pole placement to (z+s)^m.  For each core D (canonical digraphs on m nodes), the unique theta(s)
with chi_{D - theta 1^T} = (z+s)^m (if (1^T, D) observable); report cores with theta(s) >= 0 for some s <= 1/m,
together with positivity at r = s (u > 0, sum u > 1, supp(theta) light)."""
import sys, numpy as np
sys.path.insert(0, '/home/user/CLTN')
from enum_cache import load_keys, key_from_int
from graphs import from_key
from math import comb
def theta_poly(D):
    m = len(D)
    O = np.array([np.ones(m) @ np.linalg.matrix_power(D, k) for k in range(m)])
    if abs(np.linalg.det(O)) < 1e-9: return None
    cD = np.poly(D)[1:]          # chi_D coefficients c_1..c_m
    # chi_L = chi_D + sum_j theta_j q_j ; q_j coefficient vector: [1^T B_k e_j], B_k = D^k + c1 D^{k-1} + ... + ck I
    Bs = [np.eye(m)]
    for k in range(1, m): Bs.append(D @ Bs[-1] + cD[k - 1] * np.eye(m))
    Qm = np.array([np.ones(m) @ B for B in Bs])   # m x m, row k = coefficient of z^{m-1-k}
    return Qm, cD
def scan(m, keys, sgrid):
    out = []
    for v in keys:
        A = np.array(from_key(key_from_int(v, m), m), float); D = 1 - np.eye(m) - A
        if D.sum() == 0: continue
        r = theta_poly(D)
        if r is None: continue
        Qm, cD = r
        try: Qi = np.linalg.inv(Qm)
        except np.linalg.LinAlgError: continue
        for s in sgrid:
            tgt = np.array([comb(m, k) * s ** k for k in range(1, m + 1)])
            th = Qi @ (tgt - cD)
            if th.min() >= -1e-12:
                b = 1 / s * 0.999
                try: u = np.linalg.solve(D.T + np.eye(m) / b, np.ones(m))
                except np.linalg.LinAlgError: continue
                supp = th > 1e-9
                pos = u.min() > 0 and u.sum() > 1 and all(u[supp] < 1)
                out.append((s, v, pos, th)); break
    return out
if __name__ == '__main__':
    m = int(sys.argv[1]); keys = load_keys(m, '/home/user/CLTN/certificates/digraphs_n%d.keys' % m) if m == 6 else None
    if keys is None:
        from graphs import enumerate_digraphs
        from enum_cache import key_to_int
        keys = [key_to_int(k) for k in enumerate_digraphs(m)]
    sgrid = np.linspace(1e-3, 1.0 / m, 400)
    out = scan(m, keys, sgrid)
    out.sort(key=lambda x: x[0])
    print(m, 'graphs', len(keys), 'with theta(s)>=0 for some s<=1/m:', len(out), 'of which positive:', sum(1 for o in out if o[2]))
    for s, v, pos, th in [o for o in out if o[2]][:10]:
        print('  s=%.4f ceiling>=%.2f key=%d theta=%s' % (s, 1 / s, v, np.round(th, 3).tolist()))
