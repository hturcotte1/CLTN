import numpy as np, sys
from first_order import analyse
rng = np.random.default_rng(int(sys.argv[1]))
st = {}; ex = []
for trial in range(int(sys.argv[2])):
    m = int(rng.integers(3, 11))
    p = rng.uniform(0.5, 1.0)
    A = (rng.random((m, m)) < p).astype(float); np.fill_diagonal(A, 0)
    r = analyse(A)
    if r is None: continue
    if abs(r['S']) > 1e-7:
        key = ('fixed-eig', r['pos'])
    elif abs(r['k0'].imag) < 1e-9:
        key = ('real', r['pos'])
    else:
        key = ('cplx', 'pos' if r['pos'] else 'nonpos', 'ReSp<=0(good)' if r['ReSp'] <= 0 else 'ReSp>0(BAD)')
    st[key] = st.get(key, 0) + 1
    if key[0] == 'cplx' and key[1] == 'pos' and r['ReSp'] > 0 and len(ex) < 5: ex.append((A.astype(int).tolist(), r))
for k in sorted(st, key=str): print(k, st[k])
for e in ex: print(e)
