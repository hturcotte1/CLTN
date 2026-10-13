"""EXPLORATION: among n=5 graphs whose exact limit test chi_{A|perp}(z-1/2) is 'boundary only' (None, no False),
how many have alpha(A|perp) = -1/2 (float, 1e-9) and how many are tournaments / oriented?"""
import sys
sys.path.insert(0, '/home/user/CLTN')
sys.path.insert(0, '.')
import numpy as np, graphs, enum_cache
from collections import Counter
import importlib.util
spec = importlib.util.spec_from_file_location('xp', 'x_I7_prune.py')
src = open('x_I7_prune.py').read().split("keys5 = enum_cache.load_keys")[0]
ns = {}; exec(src, ns)
def helmert(m):
    Q = np.zeros((m, m-1))
    for k in range(1, m):
        Q[:k, k-1] = 1; Q[k, k-1] = -k; Q[:, k-1] /= np.sqrt(k*(k+1))
    return Q
H = helmert(5); c = Counter()
for v in enum_cache.load_keys(5, '/home/user/CLTN/certificates/digraphs_n5.keys'):
    A = graphs.from_key(enum_cache.key_from_int(v, 5), 5)
    if graphs.is_complete_bidirected(A): continue
    kill, bnd, lim, lbnd, h = ns['analyse'](A)
    if not lbnd: continue
    al = np.linalg.eigvals(H.T @ np.array(A, float) @ H).real.max()
    E = sum(map(sum, A)); orient = graphs.is_oriented(A)
    c[('alpha=-1/2' if abs(al + .5) < 1e-9 else f'alpha>-1/2' if al > -.5 else 'alpha<-1/2??',
       'tournament' if orient and E == 10 else 'oriented' if orient else 'mixed/sym')] += 1
for k, v in sorted(c.items()): print(k, v)
