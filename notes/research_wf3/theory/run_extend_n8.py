import numpy as np, json, sys, time
sys.path.insert(0, '.')
import heur_extend as he, load_hits, tlib
bases = [(r['d6'], r['A']) for r in load_hits.n7_hits()] + [('&FNFL~[udo?', tlib.from_d6('&FNFL~[udo?'))]
betas = np.concatenate([np.linspace(2.0005, 2.6, 300), np.linspace(2.61, 4.0, 140)])
w3 = np.exp(2j * np.pi / 3)
best = []; cube = []
t0 = time.time()
for bi, (d6, A) in enumerate(bases):
    B, pats = he.extensions(A)
    ev = he.alphas(B)
    a, q = he.qlim(B, ev, betas)
    for k in np.where(q > 0)[0]:
        best.append((float(q[k]), d6, pats[k].tolist(), float(a[k])))
    near = np.where(np.abs(a + 0.5) < 1e-7)[0]
    for k in near:
        if np.min(np.abs(ev[k] - w3)) < 1e-6 or np.min(np.abs(ev[k] - w3.conjugate())) < 1e-6:
            cube.append((d6, pats[k].tolist()))
    if bi % 10 == 0: print(bi, len(best), len(cube), '%.0fs' % (time.time() - t0), flush=True)
best.sort(key=lambda x: -x[0])
json.dump({'best': best[:3000], 'n_hits': len(best), 'cube': cube}, open('ext_n8.json', 'w'))
print('n8 one-node extensions with q_lim>0:', len(best), ' top:', best[:10])
print('cube-root boundary candidates:', len(cube))
