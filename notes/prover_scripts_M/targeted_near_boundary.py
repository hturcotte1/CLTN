"""HEURISTIC selection + EXACT decision. Targets the first-order mechanism of §4: n=7 positivity survivors with
alpha slightly above -1/2 (critical moving zero Re z0 in (0.44, 0.5)) and first-order coefficient R > 1/8.
Every selected graph is decided EXACTLY by the Lemma-R curve test; a 'stable_somewhere' with alpha >= -1/2 would
refute M_graph^df at n = 7."""
import sys, numpy as np, json, time
sys.path.insert(0, '/home/user/CLTN')
from zeros_R import zeros_R, compress
from curve1d import curve_test
from scan_n7pos import from_d6 as from_d6_np
from certify_M_violation import from_d6
path, start, stop, outp = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
U = compress(7)
t0 = time.time(); nsel = 0; ncurve = 0; hits = []; tot = 0
batch = []
def flush(batch, fh):
    global nsel, ncurve
    As = np.array([from_d6_np(s) for s in batch])
    P = np.einsum('ji,bjk,kl->bil', U, As, U)
    alpha = np.linalg.eigvals(P).real.max(axis=1)
    for b in np.nonzero((alpha > -0.5 + 1e-9) & (alpha < -0.44))[0]:
        zr = zeros_R(As[b])
        crit = [(z, R) for z, R, k in zr if k == 'moving' and R is not None and np.isfinite(R) and abs(z.real + alpha[b]) < 1e-7]
        if not crit: continue
        R = max(r for z, r in crit)
        if R <= 0.125: continue
        nsel += 1
        r = curve_test(from_d6(batch[b])); ncurve += 1
        rec = {"d6": batch[b], "alpha": float(alpha[b]), "R": float(R), "curve": r['verdict']}
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        if r['verdict'] != 'never': print('!!! ', rec, flush=True)
with open(outp, 'w') as fh:
    for i, line in enumerate(open(path)):
        if i < start: continue
        if i >= stop: break
        sp = line.split()
        if int(sp[1]) < 22: continue
        batch.append(sp[0]); tot += 1
        if len(batch) == 5000:
            flush(batch, fh); batch = []
            print(i, 'dense seen', tot, 'selected', nsel, round(time.time() - t0), flush=True)
    if batch: flush(batch, fh)
print('DONE dense seen', tot, 'selected (alpha in (-0.5,-0.44), R>1/8)', nsel, 'curve-tested', ncurve, round(time.time() - t0), flush=True)
