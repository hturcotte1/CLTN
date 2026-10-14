"""HEURISTIC (float): histogram of alpha(A|perp) and of the largest first-order coefficient R among critical zeros,
over the dense (edges >= 22) n=7 positivity survivors (all four parts), to size a targeted exact search."""
import sys, numpy as np, gzip
from zeros_R import compress
from scan_n7pos import from_d6
U = compress(7)
P = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/n7pos/n7_pos_all_PU.txt.gz'
edges_min = 22
hist = np.zeros(40, int); batch = []; tot = 0
sel = open('alpha_window.txt', 'w')
def flush(batch):
    global tot
    As = np.array([from_d6(s) for s in batch])
    al = np.linalg.eigvals(np.einsum('ji,bjk,kl->bil', U, As, U)).real.max(axis=1)
    idx = np.clip(((al + 1.0) / 0.05).astype(int), 0, 39)
    np.add.at(hist, idx, 1)
    for b in np.nonzero((al > -0.5 + 1e-9) & (al < -0.3))[0]:
        sel.write(f"{batch[b]} {al[b]:.6f}\n")
    tot += len(batch)
for line in gzip.open(P, 'rt'):
    sp = line.split()
    if int(sp[1]) < edges_min: continue
    batch.append(sp[0])
    if len(batch) == 20000: flush(batch); batch = []
if batch: flush(batch)
sel.close()
print('dense total', tot)
for k in range(40):
    if hist[k]: print(f"alpha in [{-1 + 0.05 * k:.2f}, {-1 + 0.05 * (k + 1):.2f}): {hist[k]}")
