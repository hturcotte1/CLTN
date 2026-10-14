"""HEURISTIC ONLY. Given base digraphs (digraph6) on n nodes, enumerate all one-node extensions (2^(2n) ways),
keep domination-free ones, and float-test: positivity at beta + legal stability somewhere; report alpha.
Flags FIN_NOT_S (legal SFP somewhere on the grid while nu_inf < 1/2, i.e. alpha > -1/2): graph-level M violation candidate."""
import sys, numpy as np
sys.path.insert(0, '/home/user/CLTN')
from graphs import is_domination_free
from scan_n7pos import from_d6, process
from to_d6 import to_d6
import io
bases = sys.argv[1:]
out = io.StringIO()
tot = 0
for s in bases:
    B = from_d6(s).astype(int); n = len(B)
    batch = []
    for bits in range(1 << (2 * n)):
        A = np.zeros((n + 1, n + 1), int); A[:n, :n] = B
        for i in range(n):
            A[i, n] = bits >> i & 1
            A[n, i] = bits >> (n + i) & 1
        if not is_domination_free(A.tolist()): continue
        batch.append(to_d6(A.tolist()))
    tot += len(batch)
    for k in range(0, len(batch), 4000):
        process(batch[k:k+4000], out)
print('extensions tested (domination-free):', tot)
print(out.getvalue())
