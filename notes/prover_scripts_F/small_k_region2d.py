"""Sanity cross-check (not a proof, relies on region2d.py's CAD): for every non-complete digraph
H on k = 2, 3, 4 nodes, the full support never satisfies (1)&(3) at any legal (beta, t).
Uses region2d with sigma = everything (so the o_k / condition (2) part is vacuous)."""
import sys, time, json
sys.path.insert(0, "/home/user/CLTN")
import graphs, region2d
from collections import Counter
out = {}
for k in (3, 4):
    keys = sorted(graphs.enumerate_digraphs(k))
    tally = Counter(); t0 = time.time()
    for key in keys:
        A = graphs.from_key(key, k)
        try:
            r = region2d.region2d(A, range(k), crosscheck=True)
            v = r["verdict"]
        except ArithmeticError as ex:
            v = "ERROR:" + str(ex)
            print("ERROR", k, A, ex, flush=True)
        tally[(graphs.is_complete_bidirected(A), v)] += 1
        if v != "never" and not graphs.is_complete_bidirected(A):
            print("NON-COMPLETE NOT NEVER:", k, A, v, flush=True)
    print(k, len(keys), dict(tally), f"{time.time()-t0:.1f}s", flush=True)
    out[k] = {str(a): b for a, b in tally.items()}
json.dump(out, open("small_k_region2d.json", "w"), indent=1)
