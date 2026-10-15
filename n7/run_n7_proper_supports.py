"""
run_n7_proper_supports.py -- the complete list of stable NON-CLIQUE PROPER supports at n = 7.

By Theorem F (notes/lemma_F_reduction.md) a proper support sigma of a 7-node graph is a stable fixed point
iff the induced subgraph on sigma has a stable positive full support AND the off-support condition (2)
holds; the only non-clique 6-node graphs with a stable positive full support are the three n = 6
counterexamples (complete map), and there are none on <= 5 nodes. So every stable non-clique proper
support at n = 7 is a size-6 support inducing one of the three n = 6 hits. This script enumerates every
7-node graph obtained from each hit (on nodes 0..5) by adding node 6 with any in/out pattern (3 x 4096,
deduplicated up to isomorphism by the project's canonical form), and certifies with region2d (exact
two-parameter decomposition including the off-support polynomial o_6) whether the support {0..5} is
stable somewhere in the legal plane. Output: certificates/n7_proper_supports_shard<k>of<K>.jsonl.
No floating point anywhere.
"""
import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from region2d import region2d
from graphs import canonical_form
from enum_cache import key_to_int
from n7.curve_test import from_d6

HITS = {"geneson_idx1000255": "&EZjtwWE"}
def load_hits():
    hits = dict(HITS)
    for f, name in [("certificates/COUNTEREXAMPLE_CANDIDATE_n6_idx1387528.json", "idx1387528"), ("certificates/COUNTEREXAMPLE_CANDIDATE_n6_idx1424629.json", "idx1424629")]:
        hits[name] = json.load(open(f))["A_ours"]
    out = {}
    for k, v in hits.items():
        out[k] = from_d6(v) if isinstance(v, str) else [[int(x) for x in row] for row in v]
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--shard", default="0/1"); ap.add_argument("--out", default="certificates")
    a = ap.parse_args(); k, K = map(int, a.shard.split("/"))
    hits = load_hits(); sigma = [0, 1, 2, 3, 4, 5]
    cases = []
    seen = set()
    for name, H in hits.items():
        for pat in range(4096):
            A = [row[:] + [(pat >> i) & 1] for i, row in enumerate(H)] + [[(pat >> (6 + j)) & 1 for j in range(6)] + [0]]
            key = key_to_int(canonical_form(A)[0])
            if key in seen:
                continue
            seen.add(key); cases.append((name, pat, key, A))
    cases.sort(key=lambda c: c[2])
    mine = [c for i, c in enumerate(cases) if i % K == k]
    print(f"{len(cases)} distinct graphs; shard {a.shard} has {len(mine)}", flush=True)
    base = os.path.join(a.out, f"n7_proper_supports_shard{k}of{K}")
    tally = {"graphs": 0, "never": 0, "stable_somewhere": 0, "escalate": 0, "errors": 0, "distinct_total": len(cases)}
    t0 = time.time()
    with open(base + ".jsonl", "w") as fh:
        for name, pat, key, A in mine:
            tally["graphs"] += 1
            try:
                c = region2d(A, sigma, crosscheck=True, replay=True)
            except Exception as ex:
                tally["errors"] += 1; fh.write(json.dumps({"base": name, "pattern": pat, "key": key, "verdict": "error", "error": repr(ex)}) + "\n"); continue
            tally[c["verdict"]] += 1
            rec = {"base": name, "pattern": pat, "key": key, "node6_in": [i for i in range(6) if A[i][6]], "node6_out": [j for j in range(6) if A[6][j]],
                   "adjacency_rows": ["".join(map(str, r)) for r in A], "support": sigma, "verdict": c["verdict"], "stable_cells": c["stable_cells"]}
            if c["verdict"] != "never":
                rec["cert"] = c
            fh.write(json.dumps(rec) + "\n")
            if tally["graphs"] % 100 == 0:
                print(f"  {tally['graphs']}/{len(mine)} {round((time.time()-t0)/tally['graphs'],2)} s/graph {tally}", flush=True); fh.flush()
    tally["seconds"] = round(time.time() - t0, 1)
    json.dump(tally, open(base + "_tally.json", "w"), indent=1); print("DONE", json.dumps(tally), flush=True)

if __name__ == "__main__":
    main()
