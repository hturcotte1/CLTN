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
    ap = argparse.ArgumentParser(); ap.add_argument("--shard", default="0/1"); ap.add_argument("--out", default="certificates"); ap.add_argument("--resume", action="store_true"); ap.add_argument("--resume-from", default=None)
    a = ap.parse_args(); k, K = map(int, a.shard.split("/"))
    hits = load_hits(); sigma = [0, 1, 2, 3, 4, 5]
    # Pairs (graph, support) up to isomorphism: the three bases are pairwise non-isomorphic, so two pairs are
    # isomorphic iff they come from the same base and their node-6 patterns are related by an automorphism of
    # the base. (Deduplicating by the 7-node graph alone would merge pairs with DIFFERENT supports.)
    import itertools
    cases = []
    for name, H in hits.items():
        auts = [perm for perm in itertools.permutations(range(6)) if all(H[perm[i]][perm[j]] == H[i][j] for i in range(6) for j in range(6))]
        seen = set()
        for pat in range(4096):
            ins = [(pat >> i) & 1 for i in range(6)]; outs = [(pat >> (6 + j)) & 1 for j in range(6)]
            canon = min(tuple(ins[perm[i]] for i in range(6)) + tuple(outs[perm[j]] for j in range(6)) for perm in auts)
            if canon in seen:
                continue
            seen.add(canon)
            A = [row[:] + [ins[i]] for i, row in enumerate(H)] + [outs + [0]]
            key = key_to_int(canonical_form(A)[0])
            cases.append((name, pat, key, A))
        print(f"{name}: |Aut| = {len(auts)}, {len(seen)} pattern classes", flush=True)
    cases.sort(key=lambda c: (c[0], c[1]))
    mine = [c for i, c in enumerate(cases) if i % K == k]
    print(f"{len(cases)} distinct graphs; shard {a.shard} has {len(mine)}", flush=True)
    base = os.path.join(a.out, f"n7_proper_supports_shard{k}of{K}")
    tally = {"graphs": 0, "never": 0, "stable_somewhere": 0, "escalate": 0, "errors": 0, "distinct_total": len(cases)}
    done = {}
    if a.resume:
        import glob as _g
        for pth in _g.glob(os.path.join(a.out, "n7_proper_supports_shard*of*.jsonl")) + _g.glob(a.resume_from or "/nonexistent"):
            for line in open(pth):
                r = json.loads(line); done[(r["base"], r["pattern"])] = line
    t0 = time.time()
    with open(base + ".jsonl", "w") as fh:
        for name, pat, key, A in mine:
            if (name, pat) in done:
                r = json.loads(done[(name, pat)]); tally["graphs"] += 1; tally[r["verdict"] if r["verdict"] in tally else "errors"] += 1
                fh.write(done[(name, pat)] if done[(name, pat)].endswith("\n") else done[(name, pat)] + "\n"); continue
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
