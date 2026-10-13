"""
run_region2d_sweep.py -- complete two-parameter certification of the full-support
positivity survivors of a region sweep.

   python3 run_region2d_sweep.py --n 6 --from certificates/region_n6_limit_shard0of3.jsonl ... --shard 0/4

Input: the JSONL records of run_region_sweep.py; every record carrying a "cert"
(i.e. positivity is possible for some legal beta) is re-analysed with
region2d.region2d, which decides the verdict on every cell of a cylindrical
algebraic decomposition of the legal (beta, t) plane. Output per graph: verdict
('never' / 'stable_somewhere' / 'escalate'), the beta-critical values
(isolating interval + irreducible minimal polynomial), the t-degrees of the
irreducible factors of Delta_{m-1}, and for every non-'never' graph the full
certificate. Any stable_somewhere on a non-complete graph goes to *_HITS.jsonl.
"""
import argparse, json, os, sys, time

from enum_cache import key_from_int
from graphs import from_key, is_complete_bidirected
from region2d import region2d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--from", dest="src", nargs="+", required=True)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--out", default="certificates")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    k, K = map(int, args.shard.split("/"))
    n = args.n
    surv = []
    for path in args.src:
        for line in open(path):
            r = json.loads(line)
            if "cert" in r:
                surv.append((r["idx"], r["key"]))
    surv.sort()
    surv = [s for i, s in enumerate(surv) if i % K == k]
    if args.limit:
        surv = surv[:args.limit]
    base = os.path.join(args.out, f"region2d_n{n}_shard{k}of{K}")
    tally = {"graphs": 0, "never": 0, "stable_somewhere": 0, "escalate": 0, "complete": 0, "errors": 0}
    t0 = time.time()
    with open(base + ".jsonl", "w") as fh, open(base + "_HITS.jsonl", "w") as hits:
        for idx, key in surv:
            A = from_key(key_from_int(key, n), n)
            complete = is_complete_bidirected(A)
            tally["graphs"] += 1
            tally["complete"] += complete
            try:
                c = region2d(A, list(range(n)), crosscheck=True)
            except Exception as e:
                tally["errors"] += 1
                fh.write(json.dumps({"key": key, "idx": idx, "verdict": "error", "error": repr(e)}) + "\n")
                print(f"  ERROR idx={idx} key={key}: {repr(e)[:200]}", flush=True)
                continue
            rec = {"key": key, "idx": idx, "complete": complete, "verdict": c["verdict"],
                   "beta_critical": c["beta_critical"], "projection": c["projection_names"],
                   "n_beta_cells": len(c["cells"]), "n_cells": sum(len(r["t_cells"]) for r in c["cells"])}
            tally[c["verdict"]] += 1
            if c["verdict"] != "never":
                rec["cert"] = c
                if not complete:
                    hits.write(json.dumps(rec) + "\n"); hits.flush()
                    print(f"\n*** CANDIDATE COUNTEREXAMPLE (2D) idx={idx} key={key} stable cells={c['stable_cells']} ***", flush=True)
            if complete and c["verdict"] != "stable_somewhere":
                print(f"\n*** COMPLETE GRAPH idx={idx} NOT STABLE -- BUG, STOP ***", flush=True)
                sys.exit(3)
            fh.write(json.dumps(rec) + "\n")
            if tally["graphs"] % 100 == 0:
                el = time.time() - t0
                print(f"  {tally['graphs']}/{len(surv)}  {el/tally['graphs']:.2f}s/graph  {tally}", flush=True)
    tally["seconds"] = round(time.time() - t0, 1)
    json.dump(tally, open(base + "_tally.json", "w"), indent=1)
    print("DONE", json.dumps(tally))


if __name__ == "__main__":
    main()
