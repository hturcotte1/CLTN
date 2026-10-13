"""
run_region_sweep.py -- Stage 2/3 driver: exact q-regions for the FULL SUPPORT of
every digraph on n nodes, at a fixed eps-slice or in the eps -> 0 limit.

   python3 run_region_sweep.py --n 6 --mode limit --shard 0/4
   python3 run_region_sweep.py --n 6 --mode slice --eps 1/100 --shard 0/4

Reads certificates/digraphs_n<n>.keys (validated against A000273). For every
graph: fastregion.beta_data + region(). Graphs whose positivity set is empty
(no legal beta with a common sign for all a_i) get a one-line 'never' record;
everything else gets the full certificate. Any 'stable_somewhere' verdict on a
NON-complete graph is a candidate counterexample: it is written to the
*_HITS.jsonl file and the run keeps going (the orchestrator stops and
hand-verifies; see CLAUDE.md §3), unless --stop-on-hit is given.

No model in the loop. Deterministic. Sharded by index modulo.
"""
import argparse, json, os, sys, time
from fractions import Fraction as F

from enum_cache import load_keys, key_from_int
from graphs import from_key, is_complete_bidirected, is_domination_free
from fastregion import beta_data, region, INF, stable_q_intervals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--mode", choices=["limit", "slice"], required=True)
    ap.add_argument("--eps", type=F, default=None)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--limit", type=int, default=None, help="process at most this many graphs of the shard")
    ap.add_argument("--out", default="certificates")
    ap.add_argument("--crosscheck-every", type=int, default=50,
                    help="run ctln.analyze_support at every sample point for 1 in N graphs (slice mode)")
    ap.add_argument("--stop-on-hit", action="store_true")
    args = ap.parse_args()
    k, K = map(int, args.shard.split("/"))
    if args.mode == "slice":
        assert args.eps is not None and 0 < args.eps < 1
        t = (1 - args.eps) / args.eps
        tag = f"eps{args.eps.numerator}-{args.eps.denominator}"
    else:
        t = INF
        tag = "limit"
    n = args.n
    keys = load_keys(n)
    os.makedirs(args.out, exist_ok=True)
    base = os.path.join(args.out, f"region_n{n}_{tag}_shard{k}of{K}")
    tally = {"graphs": 0, "complete": 0, "never_positivity": 0, "never_full": 0,
             "stable_somewhere": 0, "escalate": 0, "anomaly": 0, "crosschecked": 0,
             "domfree": 0, "domfree_positivity_possible": 0, "dominated_positivity_possible": 0}
    t0 = time.time()
    with open(base + ".jsonl", "w") as fh, open(base + "_HITS.jsonl", "w") as hits:
        for idx in range(k, len(keys), K):
            if args.limit and tally["graphs"] >= args.limit:
                break
            v = keys[idx]
            A = from_key(key_from_int(v, n), n)
            tally["graphs"] += 1
            complete = is_complete_bidirected(A)
            domfree = is_domination_free(A)
            tally["complete"] += complete
            tally["domfree"] += domfree
            bd = beta_data(A, list(range(n)))
            cc = (args.mode == "slice") and (tally["graphs"] % args.crosscheck_every == 0)
            tally["crosschecked"] += cc
            c = region(A, list(range(n)), t, bd=bd, crosscheck=cc)
            rec = {"key": v, "idx": idx, "complete": complete, "domfree": domfree, "verdict": c["verdict"]}
            if c["verdict"] == "never" and "cells" not in c:
                tally["never_positivity"] += 1
                rec["reason"] = c.get("reason")
                if c.get("anomaly"):
                    tally["anomaly"] += 1
            else:
                if domfree:
                    tally["domfree_positivity_possible"] += 1
                else:
                    tally["dominated_positivity_possible"] += 1
                rec["cert"] = c
                if c["verdict"] == "never":
                    tally["never_full"] += 1
                elif c["verdict"] == "stable_somewhere":
                    tally["stable_somewhere"] += 1
                    rec["stable_q"] = stable_q_intervals(c)
                    if not complete:
                        hits.write(json.dumps(rec) + "\n"); hits.flush()
                        print(f"\n*** CANDIDATE COUNTEREXAMPLE idx={idx} key={v} stable q={rec['stable_q']} ***", flush=True)
                        if args.stop_on_hit:
                            fh.write(json.dumps(rec) + "\n")
                            sys.exit(2)
                else:
                    tally["escalate"] += 1
                if c.get("anomaly"):
                    tally["anomaly"] += 1
            if complete and c["verdict"] != "stable_somewhere":
                print(f"\n*** COMPLETE GRAPH NOT STABLE EVERYWHERE idx={idx}: {c['verdict']} -- BUG, STOP ***", flush=True)
                sys.exit(3)
            fh.write(json.dumps(rec) + "\n")
            if tally["graphs"] % 2000 == 0:
                el = time.time() - t0
                print(f"  {tally['graphs']} graphs  {tally['graphs']/el:.1f}/s  {tally}", flush=True)
    tally["seconds"] = round(time.time() - t0, 1)
    json.dump(tally, open(base + "_tally.json", "w"), indent=1)
    print("DONE", json.dumps(tally))


if __name__ == "__main__":
    main()
