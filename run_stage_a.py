"""
run_stage_a.py -- Stage A cheap exact filter, full-support case at n=6.

Usage:  python3 run_stage_a.py --n 6 --eps 1/2 --q 11/10 3/2 2 --limit 5000

Writes one JSONL certificate per case to certificates/stageA_<n>_<eps>.jsonl.
Deterministic, no model in the loop. Safe to shard by --offset/--limit and run
many copies in parallel.

This is a FILTER, not an answer. Its job is to shrink the set that goes to
symbolic.q_region in Stage B. See CLAUDE.md section 5.
"""
import argparse, json, os, sys, time
from fractions import Fraction as F

from ctln import analyze_support, legal, q_of
from graphs import enumerate_digraphs, from_key, is_mixed, is_complete_bidirected


def frac(s):
    return F(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--eps", type=frac, default=F(1, 2))
    ap.add_argument("--q", type=frac, nargs="+", default=[F(11, 10), F(3, 2), F(2)])
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default="certificates")
    ap.add_argument("--mixed-only", action="store_true", default=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    print(f"enumerating digraphs on {args.n} nodes up to isomorphism ...", flush=True)
    t0 = time.time()
    keys = sorted(enumerate_digraphs(args.n))
    print(f"  {len(keys)} canonical graphs  ({time.time()-t0:.1f}s)", flush=True)

    sel = keys[args.offset: args.offset + args.limit if args.limit else None]
    path = os.path.join(args.out,
                        f"stageA_n{args.n}_eps{args.eps.numerator}-{args.eps.denominator}"
                        f"_off{args.offset}.jsonl")

    tally = {"decided_stable": 0, "decided_unstable": 0,
             "undecided": 0, "disagree": 0, "skipped": 0}
    t0 = time.time()
    with open(path, "w") as fh:
        for idx, key in enumerate(sel):
            A = from_key(key, args.n)
            if args.mixed_only and not is_mixed(A) and not is_complete_bidirected(A):
                tally["skipped"] += 1
                continue
            for qv in args.q:
                delta = qv * args.eps / (1 - args.eps)
                if not legal(args.eps, delta):
                    continue
                rec = analyze_support(A, list(range(args.n)), args.eps, delta)
                rec["graph_key"] = list(key)
                rec["n"] = args.n
                rec["checker"] = "ctln.py/stageA-1"
                if rec["BOUNDARY"]:
                    tally["undecided"] += 1
                elif rec["agrees"] is False:
                    tally["disagree"] += 1
                    print("\n*** DISAGREEMENT -- STOPPING ***", flush=True)
                    print(json.dumps(rec), flush=True)
                    fh.write(json.dumps(rec) + "\n")
                    print("Hand-verify before continuing. See CLAUDE.md section 3.")
                    sys.exit(2)
                elif rec["stable_fixed_point"]:
                    tally["decided_stable"] += 1
                else:
                    tally["decided_unstable"] += 1
                fh.write(json.dumps(rec) + "\n")
            if idx % 2000 == 0 and idx:
                el = time.time() - t0
                print(f"  {idx}/{len(sel)}  {idx/el:.0f} graphs/s  {tally}", flush=True)

    print(f"\ndone in {time.time()-t0:.0f}s -> {path}")
    print(json.dumps(tally, indent=2))
    print("\nStage B input = every record with BOUNDARY true, plus any graph whose")
    print("verdict changes between q values. See CLAUDE.md section 5.")


if __name__ == "__main__":
    main()
