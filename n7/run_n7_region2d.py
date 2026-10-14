"""
run_n7_region2d.py -- complete two-parameter (beta, t) certification of the
full support for the n = 7 positivity survivors (output of filter7p + resolve_U).

   python3 n7/run_n7_region2d.py --from n7/n7_pos_survivors.txt --shard 3/64 --out n7/region2d [--resume]

Input lines: "<digraph6> <edges> <P|R> <twin>" (every graph is domination-free and
has some beta > 2 with all a_i(beta) > 0; nothing else can carry a stable
positive full-support fixed point at any legal parameter, Lemmas D and A.2).
Each graph is analysed with region2d.region2d (cylindrical algebraic
decomposition of the legal (beta, t) plane on the proven boundary set
{a_i, D, Delta_{m-1}, L}, every cell decided exactly by Hurwitz determinants,
cross-checked pointwise, root counts replayed). Output per graph: verdict
'never' / 'stable_somewhere' / 'escalate', beta-critical values (isolating
interval + minimal polynomial), cell counts, and for every non-'never' graph
the full certificate; non-complete stable graphs also go to *_HITS.jsonl.
Sharding is by line index modulo K; --resume keeps the valid records of an
interrupted run. Exits 3 if the complete graph K7 is not stable (engine bug).
No floating point anywhere.
"""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from region2d import region2d
from graphs import is_complete_bidirected


def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", nargs="+", required=True)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--out", default="n7/region2d")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--no-crosscheck", action="store_true", help="skip the pointwise ctln cross-check (faster; replay still on)")
    args = ap.parse_args()
    k, K = map(int, args.shard.split("/"))
    rows = []
    li = 0
    for path in args.src:
        for line in open(path):
            parts = line.split()
            if not parts:
                continue
            if li % K == k:
                rows.append((li, parts[0], int(parts[1]), parts[2], parts[3] if len(parts) > 3 else "0"))
            li += 1
    if args.limit:
        rows = rows[:args.limit]
    os.makedirs(args.out, exist_ok=True)
    base = os.path.join(args.out, f"region2d_n7_shard{k}of{K}")
    tally = {"graphs": 0, "never": 0, "stable_somewhere": 0, "escalate": 0, "complete": 0, "errors": 0, "shard": args.shard, "assigned": len(rows)}
    t0 = time.time()
    done = set(); mode = "w"
    if args.resume and os.path.exists(base + ".jsonl"):
        good = []
        for line in open(base + ".jsonl"):
            try:
                r = json.loads(line)
            except Exception:
                break
            good.append(line if line.endswith("\n") else line + "\n")
            done.add(r["line"])
            tally["graphs"] += 1
            tally["complete"] += bool(r.get("complete"))
            tally[r["verdict"] if r["verdict"] in tally else "errors"] += 1
        with open(base + ".jsonl", "w") as fh:
            fh.writelines(good)
        mode = "a"; tally["resumed_records"] = len(good)
        print(f"resuming: {len(good)} records kept", flush=True)
    with open(base + ".jsonl", mode) as fh, open(base + "_HITS.jsonl", mode) as hits:
        for line_no, s, e, f, tw in rows:
            if line_no in done:
                continue
            A = from_d6(s)
            complete = is_complete_bidirected(A)
            tally["graphs"] += 1; tally["complete"] += complete
            try:
                c = region2d(A, list(range(len(A))), crosscheck=not args.no_crosscheck, replay=True)
            except Exception as ex:
                tally["errors"] += 1
                fh.write(json.dumps({"line": line_no, "d6": s, "verdict": "error", "error": repr(ex)}) + "\n"); fh.flush()
                print(f"  ERROR line={line_no} d6={s}: {repr(ex)[:200]}", flush=True)
                continue
            rec = {"line": line_no, "d6": s, "edges": e, "flag": f, "twin": tw, "complete": complete,
                   "verdict": c["verdict"], "beta_critical": c["beta_critical"], "projection": c["projection_names"],
                   "n_beta_cells": len(c["cells"]), "n_cells": sum(len(r["t_cells"]) for r in c["cells"])}
            tally[c["verdict"]] += 1
            if c["verdict"] != "never":
                rec["cert"] = c
                if not complete:
                    hits.write(json.dumps(rec) + "\n"); hits.flush()
                    print(f"\n*** CANDIDATE COUNTEREXAMPLE (2D, finite eps) line={line_no} d6={s} stable cells={c['stable_cells']} ***", flush=True)
            if complete and c["verdict"] != "stable_somewhere":
                print(f"\n*** COMPLETE GRAPH line={line_no} NOT STABLE -- BUG, STOP ***", flush=True)
                sys.exit(3)
            fh.write(json.dumps(rec) + "\n")
            if tally["graphs"] % 50 == 0:
                el = time.time() - t0
                print(f"  {tally['graphs']}/{len(rows)}  {el/max(1,tally['graphs']-tally.get('resumed_records',0)):.2f}s/graph  {tally}", flush=True)
                json.dump(dict(tally, seconds=round(el, 1), partial=True), open(base + "_tally.json", "w"), indent=1)
    tally["seconds"] = round(time.time() - t0, 1)
    json.dump(tally, open(base + "_tally.json", "w"), indent=1)
    print("DONE", json.dumps(tally))


if __name__ == "__main__":
    main()
