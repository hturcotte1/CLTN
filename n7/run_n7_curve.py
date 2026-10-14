"""
run_n7_curve.py -- exact Lemma-R curve stage of the finite-eps n = 7 map (sharded, resumable).

   python3 n7/run_n7_curve.py --from n7/n7_pos_survivors.txt --shard k/K --out n7/curve [--resume]

Input lines "<digraph6> <edges> <P|R|X> <twin>" (domination-free, positive somewhere: the only
graphs that can carry a stable positive full-support fixed point, Lemmas D and A.2).
For each line (shard = line index mod K):
  edges < 22  -> 'never' by the edge bound (legal stability forces B > Nn, i.e. edges >= 22 at
                 m = 7; notes/lemma_M_refuted.md section 5, re-derived in ANOMALIES.md, audited);
  otherwise   -> n7/curve_test.curve_test_fast: 'never' certifies (Lemma R) that the full support
                 is stable at NO legal (eps, delta); 'stable_somewhere' graphs go to region2d
                 (n7/run_n7_region2d.py) for the exact two-parameter region; 'escalate' is recorded.
Every non-'never' graph is written to <out>/curve_hits_<k>of<K>.jsonl with its stable t-cells on
the curve; counts to <out>/curve_stats_<k>of<K>.json; a checkpoint every 5,000 lines allows --resume.
No floating point anywhere.
"""
import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from n7.curve_test import curve_test_fast, from_d6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", required=True)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--out", default="n7/curve")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--no-edge-bound", action="store_true")
    a = ap.parse_args()
    k, K = map(int, a.shard.split("/"))
    os.makedirs(a.out, exist_ok=True)
    tag = f"{k}of{K}"
    hits_p = os.path.join(a.out, f"curve_hits_{tag}.jsonl")
    stats_p = os.path.join(a.out, f"curve_stats_{tag}.json")
    ck_p = os.path.join(a.out, f"curve_ckpt_{tag}.json")
    st = {"shard": a.shard, "lines": 0, "edges_lt_22": 0, "never": 0, "stable_somewhere": 0, "escalate": 0, "errors": 0, "last_line": -1}
    if a.resume and os.path.exists(ck_p):
        st = json.load(open(ck_p)); st["resumed_from_line"] = st["last_line"]
        print(f"resuming after line {st['last_line']}", flush=True)
    mode = "a" if a.resume and os.path.exists(hits_p) else "w"
    t0 = time.time(); n0 = st["lines"]
    with open(hits_p, mode) as hits:
        for li, line in enumerate(open(a.src)):
            if li % K != k or li <= st["last_line"]:
                continue
            parts = line.split()
            if not parts:
                continue
            d6, edges = parts[0], int(parts[1])
            st["lines"] += 1
            if edges < 22 and not a.no_edge_bound:
                st["edges_lt_22"] += 1
            else:
                try:
                    r = curve_test_fast(from_d6(d6))
                    v = r["verdict"]
                except Exception as ex:
                    v = "errors"; r = {"error": repr(ex)}
                st[v] += 1
                if v != "never":
                    rec = {"line": li, "d6": d6, "edges": edges, "flag": parts[2] if len(parts) > 2 else None,
                           "verdict": v, "n_critical": r.get("n_critical"), "deg_Delta": r.get("deg_Delta"),
                           "stable_t_cells": [(c["lo"], c["hi"], c["sample_t"]) for c in r.get("cells", []) if c["stable"]],
                           "undecided_t_cells": [(c["lo"], c["hi"]) for c in r.get("cells", []) if c["stable"] is None],
                           "critical": r.get("critical"), "error": r.get("error")}
                    hits.write(json.dumps(rec) + "\n"); hits.flush()
            st["last_line"] = li
            if st["lines"] % 5000 == 0:
                el = time.time() - t0
                st["seconds"] = round(el, 1); st["ms_per_line"] = round(1000 * el / max(1, st["lines"] - n0), 2)
                json.dump(st, open(ck_p, "w"), indent=1)
                print(f"  {tag} {st['lines']} lines  {st['ms_per_line']} ms/line  never={st['never']} stable={st['stable_somewhere']} esc={st['escalate']} lt22={st['edges_lt_22']}", flush=True)
    st["seconds"] = round(time.time() - t0, 1); st["done"] = True
    json.dump(st, open(ck_p, "w"), indent=1); json.dump(st, open(stats_p, "w"), indent=1)
    print("DONE", json.dumps(st), flush=True)


if __name__ == "__main__":
    main()
