"""
run_n8_limit.py -- exact eps -> 0 (limit) classification of the full support for the n = 8 prefilter
survivors (S = alpha < -1/2 strictly, O = overflow-undecided) produced by n8/run_n8_local.sh.
Same logic as n7/run_n7_limit.py: fastregion.region(A, V, INF, replay=True) gives the exact positivity
cells and the exact limit-stable beta-interval; every non-complete limit-stable graph is a CANDIDATE
(hand-verification required) and gets exact pointwise witnesses from ctln.analyze_support at large t.
   python3 n8/run_n8_limit.py [--parts n8/parts] [--out n8]
Resumable: graphs already in the results file are skipped. No floating point anywhere.
"""
import glob, json, sys, time, os, argparse
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fastregion import region, beta_data, INF
from ctln import analyze_support, is_clique
from graphs import is_domination_free, is_complete_bidirected
from n7.curve_test import from_d6


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--parts", default="n8/parts"); ap.add_argument("--out", default="n8")
    a = ap.parse_args(); n = 8
    rows = []
    for path in sorted(glob.glob(os.path.join(a.parts, "part*.txt"))):
        for line in open(path):
            p = line.split()
            if len(p) >= 3 and p[2] in ("S", "O"):
                rows.append((p[0], int(p[1]), p[2]))
    res_p = os.path.join(a.out, "n8_limit_results.jsonl"); hits_p = os.path.join(a.out, "n8_limit_HITS.jsonl")
    done = set()
    if os.path.exists(res_p):
        for line in open(res_p):
            try: done.add(json.loads(line)["d6"])
            except Exception: pass
    tally = {"S": sum(r[2] == "S" for r in rows), "O": sum(r[2] == "O" for r in rows), "never_positivity": 0, "never_full": 0, "stable_somewhere": 0, "escalate": 0, "complete": 0, "stable_noncomplete": 0, "skipped_done": len(done)}
    t0 = time.time()
    with open(res_p, "a") as fh, open(hits_p, "a") as hits:
        for s, e, f in rows:
            if s in done: continue
            A = from_d6(s)
            if not is_domination_free(A): raise RuntimeError("prefilter invariant violated: " + s)
            complete = is_complete_bidirected(A)
            c = region(A, list(range(n)), INF, bd=beta_data(A, list(range(n))), replay=True)
            rec = {"d6": s, "edges": e, "flag": f, "complete": complete, "verdict": c["verdict"]}
            tally["complete"] += complete
            if c["verdict"] == "never" and "cells" not in c:
                tally["never_positivity"] += 1; rec["reason"] = c.get("reason")
            else:
                rec["cert"] = c
                if c["verdict"] == "never": tally["never_full"] += 1
                elif c["verdict"] == "stable_somewhere":
                    tally["stable_somewhere"] += 1
                    if not complete:
                        tally["stable_noncomplete"] += 1
                        iv = c["stable_intervals"][0]; lo = iv["lo"]; hi = iv["hi"]
                        blo = F(lo["beta"]) if "beta" in lo else F(lo["beta_interval"][1]); bhi = F(hi["beta"]) if "beta" in hi else F(hi["beta_interval"][0])
                        beta = (blo + bhi) / 2 if bhi > blo else blo + F(1, 1000)
                        rec["witnesses"] = []
                        for t in (100, 1000, 10000, 100000):
                            eps = F(1, t + 1); delta = (beta - 1) * eps
                            w = analyze_support(A, list(range(n)), eps, delta)
                            rec["witnesses"].append({"beta": str(beta), "t": t, "eps": str(eps), "delta": str(delta), "q": str(delta * (1 - eps) / eps), "stable": w["stable"], "positive": bool(w["x"] and all(F(str(v)) > 0 for v in w["x"])), "detail": w.get("stability_detail")})
                        hits.write(json.dumps(rec) + "\n"); hits.flush()
                        print(f"*** n=8 LIMIT CANDIDATE {s} edges={e} ***", flush=True)
                else:
                    tally["escalate"] += 1
            fh.write(json.dumps(rec) + "\n")
            if (tally["never_positivity"] + tally["never_full"] + tally["stable_somewhere"] + tally["escalate"]) % 500 == 0: fh.flush()
    tally["seconds"] = round(time.time() - t0, 1)
    json.dump(tally, open(os.path.join(a.out, "n8_limit_tally.json"), "w"), indent=1); print("DONE", json.dumps(tally))


if __name__ == "__main__":
    main()
