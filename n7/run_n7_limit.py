"""
run_n7_limit.py -- Stage 6a: exact eps -> 0 (limit) classification of the full
support for all digraphs on n = 7 nodes.

Input: n7/n7_survivors_part*.txt from filter7 (geng | directg | filter7):
every digraph on 7 nodes that is domination-free (Lemma D) and whose compressed
adjacency has spectral abscissa < -1/2 ('S') or possibly == -1/2 ('B').
By Lemma I.3 a graph is limit-stable at some legal beta iff alpha(A|perp) < -1/2
strictly, so the 'B' graphs are certified never limit-stable by the filter's
exact Hurwitz test; only the 'S' graphs are analysed here, with
fastregion.region(A, V, INF, replay=True): exact positivity cells in beta and
the exact beta-interval (2, beta_c) of limit stability, cross-checked.

Output: n7/n7_limit_results.jsonl (one record per 'S' graph), n7/n7_limit_HITS.jsonl
(limit-stable non-complete graphs, each with pointwise exact witnesses from
ctln.analyze_support at large t), n7/n7_limit_tally.json.
"""
import glob, json, sys, time, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fastregion import region, beta_data, INF, stable_q_intervals
from ctln import analyze_support, legal, is_clique, targets
from graphs import canonical_form, is_complete_bidirected, is_domination_free
from enum_cache import key_to_int


def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


def main():
    n = 7
    rows = []
    for path in sorted(glob.glob("n7/n7_survivors_part*.txt")):
        for line in open(path):
            s, e, f = line.split()
            rows.append((s, int(e), f))
    nS = sum(1 for r in rows if r[2] == "S"); nB = len(rows) - nS
    print(f"{len(rows)} survivors: {nS} strict (alpha < -1/2), {nB} boundary (never limit-stable)", flush=True)
    t0 = time.time()
    tally = {"strict": nS, "boundary": nB, "never_positivity": 0, "never_full": 0, "stable_somewhere": 0,
             "escalate": 0, "complete": 0, "stable_noncomplete": 0}
    with open("n7/n7_limit_results.jsonl", "w") as fh, open("n7/n7_limit_HITS.jsonl", "w") as hits:
        for i, (s, e, f) in enumerate(r for r in rows if r[2] == "S"):
            A = from_d6(s)
            assert is_domination_free(A)
            key = key_to_int(canonical_form(A)[0])
            complete = is_complete_bidirected(A)
            bd = beta_data(A, list(range(n)))
            c = region(A, list(range(n)), INF, bd=bd, replay=True)
            rec = {"d6": s, "key": key, "edges": e, "complete": complete, "verdict": c["verdict"]}
            tally["complete"] += complete
            if c["verdict"] == "never" and "cells" not in c:
                tally["never_positivity"] += 1; rec["reason"] = c.get("reason")
            else:
                rec["cert"] = c
                if c["verdict"] == "never":
                    tally["never_full"] += 1
                elif c["verdict"] == "stable_somewhere":
                    tally["stable_somewhere"] += 1
                    rec["stable_q_limit"] = stable_q_intervals(c)
                    if not complete:
                        tally["stable_noncomplete"] += 1
                        # pointwise exact witnesses at a rational beta inside the first stable cell, large t
                        stable_cells = [cell for cell in c["cells"] if cell["positive"] and cell["off_ok"] and cell["stable"] is True]
                        beta0 = F(stable_cells[0]["sample_beta"])
                        wit = []
                        for t in (F(100), F(1000), F(10000), F(100000)):
                            eps = 1 / (t + 1); delta = (beta0 - 1) * eps
                            r = analyze_support(A, list(range(n)), eps, delta)
                            wit.append({"beta": str(beta0), "t": str(t), "eps": str(eps), "delta": str(delta), "q": r["q"], "legal": legal(eps, delta),
                                        "is_fixed_point": r["is_fixed_point"], "stable": r["stable"], "stability_detail": r["stability_detail"],
                                        "agrees_with_conjecture": r["agrees"], "x": r["x"]})
                        rec["A"] = A; rec["clique"] = is_clique(A, range(n)); rec["targets"] = targets(A, range(n))
                        rec["pointwise_witnesses"] = wit
                        hits.write(json.dumps(rec) + "\n"); hits.flush()
                        print(f"\n*** LIMIT COUNTEREXAMPLE at n=7: {s} edges={e} key={key} stable q={rec['stable_q_limit']}; witnesses stable: {[w['stable'] for w in wit]} ***", flush=True)
                else:
                    tally["escalate"] += 1
            fh.write(json.dumps(rec) + "\n")
            if (i + 1) % 100 == 0:
                print(f"  {i+1}/{nS}  {time.time()-t0:.0f}s  {tally}", flush=True)
    tally["seconds"] = round(time.time() - t0, 1)
    json.dump(tally, open("n7/n7_limit_tally.json", "w"), indent=1)
    print("DONE", json.dumps(tally))


if __name__ == "__main__":
    main()
