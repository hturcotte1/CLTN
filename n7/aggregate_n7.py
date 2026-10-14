"""
aggregate_n7.py -- consolidate the complete finite-eps n = 7 full-support map.

Inputs: n7/n7_pos_filter_stats.txt (prefilter), n7/*_resolved stats (resolver), n7/curve/curve_TOTAL.json
(curve stage), n7/region2d/region2d_n7_shard*.jsonl (two-parameter stage on the curve-stable graphs),
n7/n7_limit_HITS.jsonl (the 123 limit candidates + K7 via n7_limit_results), nauty labelg for canonical forms.
Outputs: certificates/n7_finite_summary.json (the accounting chain: every one of the 882,033,440 digraphs
is dropped by a proved necessary condition or carries a two-parameter certificate),
certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl (every non-clique graph whose full support is
stable somewhere in the legal plane, with its region2d certificate, its exact limit status and
whether it is one of the 123 limit candidates), and the D-equivalence classes (Lemma J key:
characteristic polynomial of A and the in-walk counts (A^T)^k 1, k = 1..n).
No floating point anywhere.
"""
import sys, os, json, glob, subprocess, re
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from n7.curve_test import from_d6
from fastregion import region, INF, icharpoly
from graphs import is_domination_free, is_complete_bidirected
from ctln import is_clique

LABELG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "labelg")


def canon(d6s):
    out = subprocess.run([LABELG], input="\n".join(d6s) + "\n", capture_output=True, text=True).stdout
    res = [l for l in out.splitlines() if l and not l.startswith(">")]
    assert len(res) == len(d6s), (len(res), len(d6s))
    return res


def d_key(A):
    n = len(A)
    cp = tuple(icharpoly(A))
    v = [1] * n; walks = []
    for _ in range(n):
        v = [sum(A[j][i] * v[j] for j in range(n)) for i in range(n)]   # (A^T v)_i = sum_j A[j][i] v_j : in-walks ending at i
        walks.append(tuple(v))
    per_node = tuple(sorted(tuple(w[i] for w in walks) for i in range(n)))
    return (cp, per_node)


def bidirected_in_twins(A):
    n = len(A)
    return [(u, v) for u in range(n) for v in range(u + 1, n)
            if A[u][v] and A[v][u] and all(A[k][u] == A[k][v] for k in range(n) if k not in (u, v))]


def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(base)
    pre = {}
    for line in open("n7/n7_pos_filter_stats.txt"):
        for k, v in re.findall(r"([a-zA-Z]+) (\d+)", line):
            pre[k] = pre.get(k, 0) + int(v)
    res_stats = {"U": 0, "U_resolved_positive": 0, "U_resolved_never": 0, "other_verdict": 0, "kept": 0, "P": 0}
    for p in glob.glob("n7/n7_pos_part?_resolved.txt.stats.json") + glob.glob(os.environ.get("N7_RESOLVED_STATS_GLOB", "/nonexistent/*")):
        s = json.load(open(p))
        for k in res_stats:
            res_stats[k] += s.get(k, 0)
    curve = json.load(open("n7/curve/curve_TOTAL.json"))
    recs = []
    for p in sorted(glob.glob("n7/region2d/region2d_n7_shard*.jsonl")):
        for line in open(p):
            recs.append(json.loads(line))
    tally = {"graphs": len(recs), "never": 0, "stable_somewhere": 0, "escalate": 0, "error": 0, "complete": 0}
    for r in recs:
        tally[r["verdict"] if r["verdict"] in tally else "error"] += 1
        tally["complete"] += bool(r.get("complete"))
    # known limit candidates by canonical form
    lim = [json.loads(l) for l in open("n7/n7_limit_HITS.jsonl")]
    lim_canon = dict(zip(canon([r["d6"] for r in lim]), [r["d6"] for r in lim]))
    stable = [r for r in recs if r["verdict"] == "stable_somewhere" and not r.get("complete")]
    esc = [r for r in recs if r["verdict"] == "escalate"]
    stable_canon = canon([r["d6"] for r in stable]) if stable else []
    hits = []
    for r, c in zip(stable, stable_canon):
        A = from_d6(r["d6"])
        limv = region(A, list(range(7)), INF, replay=True)
        cells = r["cert"]["stable_cells"]
        hits.append({"d6": r["d6"], "canonical": c, "edges": r["edges"], "clique": bool(is_clique(A, list(range(7)))),
                     "domination_free": bool(is_domination_free(A)), "in_limit_candidates": c in lim_canon,
                     "limit_verdict": limv["verdict"], "limit_intervals": limv.get("stable_intervals"),
                     "stable_cells": cells, "bounded_in_t_somewhere": any(sc[4] is not None for sc in cells),
                     "beta_critical": r.get("beta_critical"), "bidirected_in_twin_pairs": bidirected_in_twins(A),
                     "d_key": None, "region2d": r["cert"]})
    # D-equivalence classes
    keys = {}
    for h in hits:
        k = d_key(from_d6(h["d6"])); keys.setdefault(k, []).append(h["canonical"]); h["d_key"] = [list(k[0]), [list(x) for x in k[1]]]
    classes = sorted(keys.values(), key=lambda c: (-len(c), c))
    from collections import Counter
    sizes = sorted(Counter(len(c) for c in classes).items())
    new = [h for h in hits if not h["in_limit_candidates"]]
    missing_limit = [c for c in lim_canon if c not in set(stable_canon)]
    summary = {
        "n": 7, "scope": "complete finite-eps classification of the FULL support of every digraph on 7 nodes (all legal eps, delta); proper supports reduce to induced subgraphs by Theorem F",
        "status": "CANDIDATES -- awaiting human hand-verification (CLAUDE.md section 3)",
        "chain": {
            "digraphs": pre.get("total"), "A000273_7": 882033440,
            "dominated_dropped_LemmaD": pre.get("total", 0) - pre.get("domfree", 0),
            "domination_free": pre.get("domfree"),
            "never_positive_C_descartes_LemmaA2": pre.get("dropN"),
            "prefilter_kept": pre.get("keptP", 0) + pre.get("keptU", 0), "prefilter_P": pre.get("keptP"), "prefilter_U": pre.get("keptU"),
            "U_resolved_never_exact": res_stats["U_resolved_never"], "U_resolved_positive": res_stats["U_resolved_positive"], "U_other_kept": res_stats["other_verdict"],
            "positivity_survivors": pre.get("keptP", 0) + res_stats["U_resolved_positive"] + res_stats["other_verdict"],
            "curve_stage": curve,
            "region2d_stage": tally,
        },
        "accounting_ok": (pre.get("keptP", 0) + res_stats["U_resolved_positive"] + res_stats["other_verdict"] == curve["lines"]
                          and curve["lines"] == curve["edges_lt_22"] + curve["never"] + curve["stable_somewhere"] + curve["escalate"] + curve["errors"]
                          and curve["stable_somewhere"] + curve["escalate"] + curve["errors"] == tally["graphs"]),
        "stable_noncomplete": len(hits), "complete_graph_stable": tally["complete"] == 1 and any(r.get("complete") and r["verdict"] == "stable_somewhere" for r in recs),
        "escalations": len(esc), "escalated": [r["d6"] for r in esc],
        "limit_candidates": len(lim), "limit_candidates_all_stable_here": not missing_limit, "limit_candidates_missing": missing_limit,
        "new_finite_eps_only": len(new), "new_finite_eps_only_d6": [h["d6"] for h in new],
        "bounded_in_t_count": sum(h["bounded_in_t_somewhere"] for h in hits),
        "D_equivalence_classes": len(classes), "class_sizes": sizes, "classes": classes,
        "edge_counts": sorted(Counter(h["edges"] for h in hits).items()),
        "with_bidirected_in_twins": sum(bool(h["bidirected_in_twin_pairs"]) for h in hits),
    }
    json.dump(summary, open("certificates/n7_finite_summary.json", "w"), indent=1)
    with open("certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl", "w") as fh:
        for h in hits:
            fh.write(json.dumps(h) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("classes", "region2d")}, default=str)[:3000])


if __name__ == "__main__":
    main()
