"""
aggregate.py -- Stage 5: aggregate the n-node full-support certification into a map.

Reads  certificates/region_n<n>_limit_shard*of*.jsonl        (limit sweep, all graphs)
       certificates/region2d_n<n>_shard*of*.jsonl            (2-D CAD on positivity survivors)
       certificates/region_n<n>_eps*_shard*_only.jsonl        (eps-slices on survivors, cross-checks)
Writes certificates/n<n>_summary.json and notes/n<n>_map.md

For every graph that is stable somewhere (2-D verdict), lists: adjacency, edge
count, complement edges, automorphism group order, isomorphism to Geneson's
graph, limit data (q_c as minimal polynomial + interval, compressed-adjacency
characteristic polynomial), the 2-D stable cells, and which slices saw it.
"""
import argparse, glob, json, os, re
from fractions import Fraction as F
from itertools import permutations

from enum_cache import key_from_int
from graphs import from_key, is_complete_bidirected, is_domination_free, canonical_form
from fastregion import icharpoly, beta_data, region, INF, stable_q_intervals

GENESON = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]


def compressed_charpoly(A):
    """m * det(lam I - A|perp) = chi_A(lam) - chi_{A+J}(lam), integer coefficients (hi->lo)."""
    m = len(A)
    cA = icharpoly(A)
    cAJ = icharpoly([[A[i][j] + 1 for j in range(m)] for i in range(m)])
    diff = [a - b for a, b in zip(cA, cAJ)]
    while diff and diff[0] == 0:
        diff = diff[1:]
    return diff


def automorphisms(A):
    n = len(A)
    return sum(1 for p in permutations(range(n)) if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)))


def features(A):
    n = len(A)
    E = sum(A[i][j] for i in range(n) for j in range(n))
    nonedges = [(i, j) for i in range(n) for j in range(n) if i != j and A[i][j] == 0]
    mutual_nonadj = [(i, j) for i in range(n) for j in range(i + 1, n) if A[i][j] == 0 and A[j][i] == 0]
    bidir = [(i, j) for i in range(n) for j in range(i + 1, n) if A[i][j] == 1 and A[j][i] == 1]
    indeg = [sum(A[j][i] for j in range(n)) for i in range(n)]
    outdeg = [sum(A[i][j] for j in range(n)) for i in range(n)]
    return {"edges": E, "nonedges": nonedges, "mutually_nonadjacent_pairs": mutual_nonadj,
            "bidirected_pairs": len(bidir), "oneway_edges": E - 2 * len(bidir),
            "indeg": indeg, "outdeg": outdeg, "automorphisms": automorphisms(A),
            "domination_free": is_domination_free(A),
            "isomorphic_to_geneson": canonical_form(A)[0] == canonical_form(GENESON)[0] if n == 6 else False,
            "compressed_charpoly_hi_to_lo": compressed_charpoly(A)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    args = ap.parse_args()
    n = args.n
    cdir = "certificates"
    # The canonical limit run is the single-shard replay with the final code
    # (region_n<n>_limit_shard0of1.jsonl) when present; the earlier three-shard
    # run mixed two positivity filters (ANOMALIES 2026-10-13) and is used only
    # as a fallback. The 2-D certification covered the SUPERSET of survivors of
    # both runs, which is reported separately as 'certified_superset'.
    limit = {}
    canonical = os.path.join(cdir, f"region_n{n}_limit_shard0of1.jsonl")
    paths = [canonical] if os.path.exists(canonical) else sorted(
        p for p in glob.glob(os.path.join(cdir, f"region_n{n}_limit_shard*of*.jsonl")) if "_only" not in p)
    for path in paths:
        for line in open(path):
            try:
                r = json.loads(line)
            except Exception:
                continue
            limit[r["idx"]] = r
    two_d = {}
    for path in sorted(glob.glob(os.path.join(cdir, f"region2d_n{n}_shard*of*.jsonl"))):
        for line in open(path):
            try:
                r = json.loads(line)
            except Exception:
                continue
            two_d[r["idx"]] = r
    slices = {}
    for path in sorted(glob.glob(os.path.join(cdir, f"region_n{n}_eps*_shard*_only.jsonl"))):
        tag = re.search(r"_eps([0-9]+-[0-9]+)_", path).group(1)
        for line in open(path):
            try:
                r = json.loads(line)
            except Exception:
                continue
            slices.setdefault(tag, {})[r["idx"]] = r["verdict"]

    def count(d, f):
        return sum(1 for r in d.values() if f(r))

    summary = {
        "n": n,
        "limit": {"graphs": len(limit),
                  "never_positivity": count(limit, lambda r: r["verdict"] == "never" and "cert" not in r),
                  "positivity_possible": count(limit, lambda r: "cert" in r),
                  "limit_stable": count(limit, lambda r: r["verdict"] == "stable_somewhere"),
                  "limit_stable_noncomplete": count(limit, lambda r: r["verdict"] == "stable_somewhere" and not r["complete"]),
                  "escalate": count(limit, lambda r: r["verdict"] == "escalate"),
                  "domination_free": count(limit, lambda r: r["domfree"]),
                  "dominated_positivity_possible": count(limit, lambda r: "cert" in r and not r["domfree"])},
        "two_d": {"graphs": len(two_d),
                  "never": count(two_d, lambda r: r["verdict"] == "never"),
                  "stable_somewhere": count(two_d, lambda r: r["verdict"] == "stable_somewhere"),
                  "stable_noncomplete": count(two_d, lambda r: r["verdict"] == "stable_somewhere" and not r["complete"]),
                  "escalate": count(two_d, lambda r: r["verdict"] == "escalate"),
                  "errors": count(two_d, lambda r: r["verdict"] == "error")},
        "slices": {tag: {"graphs": len(d), "stable_somewhere": sum(1 for v in d.values() if v == "stable_somewhere"),
                         "escalate": sum(1 for v in d.values() if v == "escalate")} for tag, d in slices.items()},
        "consistency": {},
        "hits": [],
    }
    # consistency: every 2-D 'never' must be 'never' in the limit and in every slice; every limit/slice stable must be 2-D stable
    bad = []
    for idx, r in two_d.items():
        lim = limit.get(idx, {}).get("verdict")
        if r["verdict"] == "never" and lim == "stable_somewhere":
            bad.append(("limit stable but 2D never", idx))
        if r["verdict"] == "stable_somewhere" and lim not in ("stable_somewhere", "never", None):
            pass
        for tag, d in slices.items():
            if idx in d and d[idx] == "stable_somewhere" and r["verdict"] != "stable_somewhere":
                bad.append((f"slice {tag} stable but 2D {r['verdict']}", idx))
    summary["consistency"] = {"violations": bad, "survivors_without_2d": sorted(set(i for i, r in limit.items() if "cert" in r) - set(two_d))[:20],
                              "n_survivors_without_2d": len(set(i for i, r in limit.items() if "cert" in r) - set(two_d)),
                              "limit_run_used": [os.path.basename(p) for p in paths],
                              "certified_superset_2d": len(two_d)}
    summary["limit"]["note"] = ("never_positivity = no legal beta with all a_i > 0 (and all o_k < 0); by Lemma A.2 such a graph has no "
                                "STABLE positive fixed point -- it may still have a positive but unstable one")

    # a stable record must have no undecided cell (verdict precedence favours 'stable_somewhere')
    undecided_in_stable = [i for i, r in two_d.items() if r["verdict"] == "stable_somewhere" and
                           any(c["stable_fp"] is None for row in r.get("cert", {}).get("cells", []) for c in row["t_cells"])]
    summary["consistency"]["stable_records_with_undecided_cells"] = undecided_in_stable
    hit_idx = sorted(set(i for i, r in two_d.items() if r["verdict"] == "stable_somewhere" and not r["complete"]) |
                     set(i for i, r in limit.items() if r["verdict"] == "stable_somewhere" and not r["complete"]))
    for idx in hit_idx:
        key = (two_d.get(idx) or limit[idx])["key"]
        A = from_key(key_from_int(key, n), n)
        h = {"idx": idx, "key": key, "A": A, "features": features(A)}
        lr = limit.get(idx)
        if lr and "cert" in lr:
            h["limit_verdict"] = lr["verdict"]
            h["limit_stable_q"] = lr.get("stable_q")
            c = lr["cert"]
            h["limit_stable_intervals"] = c.get("stable_intervals")
        tr = two_d.get(idx)
        if tr:
            h["two_d_verdict"] = tr["verdict"]
            h["two_d_stable_cells"] = tr.get("cert", {}).get("stable_cells")
            h["two_d_beta_critical"] = tr.get("beta_critical")
        h["slices"] = {tag: d.get(idx) for tag, d in slices.items()}
        summary["hits"].append(h)
    os.makedirs("notes", exist_ok=True)
    json.dump(summary, open(os.path.join(cdir, f"n{n}_summary.json"), "w"), indent=1)

    L = [f"# n = {n} full-support map (auto-generated by aggregate.py)\n"]
    L.append(f"Limit sweep: {summary['limit']}\n")
    L.append(f"2-D certification: {summary['two_d']}\n")
    L.append(f"Slices: {summary['slices']}\n")
    L.append(f"Consistency: {summary['consistency']}\n")
    L.append(f"\n## Stable non-complete full supports: {len(summary['hits'])}\n")
    for h in summary["hits"]:
        f = h["features"]
        L.append(f"\n### idx {h['idx']} key {h['key']}  edges {f['edges']}  bidirected pairs {f['bidirected_pairs']}  "
                 f"mutually non-adjacent pairs {f['mutually_nonadjacent_pairs']}  |Aut| {f['automorphisms']}  "
                 f"domination-free {f['domination_free']}  isomorphic to Geneson: {f['isomorphic_to_geneson']}\n")
        L.append("```\n" + "\n".join("".join(str(v) for v in row) for row in h["A"]) + "\n```\n")
        L.append(f"- compressed adjacency charpoly (m*det(lam I - A|perp), hi->lo): {f['compressed_charpoly_hi_to_lo']}\n")
        L.append(f"- limit verdict: {h.get('limit_verdict')}  stable q (limit): {h.get('limit_stable_q')}\n")
        L.append(f"- 2-D verdict: {h.get('two_d_verdict')}  stable cells, each as (beta-cell lo, beta-cell hi, sample beta, t-range at the sample beta: lo, hi) -- "
                 f"the t-range is the slice at the sample beta (t_lo = 1/(sample beta - 2) means 'every legal t'), NOT a box: {h.get('two_d_stable_cells')}\n")
        L.append(f"- slices: {h['slices']}\n")
    open(os.path.join("notes", f"n{n}_map.md"), "w").write("".join(L))
    print(json.dumps({k: v for k, v in summary.items() if k != "hits"}, indent=1))
    print(f"hits: {len(summary['hits'])}")


if __name__ == "__main__":
    main()
