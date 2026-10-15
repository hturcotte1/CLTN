"""n8_instances_region2d.py -- exact two-parameter regions (region2d, m = 8) for three 8-node instances:
E_1 (finite-eps-only family, a = 1), F_2 (Geneson clique-class family, k = 2) and the incidental n = 8
limit candidate found by the overflow-guarded spectral filter. Candidates (n = 8 is not hand-verified).
Output: certificates/n8_instances_region2d.json. No floating point anywhere."""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from region2d import region2d
from fastregion import region, INF
from ctln import analyze_support, is_clique
from graphs import is_domination_free
from fractions import Fraction as F

base = "notes/research_wf3/theory/"
E = [r for r in json.load(open(base + "family_E_witnesses.json")) if r["a"] == 1][0]
Fk = [r for r in json.load(open(base + "F_family_witnesses.json")) if r["k"] == 2][0]
N8 = json.load(open("notes/research_wf3/filter8/CANDIDATE_n8_limit_GLjlNZ.json"))
def rows_to_A(rows): return [[int(c) for c in row] for row in rows]
inst = [("E_1 (n=8, finite-eps-only family)", rows_to_A(E["A_rows"]), E.get("witness")),
        ("F_2 (n=8, Geneson + clique class of 2)", rows_to_A(Fk["A_rows"]), {"eps": Fk["eps"], "delta": Fk["delta"]}),
        ("n8 limit candidate &GLjlNZ}NzH\\o", N8["A_as_found"] if isinstance(N8["A_as_found"][0], list) else rows_to_A(N8["A_as_found"]), None)]
out = []
for name, A, w in inst:
    m = len(A); t0 = time.time()
    lim = region(A, list(range(m)), INF, replay=True)
    r2 = region2d(A, list(range(m)), crosscheck=True, replay=True)
    rec = {"name": name, "n": m, "adjacency_rows": ["".join(map(str, r)) for r in A], "clique": bool(is_clique(A, list(range(m)))),
           "domination_free": bool(is_domination_free(A)), "limit_verdict": lim["verdict"], "limit_intervals": lim.get("stable_intervals"),
           "region2d_verdict": r2["verdict"], "stable_cells": r2["stable_cells"], "bounded_in_t_somewhere": any(c[4] is not None for c in r2["stable_cells"]),
           "seconds": round(time.time() - t0, 1), "region2d": r2}
    if w:
        s = analyze_support(A, list(range(m)), F(str(w["eps"])), F(str(w["delta"])))
        rec["stored_witness"] = {"eps": str(w["eps"]), "delta": str(w["delta"]), "stable": s["stable"], "positive": bool(s["x"] and all(F(str(v)) > 0 for v in s["x"]))}
    out.append(rec)
    print(name, "| limit", lim["verdict"], "| region2d", r2["verdict"], r2["stable_cells"], "| bounded", rec["bounded_in_t_somewhere"], "|", rec["seconds"], "s", flush=True)
    json.dump({"status": "CANDIDATES (n = 8; not hand-verified)", "instances": out}, open("certificates/n8_instances_region2d.json", "w"), indent=1)
