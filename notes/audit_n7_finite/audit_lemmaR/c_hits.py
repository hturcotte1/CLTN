"""curve_test / curve_test_fast on every n=7 limit hit (n7/n7_limit_HITS.jsonl), K7, the FNFL finite candidate,
and the four graphs named in notes/lemma_M_refuted.md sections 2.2 and 2.4. All must be stable_somewhere
(Lemma R '<=' direction applied to known legal-stable witnesses). Every stable cell gets an exact legal witness
with ctln.analyze_support at beta = 2 + 1/t + 10^-6. Also checks: the FNFL certificate's own witness points,
from_d6 vs the certificate's adjacency rows, and that the witness point projects (Lemma R '=>') to a stable curve
point lying in a stable cell. Output: c_hits.json"""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
from n7.curve_test import curve_test, curve_test_fast, from_d6, coeff_polys, coeffs_at
import ctln, sympy as sp

def witness(A, t, eta):
    beta = 2 + 1 / t + eta; eps = 1 / (t + 1); delta = (beta - 1) * eps
    assert ctln.legal(eps, delta)
    return ctln.analyze_support(A, list(range(7)), eps, delta)["stable"]

graphs = []
for line in open('/home/user/CLTN/n7/n7_limit_HITS.jsonl'):
    d = json.loads(line); graphs.append(("limit_hit", d["d6"]))
K7 = [[0 if i == j else 1 for j in range(7)] for i in range(7)]
cert = json.load(open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))
fn = cert["graph"]["digraph6_as_enumerated"]
rows = [[int(ch) for ch in r] for r in cert["graph"]["adjacency_rows"]]
out = {"fnfl_rows_match_d6": from_d6(fn) == rows}
graphs += [("K7", None), ("FNFL", fn), ("sec2.2", "&FK@dyUfXp?"), ("sec2.2b", "&FMEHz\\vpv?"),
           ("sec2.4a", "&FNALh\\Uwp?"), ("sec2.4b", "&FM?_~LVd~?")]
res = []; bad = []
for tag, s in graphs:
    A = K7 if s is None else from_d6(s)
    t0 = time.time(); P = coeff_polys(A); r = curve_test(A, P); rf = curve_test_fast(A)
    wit = []
    for c in r["cells"]:
        if c["stable"]:
            t = F(c["sample_t"]); w = witness(A, t, F(1, 10 ** 6)); k = 6
            while w is not True and k < 15:          # openness: some eta > 0 works; 1e-6 is not guaranteed
                k += 1; w = witness(A, t, F(1, 10 ** k))
            wit.append((c["sample_t"], w, k))
    rec = {"tag": tag, "d6": s, "verdict": r["verdict"], "fast": rf["verdict"], "n_critical": r["n_critical"],
           "stable_cells": [(c["lo"], c["hi"]) for c in r["cells"] if c["stable"]],
           "escalated_cells": sum(1 for c in r["cells"] if c["stable"] is None),
           "witness": wit, "witness_needed_eta_below_1e-6": any(k > 6 for _, _, k in wit), "sec": round(time.time() - t0, 3)}
    if r["verdict"] != "stable_somewhere" or rf["verdict"] != "stable_somewhere" or not all(w is True for _, w, _k in wit):
        bad.append(rec)
    res.append(rec)
# FNFL certificate witnesses: re-run analyze_support, then project along the ray to the curve (Lemma R '=>')
A = from_d6(fn); P = coeff_polys(A); r = curve_test(A, P)
proj = []
for w in cert["graph"]["witnesses"]["stable_positive_legal"]:
    eps, delta = F(w["eps"]), F(w["delta"]); t = (1 - eps) / eps; beta = 1 + delta / eps
    a = ctln.analyze_support(A, list(range(7)), eps, delta)
    # Lemma R '=>': t' = lambda* t with lambda* = (1 + sqrt(1 + beta/t))/beta, i.e. t' = t(1 + sqrt(1 + beta/t))/beta.
    # Exact rational bracketing: r < t'  <=>  r beta/t - 1 < sqrt(1 + beta/t)  <=>  (lhs < 0) or lhs^2 < 1 + beta/t.
    def below(rq):
        lhs = rq * beta / t - 1
        return lhs < 0 or lhs * lhs < 1 + beta / t
    lo_q = F(0); hi_q = F(1)
    while below(hi_q): hi_q *= 2
    while hi_q - lo_q > F(1, 10 ** 6):
        mid = (lo_q + hi_q) / 2
        if below(mid): lo_q = mid
        else: hi_q = mid
    assert below(lo_q) and not below(hi_q)
    inside = [c for c in r["cells"] if c["stable"] and F(c["lo"]) < F(str(lo_q)) and (c["hi"] is None or F(str(hi_q)) < F(c["hi"]))]
    proj.append({"eps": w["eps"], "delta": w["delta"], "stable_fp": a["stable_fixed_point"], "stable": a["stable"],
                 "curve_t_bracket": [str(lo_q), str(hi_q)], "bracket_inside_a_stable_cell": bool(inside)})
out.update({"n_graphs": len(res), "bad": bad, "fnfl_witness_projection": proj, "results": res})
json.dump(out, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_lemmaR/c_hits.json', 'w'), indent=1)
print(json.dumps({"n_graphs": len(res), "n_bad": len(bad), "fnfl_rows_match_d6": out["fnfl_rows_match_d6"], "fnfl_projection": proj}, indent=1))
print("graphs needing eta < 1e-6 for the witness:", [(x["d6"], x["witness"]) for x in res if x["witness_needed_eta_below_1e-6"]])
print("FNFL/2.2/2.4:", [(x["tag"], x["verdict"], x["stable_cells"]) for x in res if x["tag"] != "limit_hit"])
