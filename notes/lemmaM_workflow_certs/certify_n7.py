"""EXACT certificate (ctln.analyze_support; Fractions only) for the n=7 candidate found by the float search
ext_scan.py (one-node extension of an n=6 stable graph).  Claims to be checked:
  (1) at beta = 100/49 and t = t1 (legal), the full support is a STABLE FIXED POINT (x > 0, Hurwitz stable)
  (2) at the same beta and larger t (legal), K is NOT stable (decided False)  -> strong form M+_s refuted
  (3) alpha(A|perp) >= -1/2 > -(1 - 1/beta): limit-unstable at this beta    -> weak form M+_w refuted
Re-run: python3 certify_n7.py   Writes cert_n7_M_violation.json."""
import sys, json
sys.path.insert(0, "/home/user/CLTN")
sys.path.insert(0, "/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf/lemmaM_construct")
from fractions import Fraction as F
from ctln import analyze_support, legal, q_of, is_clique, targets
from exact_ends import large_end_stable, positivity, chi_perp, shift_poly
from ctln import hurwitz_stable

ROWS = ["0110000", "0000111", "0001111", "1100100", "1001011", "1111101", "1111110"]
A = [[int(c) for c in r] for r in ROWS]
m = len(A); sig = list(range(m))
for r in ROWS:
    print(r)
beta = F(100, 49)
T1 = [F(20), F(30), F(35), F(40), F(45)]
T2 = [F(60), F(100), F(1000), F(10 ** 6)]
rec = {"adjacency_rows": ROWS, "beta": str(beta), "points": []}
for t in T1 + T2:
    eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
    r = analyze_support(A, sig, eps, delta)
    p = {"t": str(t), "eps": str(eps), "delta": str(delta), "q": str(q_of(eps, delta)), "legal": legal(eps, delta),
         "stable": r["stable"], "detail": r["stability_detail"], "positive": r["positive"],
         "stable_fixed_point": r["stable_fixed_point"], "x": r["x"]}
    rec["points"].append(p)
    print("t=%-8s eps=%-10s delta=%-14s q=%-12s legal=%s positive=%s stable=%s SFP=%s (%s)" %
          (t, eps, delta, q_of(eps, delta), p["legal"], p["positive"], p["stable"], p["stable_fixed_point"], p["detail"]))
b = 1 / beta
rec["limit_stable_at_beta"] = large_end_stable(A, b)
# alpha(A|perp) vs -1/2 exactly: chi_perp shifted by +1/2 (roots lam + 1/2), Hurwitz?
rec["alpha_perp_lt_minus_half"] = hurwitz_stable(shift_poly(chi_perp(A), F(-1, 2)))[0]
rec["chi_perp_hi_to_lo"] = [str(c) for c in chi_perp(A)]
pos, a = positivity(A, b)
rec["positivity"] = {"all_positive": pos, "a_scaled": a}
rec["clique"] = is_clique(A, sig)
print("limit-stable at beta:", rec["limit_stable_at_beta"], "| alpha(A|perp) < -1/2:", rec["alpha_perp_lt_minus_half"],
      "| chi_perp:", rec["chi_perp_hi_to_lo"], "| positivity:", pos, "| clique:", rec["clique"])
json.dump(rec, open("cert_n7_M_violation.json", "w"), indent=1)
