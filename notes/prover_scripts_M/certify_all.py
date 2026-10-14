"""
EXACT re-runnable certificate for every counterexample claimed in notes_lemmaM.md. No floating point in any decision.
Writes certificates_lemmaM.json and exits non-zero if any check fails.

  C1  pure-stability Lemma M (M_s and M_w) fails, m = 3                         (notes §2.1)
  C2  Lemma M with positivity a(beta) > 0 fails, n = 7, &FK@dyUfXp?              (notes §2.2)  [FLAG: SFP of a known candidate]
  C3  graph-level pure form fails, n = 6, two dominated graphs                   (notes §2.3)
  C4  graph-level domination-free pure form fails, n = 7, two graphs; neither is an SFP anywhere (region2d)  (notes §2.4)
  C5  Lemma AS check: SFP beta-interval at finite t extends beyond beta_c (repo engine fastregion)          (notes §4)
"""
import sys, json
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, '/home/user/CLTN')
import ctln, graphs, fastregion, region2d
from certify_M_violation import from_d6, hurwitz_minors
from certify_graph_level_pure import omega_at

OUT = {}
FAIL = []

def check(name, cond):
    print(('PASS ' if cond else 'FAIL ') + name)
    if not cond:
        FAIL.append(name)

def point(A, beta, t):
    eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
    r = ctln.analyze_support(A, list(range(len(A))), eps, delta)
    return {"beta": str(beta), "t": str(t), "eps": str(eps), "delta": str(delta), "q": str(ctln.q_of(eps, delta)),
            "legal": ctln.legal(eps, delta), "positive": r["positive"], "stable": r["stable"],
            "stability_detail": r["stability_detail"], "SFP": r["stable_fixed_point"]}

def omega_signs(A, beta):
    om = omega_at(A, F(beta))
    H = hurwitz_minors(om)
    return [(h > 0) - (h < 0) for h in H], om[-1]

def strictly_not_hurwitz(signs, s0):
    # some minor < 0  =>  not Hurwitz ; Delta_{deg-1} != 0 and omega(0) = s != 0  =>  no root on iR  =>  a root with Re > 0
    return any(x < 0 for x in signs) and signs[-2] != 0 and s0 != 0

# ---------------------------------------------------------------- C1
A1 = [[0, 1, 1], [0, 0, 1], [1, 0, 0]]
b1 = F(22, 5)
p_lo, p_mid, p_hi = point(A1, b1, F(11, 25)), point(A1, b1, F(1)), point(A1, b1, F(100))
sg, s0 = omega_signs(A1, b1)
OUT["C1"] = {"adjacency": A1, "points": [p_lo, p_mid, p_hi], "omega_minor_signs_at_beta": sg, "s": str(s0)}
check("C1 stable at legal t=11/25", p_lo["legal"] and p_lo["stable"] is True)
check("C1 unstable at legal t=1 and t=100 (same beta)", p_mid["legal"] and p_mid["stable"] is False and p_hi["legal"] and p_hi["stable"] is False)
check("C1 limit not stable at beta, strictly (alpha > -(1-1/beta))", strictly_not_hurwitz(sg, s0))

# ---------------------------------------------------------------- C2
s2 = '&FK@dyUfXp?'; A2 = from_d6(s2); b2 = F(10561, 5000)
q_lo, q_hi = point(A2, b2, F(10)), point(A2, b2, F(1000))
sg, s0 = omega_signs(A2, b2)
OUT["C2"] = {"d6": s2, "adjacency": A2, "domination_free": graphs.is_domination_free(A2), "points": [q_lo, q_hi],
             "omega_minor_signs_at_beta": sg, "s": str(s0),
             "FLAG": "t=10 point is a legal finite-eps non-clique SFP of a graph already in COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl"}
check("C2 domination-free, not complete", graphs.is_domination_free(A2) and not graphs.is_complete_bidirected(A2))
check("C2 SFP (positive and stable) at legal t=10", q_lo["legal"] and q_lo["SFP"] is True)
check("C2 positive but unstable at legal t=1000 (same beta)", q_hi["legal"] and q_hi["positive"] and q_hi["stable"] is False)
check("C2 limit not stable at beta, strictly", strictly_not_hurwitz(sg, s0))

# ---------------------------------------------------------------- C3
OUT["C3"] = []
for s in ['&EODUpt{', '&EIKDw\\{']:
    A = from_d6(s); p = point(A, F(9, 4), F(11, 2)); sg, s0 = omega_signs(A, 2)
    OUT["C3"].append({"d6": s, "adjacency": A, "dominated_nodes": graphs.dominated_nodes(A), "point": p,
                      "omega_minor_signs_at_beta2": sg, "s2": str(s0)})
    check(f"C3 {s} stable at legal (9/4, 11/2)", p["legal"] and p["stable"] is True)
    check(f"C3 {s} alpha > -1/2 strictly", strictly_not_hurwitz(sg, s0))

# ---------------------------------------------------------------- C4
OUT["C4"] = []
for s, t in [('&FNALh\\Uwp?', F(19, 2)), ('&FM?_~LVd~?', F(21, 2))]:
    A = from_d6(s); beta = 2 + 1 / t + F(1, 10**6)
    p = point(A, beta, t); sg, s0 = omega_signs(A, 2)
    cert = region2d.region2d(A, list(range(len(A))), crosscheck=True)
    OUT["C4"].append({"d6": s, "adjacency": A, "domination_free": graphs.is_domination_free(A), "point": p,
                      "omega_minor_signs_at_beta2": sg, "s2": str(s0), "region2d_verdict": cert["verdict"],
                      "region2d_positivity_cells": cert.get("positivity_cells")})
    check(f"C4 {s} domination-free", graphs.is_domination_free(A))
    check(f"C4 {s} stable at a legal point", p["legal"] and p["stable"] is True)
    check(f"C4 {s} alpha > -1/2 strictly", strictly_not_hurwitz(sg, s0))
    check(f"C4 {s} region2d: no SFP anywhere in the legal plane", cert["verdict"] == "never")

# ---------------------------------------------------------------- C5
x = sp.symbols('x')
def right_end(A, t):
    cert = fastregion.region(A, list(range(len(A))), t)
    c = cert['stable_intervals'][-1]['hi']
    lo, hi = sp.Rational(c['beta_interval'][0]), sp.Rational(c['beta_interval'][1])
    ivs = [iv for iv, _ in sp.Poly(c['beta_minpoly_hi_to_lo'], x).intervals(eps=sp.Rational(1, 10**13)) if iv[0] >= lo and iv[1] <= hi]
    assert len(ivs) == 1
    return ivs[0]
bc = right_end(A2, fastregion.INF)
OUT["C5"] = {"d6": s2, "beta_c_interval": [str(bc[0]), str(bc[1])], "finite_t": {}}
for t in [F(10), F(100), F(1000)]:
    iv = right_end(A2, t)
    OUT["C5"]["finite_t"][str(t)] = [str(iv[0]), str(iv[1])]
    check(f"C5 right end of SFP beta-interval at t={t} > beta_c", iv[0] > bc[1])

json.dump(OUT, open('certificates_lemmaM.json', 'w'), indent=1, default=str)
print('FAILED:' if FAIL else 'ALL CHECKS PASS', FAIL)
sys.exit(1 if FAIL else 0)
