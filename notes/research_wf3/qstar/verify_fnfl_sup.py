"""
Independent exact checks of the FNFL supremum produced by qstar.py --fnfl --sharpen.
(The FNFL graph is a CANDIDATE pending human hand-verification; so is everything here.)

1. Approach from inside: at rational beta_k -> beta_hi (inside the stable beta-cell)
   recompute the top t-root tau_k of the stable t-cell with region2d's own isolation,
   form q_k = (beta_k - 1) tau_k/(tau_k + 1) exactly, check q_k < q_sup and that
   q_sup - q_k shrinks; certify a stable point (beta_k, t_k), t_k just below tau_k,
   with ctln.analyze_support (exact).  This does not use the resultant/candidate
   machinery at all.
2. Level tests at the endpoints of a 1e-15-wide isolating interval of q_sup:
   level_test(a) must be True (with a certified witness), level_test(b) False.
"""
import json, sys
from fractions import Fraction as F
sys.path.insert(0, "/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/qstar")
import qstar
from qstar import Alg, cmp, CellCtx, analyse_region2d, dec

R = json.load(open(sys.argv[1]))
d = json.load(open("/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json"))["graph"]
A = [[int(c) for c in r] for r in d["adjacency_rows"]]
an = analyse_region2d(d["region2d"], "FNFL")
cell = an["cells"][0]
ctx = CellCtx(A, range(7), d["region2d"], cell, "FNFL")
s = R["n7_finite_FNFL_CANDIDATE"]["cells"][0]["sharpened"]["sup"]
qs = Alg(s["minpoly_hi_to_lo"], F(s["isolating_interval"][0]), F(s["isolating_interval"][1]))
bh = Alg(cell["beta_hi"].p, cell["beta_hi"].a, cell["beta_hi"].b)
out = {"approach_from_inside": [], "level_tests": {}}
prev_gap = None
for k in (3, 5, 7, 9, 11):
    bh.refine_to(F(1, 10 ** k))
    b = bh.a
    if out["approach_from_inside"] and out["approach_from_inside"][-1]["beta_k"] == str(b):
        continue                      # interval already narrower than 10^-k: same beta_k
    low, up = ctx.section(b)
    qk = qstar.q_of_tau(up, b - 1)
    c = cmp(qk, qs)
    assert c < 0, "q at an inner slice exceeds the claimed supremum"
    qk2 = Alg(qk.p, qk.a, qk.b); qk2.refine_to(F(1, 10 ** 20))
    qs2 = Alg(qs.p, qs.a, qs.b); qs2.refine_to(F(1, 10 ** 20))
    gap = qs2.a - qk2.b
    # certified stable point just below the top root
    up.refine_to(F(1, 10 ** 12))
    lowv = low.hi() if low.is_rat() else (low.refine_to(F(1, 10 ** 12)) or low.b)
    t = up.a - F(1, 10 ** 9)
    assert t > lowv
    ok, w = ctx.certify_point(b, t)
    assert ok, w
    out["approach_from_inside"].append({"beta_k": str(b), "beta_hi_minus_beta_k_lt": f"1e-{k}",
                                        "q_k_approx": dec((qk2.a + qk2.b) / 2, 15),
                                        "q_sup_minus_q_k_approx": dec(gap, 15),
                                        "certified_stable_point_q_approx": dec(F(w["q"]), 15)})
    if prev_gap is not None:
        assert gap < prev_gap, "q at inner slices does not approach the supremum"

    prev_gap = gap
qs.refine_to(F(1, 10 ** 15))
ta, wa = ctx.level_test(qs.a)
tb, wb = ctx.level_test(qs.b)
assert ta is True and tb is False, (ta, tb)
out["level_tests"] = {"below": str(qs.a), "below_result": ta, "below_witness": wa, "above": str(qs.b), "above_result": tb}
print(json.dumps(out, indent=1))
print("FNFL SUP CHECKS PASSED")
