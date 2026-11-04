"""
certificate_check.py -- does the n = 7 proper-support classification
(certificates/n7_proper_supports_all.jsonl, produced by n7/run_n7_proper_supports.py
with region2d.region2d) carry the off-support polynomial o_6 as a wall, and does it
classify the example of Geneson's 2026-11-04 report correctly?

Steps (all exact):
 1. find an isomorphism phi from G_1 of Table 1 (paper labelling) to the
    stage-0 labelling of Geneson's graph (the base 'geneson_idx1000255' of the run);
 2. map the report's in-neighbourhood {2,3,5,6} (1-based) of the outside vertex
    through phi, find the record with that node6_in and empty node6_out;
 3. check (a) the record's positivity_critical list contains the irreducible
    o_6 = -b^5+3b^4-3b^3+3b^2-2b+1 (up to sign), whose isolating interval contains the
    unique root beta* of o_6 in [52/25, 209/100]; (b) 209/100 lies in a positivity cell and 52/25 in none;
    (c) the record has a stable (beta-cell x t-cell) containing (209/100, 100)
    and no stable cell has a beta-cell containing 52/25;
 4. over ALL 12,288 records: for each, recompute o_6 with fastregion.beta_data and
    check that every real root of o_6 in (2, inf) appears in the record's
    positivity_critical list (as a factor's isolating interval) -- i.e. the
    off-support curve is a wall in every certificate.
Run: cd /home/user/CLTN && nice -n 10 python3 notes/referee_external/fixes_2026-11-04/prover/certificate_check.py
"""
import itertools
import json
import os
import sys
from fractions import Fraction as F

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, ROOT)
from fastregion import beta_data  # noqa: E402

b = sp.symbols("beta")
G1 = [[int(c) for c in r] for r in ["011000", "000011", "000111", "110011", "110101", "101110"]]
S0 = json.load(open(os.path.join(ROOT, "certificates", "stage0_geneson_2607.21396.json")))["A_ours"]
REPORT_IN = [2, 3, 5, 6]
O7 = sp.Poly(-b**5 + 3*b**4 - 3*b**3 + 3*b**2 - 2*b + 1, b)


def isos(A, B):
    return [p for p in itertools.permutations(range(6))
            if all(A[i][j] == B[p[i]][p[j]] for i in range(6) for j in range(6))]


def poly_from_hi(c):
    return sp.Poly([sp.Integer(x) for x in c], b)


def in_interval(x, iv):
    return F(iv[0]) <= x <= F(iv[1])


def main():
    out = {}
    P = isos(G1, S0)
    out["isomorphisms_G1_to_stage0"] = [list(p) for p in P]
    assert P
    phi = P[0]
    target_in = sorted(phi[v - 1] for v in REPORT_IN)
    out["report_in_mapped_0based_stage0"] = target_in
    recs = [json.loads(l) for l in open(os.path.join(ROOT, "certificates", "n7_proper_supports_all.jsonl"))]
    out["n_records"] = len(recs)
    hit = [r for r in recs if r["base"] == "geneson_idx1000255" and r["node6_in"] == target_in and r["node6_out"] == []]
    assert len(hit) == 1, len(hit)
    r = hit[0]
    c = r["cert"]
    # (a)
    pc = c["positivity_critical"]
    o_entries = [e for e in pc if sp.Poly(poly_from_hi(e["minpoly_hi_to_lo"]).as_expr() + O7.as_expr(), b).is_zero
                 or sp.Poly(poly_from_hi(e["minpoly_hi_to_lo"]).as_expr() - O7.as_expr(), b).is_zero]
    out["record"] = {"pattern": r["pattern"], "adjacency_rows": r["adjacency_rows"], "verdict": r["verdict"],
                     "positivity_cells": c["positivity_cells"], "positivity_critical": pc,
                     "stable_cells": r["stable_cells"], "projection_names": c.get("projection_names")}
    out["o6_in_positivity_critical"] = o_entries
    seg = (sp.Rational(52, 25), sp.Rational(209, 100))
    out["o6_roots_in_segment"] = int(O7.count_roots(*seg))
    out["o6_isolating_interval_contains_segment_root"] = bool(o_entries) and all(
        O7.count_roots(max(sp.Rational(e["interval"][0]), seg[0]), min(sp.Rational(e["interval"][1]), seg[1])) == 1
        for e in o_entries)
    # (b)
    def in_pos(x):
        return any(F(lo) < x and (hi is None or x < F(hi)) for lo, hi in c["positivity_cells"])
    out["209/100_in_positivity_cell"] = in_pos(F(209, 100))
    out["52/25_in_positivity_cell"] = in_pos(F(52, 25))
    # (c)
    def stable_at(x, tt):
        for (blo, bhi, _bs, tlo, thi) in r["stable_cells"]:
            if F(blo) < x and (bhi is None or x < F(bhi)) and F(tlo) < tt and (thi is None or tt < F(thi)):
                return True
        return False
    out["stable_cell_contains_(209/100,100)"] = stable_at(F(209, 100), F(100))
    out["some_stable_beta_cell_contains_52/25"] = any(
        F(blo) < F(52, 25) and (bhi is None or F(52, 25) < F(bhi)) for (blo, bhi, *_r) in r["stable_cells"])
    # 4. every STABLE record (never-verdict records store no cert): o_6 real roots > 2 are walls
    bad, checked, nroots = [], 0, 0
    for rr in recs:
        if "cert" not in rr:
            continue
        A = [[int(ch) for ch in row] for row in rr["adjacency_rows"]]
        bd = beta_data(A, list(range(6)))
        o6 = bd["o"][6]
        checked += 1
        if not o6 or len(o6) < 2:
            continue
        po = sp.Poly(list(reversed([sp.Integer(x) for x in o6])), b)
        crit = rr["cert"]["positivity_critical"]
        for fac, _mult in sp.factor_list(po.as_expr())[1]:
            fp = sp.Poly(fac, b)
            if fp.degree() < 1:
                continue
            for (lo, hi), _ in fp.intervals(inf=2):
                lo, hi = sp.Rational(lo), sp.Rational(hi)
                if lo == hi == 2:
                    continue
                nroots += 1
                found = False
                for e in crit:
                    g = poly_from_hi(e["minpoly_hi_to_lo"])
                    if sp.Poly(g.as_expr() - fp.as_expr(), b).is_zero or sp.Poly(g.as_expr() + fp.as_expr(), b).is_zero:
                        elo, ehi = sp.Rational(e["interval"][0]), sp.Rational(e["interval"][1])
                        ilo, ihi = max(lo, elo), min(hi, ehi)
                        if ilo <= ihi and fp.count_roots(ilo, ihi) >= 1:
                            found = True
                            break
                if not found:
                    bad.append({"base": rr["base"], "pattern": rr["pattern"], "factor": str(fac), "interval": [str(lo), str(hi)]})
    out["stable_records_checked"] = checked
    out["o6_real_roots_above_2_total"] = nroots
    out["o6_roots_missing_from_critical"] = bad[:50]
    out["n_missing"] = len(bad)
    json.dump(out, open(os.path.join(HERE, "certificate_check.json"), "w"), indent=1)
    for k in ["isomorphisms_G1_to_stage0", "report_in_mapped_0based_stage0", "o6_in_positivity_critical",
              "o6_roots_in_segment", "o6_isolating_interval_contains_segment_root", "209/100_in_positivity_cell", "52/25_in_positivity_cell",
              "stable_cell_contains_(209/100,100)", "some_stable_beta_cell_contains_52/25",
              "stable_records_checked", "o6_real_roots_above_2_total", "n_missing"]:
        print(k, "=", out[k])
    print("record:", json.dumps({k: out["record"][k] for k in ["pattern", "adjacency_rows", "verdict", "positivity_cells", "stable_cells", "projection_names"]}))


if __name__ == "__main__":
    main()
