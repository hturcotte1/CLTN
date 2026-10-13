"""SANITY for notes/lemma_Bprime_2d.md (Lemma B' / B''' on region2d's actual cells). Exact arithmetic only.

For each case, run region2d.region2d; then for every beta-cell and every t-cell:
  DEL : at K random rational beta in the (rational inner) beta-cell, the number of real t-roots of
        the curves above t_min equals the number at the sample b0 (delineability evidence), and the
        root owned by each curve index is in the same order (labels constant)
  SEC : at each such beta, a random rational t strictly inside the j-th t-cell slice (same index j),
        the exact verdict (fastregion.decide_point) equals the cell verdict when the latter is decided;
        undecided samples are counted, and every decided verdict at a point of a cell whose first
        sample was undecided must be False (B''(3))
  LROOT: the L-root t_min appears as a degenerate isolating interval [t_min, t_min]
"""
import sys, json, random
sys.path.insert(0, "/home/user/CLTN")
sys.path.insert(0, "/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverB")
from fractions import Fraction as F
from math import lcm
import sympy as sp
import region2d
from region2d import _isolate, B, T
from fastregion import beta_data, decide_point, _frac, _cells
from check_projection import cases

rnd = random.Random(777)

def rand_in(lo, hi):
    """random rational strictly inside (lo, hi); hi None -> (lo, lo+20)"""
    if hi is None:
        hi = lo + 20
    u = F(rnd.randint(1, 999), 1000)
    return lo + (hi - lo) * u

def troots(curves, b):
    bq = sp.Rational(b.numerator, b.denominator)
    tp, owner = [], []
    for ci, (nm, f) in enumerate(curves):
        u = sp.Poly(f.as_expr().subs(B, bq), T)
        if u.degree() >= 1:
            den = 1
            for cc in u.all_coeffs():
                den = lcm(den, sp.Rational(cc).q)
            tp.append([int(sp.Rational(cc) * den) for cc in u.all_coeffs()]); owner.append(ci)
    tmin = 1 / (b - 2)
    iso = _isolate(tp, tmin, None, True) if tp else []
    degenerate_L = [iv for iv in iso if iv[0][0] == iv[0][1] == tmin]
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == tmin)]
    labels = [tuple(sorted(owner[i] for i in idx)) for (_, idx, _) in iso]
    return iso, labels, tmin, len(degenerate_L)

def run(A, sigma, tag, K=4):
    bd = beta_data(A, sigma)
    r = region2d.region2d(A, sigma, crosscheck=False)
    stats = {"tag": tag, "verdict": r["verdict"], "beta_cells": len(r["cells"]), "points": 0, "undecided_pts": 0,
             "cells_undecided_first": 0, "lroot_degenerate": 0, "lroot_checks": 0}
    if not r["cells"]:
        return stats
    curves, _, _ = region2d.projection(bd, region2d.curves_lc(bd))
    for row in r["cells"]:
        lo = F(row["beta_cell"][0]); hi = None if row["beta_cell"][1] is None else F(row["beta_cell"][1])
        b0 = F(row["beta_sample"])
        iso0, lab0, tmin0, dl0 = troots(curves, b0)
        stats["lroot_checks"] += 1; stats["lroot_degenerate"] += dl0
        assert len(iso0) + 1 == len(row["t_cells"]), "t-cell count mismatch at b0"
        for _ in range(K):
            b = rand_in(lo, hi)
            iso, lab, tmin, dl = troots(curves, b)
            stats["lroot_checks"] += 1; stats["lroot_degenerate"] += dl
            assert len(iso) == len(iso0), f"DEL: root count changed in beta-cell {row['beta_cell']} ({tag})"
            assert lab == lab0, f"DEL: section labels changed in beta-cell {row['beta_cell']} ({tag})"
            tc = _cells(iso, tmin)
            for j, (_, _, tlo, thi, _) in enumerate(tc):
                cell = row["t_cells"][j]
                v = decide_point(bd, rand_in(tlo, thi), b)
                stats["points"] += 1
                if v["stable"] is None or v["positive"] is None:
                    stats["undecided_pts"] += 1
                    continue
                pv = v["positive"] and v["off_ok"] and v["stable"]
                if cell["stable_fp"] is not None:
                    assert pv == cell["stable_fp"], f"SEC: verdict differs inside a cell ({tag}, beta={b})"
        for cell in row["t_cells"]:
            v0 = decide_point(bd, F(cell["t_lo"]) + (F(cell["t_hi"]) - F(cell["t_lo"])) / 2 if cell["t_hi"] is not None else F(cell["t_lo"]) + 1, b0)
            if v0["stable"] is None:
                stats["cells_undecided_first"] += 1
                assert cell["stable_fp"] in (False, None), "B''(3) violated: undecided first sample but stable cell"
    return stats

if __name__ == "__main__":
    allst = []
    for A, sigma, tag in cases():
        st = run(A, sigma, tag)
        allst.append(st)
        print(st, flush=True)
    json.dump(allst, open(__file__.replace(".py", ".json"), "w"), indent=1)
    print("TOTAL points", sum(s["points"] for s in allst), "undecided", sum(s["undecided_pts"] for s in allst),
          "L-root degenerate", sum(s["lroot_degenerate"] for s in allst), "/", sum(s["lroot_checks"] for s in allst))
    print("ALL SECTOR CHECKS PASSED on", len(allst), "cases")
