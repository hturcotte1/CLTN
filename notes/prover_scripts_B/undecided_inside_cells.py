"""SANITY for Corollary B''(3): a point of a sign-invariant cell at which the exact Hurwitz test is
undecided (here: a zero coefficient c_l, l <= m-1, hit exactly -- c_l is linear in t) lies in a NON-stable
cell. For every region2d beta-cell sample b0 and every l = 1..m-1, put t* = -phi_l(b0)/omega_l(b0); if t* is
legal and not on a section, decide_point there must be undecided (or decided False), and the cell
containing (b0, t*) must have verdict False. A stable cell containing such a point would refute B''(3)."""
import sys, json
sys.path.insert(0, "/home/user/CLTN")
sys.path.insert(0, "/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverB")
from fractions import Fraction as F
import region2d
from fastregion import beta_data, decide_point, peval, _cells
from check_projection import cases
from sector_constancy import troots

tot = {"hits": 0, "undecided": 0, "in_stable_cell": 0, "cells_with_hit": 0, "decided_false": 0}
for A, sigma, tag in cases():
    bd = beta_data(A, sigma); m = bd["m"]
    r = region2d.region2d(A, sigma, crosscheck=False)
    if not r["cells"]:
        continue
    curves, _, _ = region2d.projection(bd, region2d.curves_lc(bd))
    for row in r["cells"]:
        b0 = F(row["beta_sample"])
        iso, lab, tmin, _ = troots(curves, b0)
        roots_lo_hi = [(F(int(a.numerator), int(a.denominator)), F(int(b.numerator), int(b.denominator))) for (a, b), _, _ in iso]
        for l in range(1, m):
            om = peval(bd["omega"][l], b0); ph = peval(bd["phi"][l], b0)
            if om == 0:
                continue
            ts = -ph / om
            if ts <= tmin:
                continue
            # which t-cell? (must be strictly between isolating intervals)
            j = sum(1 for (a, b) in roots_lo_hi if b < ts)
            if any(a <= ts <= b for (a, b) in roots_lo_hi):
                continue   # inside an isolating interval: cannot place it without refinement; skip
            cell = row["t_cells"][j]
            v = decide_point(bd, ts, b0)
            tot["hits"] += 1
            if v["stable"] is None:
                tot["undecided"] += 1
            elif v["stable"] is False or not (v["positive"] and v["off_ok"]):
                tot["decided_false"] += 1
            assert v["stable"] is not True, "zero coefficient but Hurwitz?!"
            if cell["stable_fp"]:
                tot["in_stable_cell"] += 1
                raise AssertionError(f"B''(3) refuted? {tag} b0={b0} t*={ts}")
print(tot)
print("B''(3) SANITY PASSED")
