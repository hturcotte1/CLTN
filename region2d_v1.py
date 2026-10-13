"""
region2d.py -- the exact TWO-parameter stability region of one (graph, support).

Parameters (beta, t) with beta = 1 + delta/eps > 2, t = (1-eps)/eps > 0, legal
iff L := t (beta - 2) - 1 > 0. Everything is a polynomial in (beta, t) because
I - W_sigma = eps (C(beta) + t J):  see fastregion.py.

Method: cylindrical algebraic decomposition, projecting onto beta.
  boundary polynomials      a_i(beta), o_k(beta)      (univariate, "vertical" boundaries)
                            D = d + t s               (linear in t)
                            L = t (beta-2) - 1        (legal boundary)
                            Delta_{m-1}(beta, t)      (degree <= m-1 in t)
  projection set in beta    every univariate boundary polynomial; for each
                            irreducible factor f of Delta with deg_t f >= 1:
                            its leading coefficient in t and (if deg_t f >= 2)
                            its discriminant in t; the resultant in t of every
                            pair among {factors of Delta, D, L}; lc_t(D) = s.
Over each open beta-cell between consecutive real roots of the projection set,
every boundary polynomial has a constant number of distinct real roots in t,
these roots vary continuously and never cross (nonzero discriminants and
resultants, nonzero leading coefficients), so the sign pattern of all boundary
polynomials -- hence the verdict -- is constant on each (beta-cell x t-cell).
Each such cell is decided at one rational sample point with
fastregion.decide_point (exact), and cross-checked with ctln.analyze_support.
"""
from fractions import Fraction as F
import sympy as sp
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.domains import ZZ

import ctln
from fastregion import beta_data, decide_point, peval, hurwitz_det_poly, _frac, _cells, _replay_isolation

B, T = sp.symbols("beta t")


def _poly_B(coeffs_low_to_high):
    return sum(sp.Integer(c) * B ** k for k, c in enumerate(coeffs_low_to_high))


def bivariate_data(bd):
    """sympy expressions in (beta, t) for D, Delta, and the stability coefficients."""
    m = bd["m"]
    P = [_poly_B(bd["phi"][k]) + T * _poly_B(bd["omega"][k]) for k in range(m + 1)]   # highest degree first
    deg = m

    def ent(i, j):
        idx = 2 * i - j
        return P[idx] if 0 <= idx <= deg else sp.Integer(0)

    k = deg - 1
    H = sp.Matrix(k, k, lambda i, j: ent(i + 1, j + 1))
    Delta = sp.Poly(sp.expand(H.det(method="berkowitz")), B, T)
    D = sp.Poly(_poly_B(bd["d"]) + T * _poly_B(bd["s"]), B, T)
    L = sp.Poly(T * (B - 2) - 1, B, T)
    return {"P": P, "Delta": Delta, "D": D, "L": L}


def projection_polys(bd, bv):
    """Integer polynomials in beta whose real roots delimit the beta-cells."""
    polys = []
    names = []
    for i, a in enumerate(bd["a"]):
        polys.append(sp.Poly(_poly_B(a), B)); names.append(f"a{i}")
    for kk, o in bd["o"].items():
        polys.append(sp.Poly(_poly_B(o), B)); names.append(f"o{kk}")
    # factors of Delta
    _, facs = sp.factor_list(bv["Delta"].as_expr(), B, T)
    tfacs = []
    for f, mult in facs:
        pf = sp.Poly(f, B, T)
        if pf.degree(T) == 0:
            polys.append(sp.Poly(f, B)); names.append("Delta_vert_factor")
        else:
            tfacs.append(pf)
    curves = [("D", bv["D"]), ("L", bv["L"])] + [(f"Delta_f{i}", f) for i, f in enumerate(tfacs)]
    # Make the curve set pairwise coprime: a common factor h of two curves is split
    # off as its own curve (its roots bound both), so every resultant below is nonzero.
    changed = True
    while changed:
        changed = False
        for i in range(len(curves)):
            for j in range(i + 1, len(curves)):
                g = sp.gcd(curves[i][1].as_expr(), curves[j][1].as_expr())
                gp = sp.Poly(g, B, T)
                if gp.total_degree() >= 1:
                    ni, fi = curves[i]; nj, fj = curves[j]
                    qi = sp.Poly(sp.cancel(fi.as_expr() / g), B, T)
                    qj = sp.Poly(sp.cancel(fj.as_expr() / g), B, T)
                    new = [(n, f) for k, (n, f) in enumerate(curves) if k not in (i, j)]
                    for n, f in ((ni, qi), (nj, qj), (f"gcd_{ni}_{nj}", gp)):
                        if f.total_degree() >= 1:
                            if f.degree(T) == 0:
                                polys.append(sp.Poly(f.as_expr(), B)); names.append(f"vert_{n}")
                            else:
                                new.append((n, f))
                    curves = new
                    changed = True
                    break
            if changed:
                break
    for nm, f in curves:
        lc = sp.Poly(f.as_expr().coeff(T, f.degree(T)), B)
        if lc.degree() >= 1:
            polys.append(lc); names.append(f"lc_{nm}")
        if f.degree(T) >= 2:
            disc = sp.Poly(sp.discriminant(f.as_expr(), T), B)
            assert not disc.is_zero, f"{nm}: zero discriminant in t (repeated factor) -- unexpected after factoring"
            polys.append(disc); names.append(f"disc_{nm}")
    for i in range(len(curves)):
        for j in range(i + 1, len(curves)):
            r = sp.Poly(sp.resultant(curves[i][1].as_expr(), curves[j][1].as_expr(), T), B)
            if r.is_zero:
                raise ArithmeticError(f"common factor between {curves[i][0]} and {curves[j][0]}")
            if r.degree() >= 1:
                polys.append(r); names.append(f"res_{curves[i][0]}_{curves[j][0]}")
    return polys, names, curves


def region2d(A, sigma, crosscheck=True, verbose=False, replay=True):
    sigma = sorted(sigma)
    bd = beta_data(A, sigma)
    bv = bivariate_data(bd)
    polys, names, curves = projection_polys(bd, bv)
    b_min = F(2)
    intp = [[int(c) for c in p.all_coeffs()] for p in polys if p.degree() >= 1]
    iso = dup_isolate_real_roots_list(intp, ZZ, inf=b_min, strict=True, basis=True)
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == b_min)]
    if replay:
        _replay_isolation(iso, intp, b_min)
    bcells = _cells(iso, b_min)
    crit_beta = [{"interval": [str(_frac(a)), str(_frac(b))], "minpoly_hi_to_lo": [int(c) for c in f]} for (a, b), _, f in iso]

    out = []
    stable_cells = []
    for lo_i, hi_i, lo, hi, b0 in bcells:
        t_min = 1 / (b0 - 2)
        # univariate polynomials in t at beta = b0 (integer coefficients after clearing denominators)
        tpolys = []
        for nm, f in curves:
            u = sp.Poly(f.as_expr().subs(B, sp.Rational(b0.numerator, b0.denominator)), T)
            u = sp.Poly(u.as_expr() * sp.ilcm(*[sp.Rational(c).q for c in u.all_coeffs()]), T) if u.degree() >= 0 else u
            tpolys.append((nm, [int(c) for c in sp.Poly(u, T).all_coeffs()] if u.degree() >= 1 else None))
        tiso = dup_isolate_real_roots_list([c for _, c in tpolys if c], ZZ, inf=t_min, strict=True, basis=True)
        tiso = [iv for iv in tiso if not (iv[0][0] == iv[0][1] == t_min)]
        if replay:
            _replay_isolation(tiso, [c for _, c in tpolys if c], t_min)
        tcells = _cells(tiso, t_min)
        row = {"beta_cell": (lo_i, hi_i), "beta_sample": str(b0), "t_cells": []}
        for tlo_i, thi_i, tlo, thi, t0 in tcells:
            v = decide_point(bd, t0, b0)
            if v["stable"] is None or v["positive"] is None:
                # resample inside the same cell
                for s0 in ([tlo + (thi - tlo) / 3, tlo + 2 * (thi - tlo) / 3] if thi is not None else [tlo + 2, tlo + 5]):
                    v = decide_point(bd, s0, b0)
                    if v["stable"] is not None and v["positive"] is not None:
                        t0 = s0
                        break
            rec = {"t_cell": (tlo_i, thi_i), "t_lo": str(tlo), "t_hi": None if thi is None else str(thi),
                   "t_sample": str(t0), **v}
            if crosscheck:
                eps = 1 / (t0 + 1)
                delta = (b0 - 1) * eps
                r = ctln.analyze_support(A, sigma, eps, delta)
                agree = (r["positive"] == v["positive"] and r["off_support_ok"] == v["off_ok"] and r["stable"] == v["stable"])
                rec["crosscheck"] = agree
                if not agree:
                    raise AssertionError(f"region2d disagrees with analyze_support at beta={b0}, t={t0}")
            verdict = (v["positive"] and v["off_ok"] and v["stable"]) if v["stable"] is not None else None
            rec["stable_fp"] = verdict
            row["t_cells"].append(rec)
            if verdict:
                stable_cells.append((str(lo), None if hi is None else str(hi), str(b0), rec["t_lo"], rec["t_hi"]))
        row["t_critical"] = [{"interval": [str(_frac(a)), str(_frac(b))], "minpoly_hi_to_lo": [int(c) for c in f]} for (a, b), _, f in tiso]
        out.append(row)
        if verbose:
            print(f"beta in cell {lo_i}..{hi_i} (sample {b0}): " + ", ".join(
                f"t in ({c['t_lo']},{c['t_hi']}):{'S' if c['stable_fp'] else ('?' if c['stable_fp'] is None else '-')}" for c in row["t_cells"]))
    return {"sigma": sigma, "beta_critical": crit_beta, "projection_names": names, "cells": out,
            "stable_cells": stable_cells, "Delta": str(bv["Delta"].as_expr()) if len(str(bv["Delta"].as_expr())) < 4000 else "(long)",
            "verdict": "stable_somewhere" if stable_cells else ("escalate" if any(c["stable_fp"] is None for r in out for c in r["t_cells"]) else "never")}


if __name__ == "__main__":
    import time, json
    A = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
    t0 = time.time()
    r = region2d(A, range(6), verbose=True)
    print("verdict:", r["verdict"], f"{time.time()-t0:.1f}s")
    print("beta critical values:", [(c["interval"], len(c["minpoly_hi_to_lo"]) - 1) for c in r["beta_critical"]])
    print("projection polys:", r["projection_names"])
    print("stable cells (beta_lo, beta_hi, beta_sample, t_lo, t_hi):")
    for c in r["stable_cells"]:
        print("   ", c)
