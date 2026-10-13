"""
region2d.py -- the exact TWO-parameter stability region of one (graph, support).

Parameters (beta, t): beta = 1 + delta/eps > 2, t = (1-eps)/eps > 0, legal iff
L := t (beta - 2) - 1 > 0. Everything is polynomial in (beta, t) (fastregion.py).

Verdict at a point (Lemma A.2, notes/lemmas.md): sigma supports a stable fixed
point iff all a_i(beta) > 0, all o_k(beta) < 0, and p(mu) = det(mu I + K) is
Hurwitz. By Lemma B the verdict can only change on the zero sets of a_i, o_k,
D = det K (the constant term of p) and the Hurwitz minor Delta_{m-1} of p.

Method (cylindrical algebraic decomposition, projecting onto beta), restricted
to the positivity set P = {beta > 2 : a_i > 0, o_k < 0} (outside P nothing can
be stable):
  1. isolate the roots of the a_i, o_k on (2, inf); P = union of cells whose
     sample satisfies the sign conditions (exact).
  2. boundary curves in (beta, t): D (linear in t), Delta_{m-1} (degree m-1 in
     t) and L. Projection set in beta: leading coefficients in t, discriminants
     in t of the irreducible factors of Delta_{m-1}, pairwise resultants in t of
     all curves (made pairwise coprime first), plus the a_i, o_k.
  3. inside each P-interval, isolate the projection roots; on each open
     beta-cell take a rational sample beta0, isolate the roots in t of the
     curves at beta0 above t_min = 1/(beta0 - 2), and decide each (beta-cell x
     t-cell) at a rational sample with fastregion.decide_point (full exact
     Hurwitz test + positivity), cross-checked with ctln.analyze_support.
  Sign-invariance of all boundary polynomials on each cell (nonzero
  discriminants, resultants and leading coefficients on the open beta-cell)
  makes the sample verdict the verdict of the whole cell.
Replay (independent of sympy's VAS isolator): every isolating interval shows a
sign change of its irreducible factor, and the number of distinct real roots of
every projection polynomial in the considered beta-range is recounted with a
Descartes-bisection counter on the squarefree part (integer arithmetic).
"""
from fractions import Fraction as F
from math import gcd, lcm
import sympy as sp
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.sqfreetools import dup_sqf_part
from sympy.polys.domains import ZZ

import ctln
from fastregion import beta_data, decide_point, peval, _frac, _cells

B, T = sp.symbols("beta t")


def _poly_B(coeffs_low_to_high):
    return sum(sp.Integer(c) * B ** k for k, c in enumerate(coeffs_low_to_high))

# ------------------------------------------------------------- curves

def curves_lc(bd):
    """Boundary curves (Lemma B, notes/lemmas.md): D = det K (= constant term c_m of
    p, linear in t), the Hurwitz minor Delta_{m-1} of p (degree m-1 in t), and the
    legality boundary L. The verdict can change only on their zero sets (and on the
    vertical lines a_i = 0, o_k = 0). NOTE: Lienard-Chipart does not remove
    Delta_{m-1} (its halved minor sets always contain Delta_{m-1} or Delta_m =
    c_m Delta_{m-1}); an earlier draft that used Delta_{m-2} was wrong."""
    m = bd["m"]
    c = [_poly_B(bd["phi"][k]) + T * _poly_B(bd["omega"][k]) for k in range(m + 1)]   # c[0] = 1

    def ent(i, j):
        idx = 2 * i - j
        return c[idx] if 0 <= idx <= m else sp.Integer(0)

    k = m - 1
    H = sp.Matrix(k, k, lambda i, j: ent(i + 1, j + 1))
    curves = [("D", sp.Poly(c[m], B, T)),
              (f"Delta{k}", sp.Poly(sp.expand(H.det(method="berkowitz")), B, T)),
              ("L", sp.Poly(T * (B - 2) - 1, B, T))]
    return curves


def _coprime(curves, polys, names):
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
    return curves


def projection(bd, curves):
    polys, names = [], []
    # split every curve into irreducible factors (so discriminants are nonzero)
    split = []
    for nm, f in curves:
        if f.degree(T) <= 1:
            split.append((nm, f))
            continue
        _, facs = sp.factor_list(f.as_expr(), B, T)
        for i, (g, mult) in enumerate(facs):
            gp = sp.Poly(g, B, T)
            if gp.degree(T) == 0:
                if gp.degree(B) >= 1:
                    polys.append(sp.Poly(g, B)); names.append(f"vert_{nm}")
            else:
                split.append((f"{nm}_f{i}", gp))
    curves = _coprime(split, polys, names)
    for nm, f in curves:
        lc = sp.Poly(f.as_expr().coeff(T, f.degree(T)), B)
        if lc.degree() >= 1:
            polys.append(lc); names.append(f"lc_{nm}")
        if f.degree(T) >= 2:
            disc = sp.Poly(sp.discriminant(f.as_expr(), T), B)
            assert not disc.is_zero, f"{nm}: zero discriminant after factoring"
            if disc.degree() >= 1:
                polys.append(disc); names.append(f"disc_{nm}")
    for i in range(len(curves)):
        for j in range(i + 1, len(curves)):
            r = sp.Poly(sp.resultant(curves[i][1].as_expr(), curves[j][1].as_expr(), T), B)
            if r.is_zero:
                raise ArithmeticError(f"common factor between {curves[i][0]} and {curves[j][0]} after coprime split")
            if r.degree() >= 1:
                polys.append(r); names.append(f"res_{curves[i][0]}_{curves[j][0]}")
    return curves, polys, names

# ------------------------------------------------------------- Descartes-bisection root counting (replay)

def _shift1(p):
    """p(x+1), coefficients low->high, in place copy."""
    c = list(p)
    d = len(c) - 1
    for i in range(d):
        for j in range(d - 1, i - 1, -1):
            c[j] += c[j + 1]
    return c


def _variations(c):
    s = [v for v in c if v != 0]
    return sum(1 for i in range(len(s) - 1) if (s[i] > 0) != (s[i + 1] > 0))


def _div_x_minus_1(c):
    """exact division of c (low->high) by (x - 1); c(1) must be 0."""
    hi = list(reversed(c))
    out, acc = [], 0
    for v in hi:
        acc = acc + v
        out.append(acc)
    assert out[-1] == 0
    return list(reversed(out[:-1]))


def _count_open01(c, depth=0):
    """distinct roots of squarefree c (low->high ints) in the OPEN interval (0,1)."""
    while c and c[-1] == 0:
        c = c[:-1]
    if len(c) <= 1:
        return 0
    v = _variations(_shift1(list(reversed(c))))      # Descartes bound for (0,1)
    if v == 0:
        return 0
    if v == 1:
        return 1
    if depth > 400:
        raise ArithmeticError("Descartes bisection did not terminate (non-squarefree input?)")
    d = len(c) - 1
    left = [v_ << (d - k) for k, v_ in enumerate(c)]      # 2^d c(x/2)   roots of c in (0,1/2) -> (0,1)
    right = _shift1(left)                                   # 2^d c((x+1)/2)
    cnt = 0
    if right[0] == 0:                                       # c(1/2) == 0
        cnt = 1
        left = _div_x_minus_1(left)
        right = right[1:]
    return cnt + _count_open01(left, depth + 1) + _count_open01(right, depth + 1)


def count_roots_closed(p_hi_to_lo, a, b):
    """distinct real roots of integer polynomial p in [a, b] (a < b rationals); b None = +inf.
    Uses the squarefree part (sympy gcd) and Descartes bisection (own code)."""
    s = dup_sqf_part([int(v) for v in p_hi_to_lo], ZZ)
    c = [int(v) for v in reversed(s)]                        # low->high
    d = len(c) - 1
    if d < 1:
        return 0
    a = F(a)
    if b is None:
        bound = 1 + max(abs(F(v)) for v in c[:-1]) / abs(c[-1])      # Cauchy: all roots have |x| < bound
        b = max(F(int(bound) + 1), a + 1)
    b = F(b)
    # x = (A + (Bn - A) z) / Dn  with common denominator Dn
    Dn = lcm(a.denominator, b.denominator)
    A_, Bn = a.numerator * (Dn // a.denominator), b.numerator * (Dn // b.denominator)
    # q(z) = Dn^d c(x) = sum_k c_k Dn^(d-k) (A + (Bn-A) z)^k
    q = [0] * (d + 1)
    powA = [1]
    span = Bn - A_
    binom = [1]
    for k, ck in enumerate(c):
        if k > 0:
            binom = [1] + [binom[i] + binom[i + 1] for i in range(len(binom) - 1)] + [1]
        # (A + span z)^k = sum_j C(k,j) A^(k-j) span^j z^j
        scale = ck * Dn ** (d - k)
        if scale:
            for j in range(k + 1):
                q[j] += scale * binom[j] * A_ ** (k - j) * span ** j
    cnt = 0
    if q[0] == 0:                       # root at a
        cnt += 1
        q = q[1:]
    if sum(q) == 0:                     # root at b: q(1) == 0
        cnt += 1
        q = _div_x_minus_1(q)
    return cnt + _count_open01(q)


def replay_counts(iso, polys_hi_to_lo, lo, hi):
    """Check sympy's isolation of `polys` on [lo, hi]: sign change (or exact root) in
    every interval, and per-polynomial distinct-root counts equal to our own."""
    per = {}
    for (a_, c_), idx, f in iso:
        a_, c_ = _frac(a_), _frac(c_)
        flo = [F(int(v)) for v in reversed(f)]
        if a_ == c_:
            assert peval(flo, a_) == 0, "degenerate isolating interval is not a root"
        else:
            assert peval(flo, a_) * peval(flo, c_) < 0, "isolating interval shows no sign change"
        for i in idx:
            per[i] = per.get(i, 0) + 1
    for i, p in enumerate(polys_hi_to_lo):
        if len(p) <= 1:
            assert per.get(i, 0) == 0
            continue
        n_true = count_roots_closed(p, lo, hi)
        assert n_true == per.get(i, 0), f"replay: polynomial {i} has {n_true} roots in [{lo},{hi}] but {per.get(i, 0)} isolating intervals"

# ------------------------------------------------------------- main

def _isolate(polys_hi_to_lo, lo, hi, replay):
    iso = dup_isolate_real_roots_list(polys_hi_to_lo, ZZ, inf=lo, sup=hi, strict=True, basis=True)
    if replay:
        replay_counts(iso, polys_hi_to_lo, lo, hi)
    return iso


def positivity_intervals(bd, replay=True):
    """Cells of beta > 2 where all a_i > 0 and all o_k < 0, as rational hulls
    [lo, hi] (hi None = unbounded), plus the raw isolation for the certificate."""
    polys = [list(reversed(a)) for a in bd["a"]] + [list(reversed(o)) for o in bd["o"].values()]
    polys = [[int(v) for v in p] for p in polys if len(p) > 1]
    b_min = F(2)
    iso = _isolate(polys, b_min, None, replay) if polys else []
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == b_min)]
    cells = _cells(iso, b_min)
    P = []
    for lo_i, hi_i, lo, hi, s in cells:
        if all(peval(a, s) > 0 for a in bd["a"]) and all(peval(o, s) < 0 for o in bd["o"].values()):
            P.append((lo_i, hi_i, lo, hi))
    return P, iso


def region2d(A, sigma, crosscheck=True, verbose=False, replay=True):
    sigma = sorted(sigma)
    bd = beta_data(A, sigma)
    m = bd["m"]
    P, isoP = positivity_intervals(bd, replay)
    out = {"sigma": sigma, "method": "CAD on {a_i, o_k, D, Delta_{m-1}, L} restricted to the positivity set; Descartes replay",
           "positivity_cells": [(str(lo), None if hi is None else str(hi)) for _, _, lo, hi in P],
           "positivity_critical": [{"interval": [str(_frac(a)), str(_frac(b))], "minpoly_hi_to_lo": [int(c) for c in f]} for (a, b), _, f in isoP],
           "cells": [], "stable_cells": [], "projection_names": None, "beta_critical": []}
    if not P:
        out["verdict"] = "never"
        out["reason"] = "positivity set empty"
        return out
    curves0 = curves_lc(bd)
    curves, polys, names = projection(bd, curves0)
    out["projection_names"] = names
    pos_polys = [[int(v) for v in reversed(a)] for a in bd["a"]] + [[int(v) for v in reversed(o)] for o in bd["o"].values()]
    pos_polys = [p for p in pos_polys if len(p) > 1]
    npos = len(pos_polys)
    intp = pos_polys + [[int(c) for c in p.all_coeffs()] for p in polys if p.degree() >= 1]
    escal = False
    b_min = F(2)
    for (plo_i, phi_i, plo, phi) in P:
        lo_r = b_min if plo_i == -1 else _frac(isoP[plo_i][0][0])
        hi_r = None if phi_i is None else _frac(isoP[phi_i][0][1])
        iso = _isolate(intp, lo_r, hi_r, replay)
        if plo_i == -1:
            iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == b_min)]
        ivs = [(_frac(a), _frac(b), set(idx), f) for (a, b), idx, f in iso]
        is_pos = [any(i < npos for i in idx) for (a, b, idx, f) in ivs]
        pos_pos = [k for k, flag in enumerate(is_pos) if flag]
        # the positivity-cell endpoints are the first (unless the cell starts at 2) and last
        # (unless unbounded) positivity roots inside the hull
        k_start = -1 if plo_i == -1 else pos_pos[0]
        k_end = len(ivs) if phi_i is None else pos_pos[-1]
        assert (plo_i == -1 or len(pos_pos) >= 1) and (phi_i is None or len(pos_pos) >= 1), "positivity endpoint not found in joint isolation"
        for a_, b_, idx, f in ivs:
            out["beta_critical"].append({"interval": [str(a_), str(b_)], "minpoly_hi_to_lo": [int(c) for c in f],
                                         "positivity_root": any(i < npos for i in idx)})
        prev_hi = b_min if k_start == -1 else ivs[k_start][1]
        ivs = [list(iv) for iv in ivs]
        for k in range(k_start + 1, k_end + 1):
            if k < len(ivs):
                # an isolating interval may start exactly at the hull end (its root is
                # strictly above it): refine it until the cell below it is nonempty
                guard = 0
                while ivs[k][0] <= prev_hi and ivs[k][0] != ivs[k][1]:
                    from sympy.polys.rootisolation import dup_refine_real_root
                    a_n, b_n = dup_refine_real_root(ivs[k][3], ivs[k][0], ivs[k][1], ZZ, steps=1)
                    ivs[k][0], ivs[k][1] = _frac(a_n), _frac(b_n)
                    guard += 1
                    assert guard < 500, "could not refine interval above the hull end"
                cell_lo, cell_hi = prev_hi, ivs[k][0]
                assert cell_lo < cell_hi, "empty beta-cell"
                b0 = (cell_lo + cell_hi) / 2
            else:
                cell_lo, cell_hi = prev_hi, None
                b0 = prev_hi + 1
            tpolys = []
            bq = sp.Rational(b0.numerator, b0.denominator)
            for nm, f in curves:
                u = sp.Poly(f.as_expr().subs(B, bq), T)
                if u.degree() >= 1:
                    den = 1
                    for cc in u.all_coeffs():
                        den = lcm(den, sp.Rational(cc).q)
                    tpolys.append([int(sp.Rational(cc) * den) for cc in u.all_coeffs()])
            t_min = 1 / (b0 - 2)
            tiso = _isolate(tpolys, t_min, None, replay) if tpolys else []
            tiso = [iv for iv in tiso if not (iv[0][0] == iv[0][1] == t_min)]
            tcells = _cells(tiso, t_min)
            row = {"beta_cell": (str(cell_lo), None if cell_hi is None else str(cell_hi)), "beta_sample": str(b0), "t_cells": [],
                   "t_critical": [{"interval": [str(_frac(a)), str(_frac(b))], "minpoly_hi_to_lo": [int(c) for c in f]} for (a, b), _, f in tiso]}
            for tlo_i, thi_i, tlo, thi, t0 in tcells:
                v = decide_point(bd, t0, b0)
                if v["stable"] is None or v["positive"] is None:
                    for s0 in ([tlo + (thi - tlo) / 3, tlo + 2 * (thi - tlo) / 3] if thi is not None else [tlo + 2, tlo + 5]):
                        v = decide_point(bd, s0, b0)
                        if v["stable"] is not None and v["positive"] is not None:
                            t0 = s0
                            break
                rec = {"t_lo": str(tlo), "t_hi": None if thi is None else str(thi), "t_sample": str(t0), **v}
                if crosscheck:
                    eps = 1 / (t0 + 1); delta = (b0 - 1) * eps
                    r = ctln.analyze_support(A, sigma, eps, delta)
                    agree = (r["positive"] == v["positive"] and r["off_support_ok"] == v["off_ok"] and r["stable"] == v["stable"])
                    rec["crosscheck"] = agree
                    if not agree:
                        raise AssertionError(f"region2d disagrees with analyze_support at beta={b0}, t={t0}")
                verdict = (v["positive"] and v["off_ok"] and v["stable"]) if v["stable"] is not None else None
                rec["stable_fp"] = verdict
                if verdict is None:
                    escal = True
                row["t_cells"].append(rec)
                if verdict:
                    out["stable_cells"].append((str(cell_lo), None if cell_hi is None else str(cell_hi), str(b0), rec["t_lo"], rec["t_hi"]))
            out["cells"].append(row)
            if verbose:
                print(f"beta in ({cell_lo},{cell_hi}) sample {b0}: " + ", ".join(
                    f"t in ({c['t_lo']},{c['t_hi']}):{'S' if c['stable_fp'] else ('?' if c['stable_fp'] is None else '-')}" for c in row["t_cells"]))
            if k < len(ivs):
                prev_hi = ivs[k][1]
    out["verdict"] = "stable_somewhere" if out["stable_cells"] else ("escalate" if escal else "never")
    return out


if __name__ == "__main__":
    import time
    G = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
    t0 = time.time(); r = region2d(G, range(6), verbose=True)
    print("verdict:", r["verdict"], f"{time.time()-t0:.1f}s")
    print("positivity cells:", r["positivity_cells"])
    print("stable cells:", r["stable_cells"])
    print("beta critical:", [(c["interval"], len(c["minpoly_hi_to_lo"]) - 1) for c in r["beta_critical"]])
