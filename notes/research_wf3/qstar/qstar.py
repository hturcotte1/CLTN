"""
qstar.py -- EXACT supremum of q over the stable regions of non-clique full supports,
read from region2d certificates (CLAUDE.md conventions; no floating point anywhere
in a decision; decimals are printed only as labelled approximations of exact
rational enclosures).

Coordinates (fastregion.py / region2d.py): beta = 1 + delta/eps > 2,
t = (1-eps)/eps > 0, q = (beta-1) t/(t+1), legal iff t(beta-2) > 1.
For fixed beta, q is strictly increasing in t; for fixed t it is increasing in beta;
q < beta - 1 always.

WHAT IS COMPUTED, AND WHY IT IS SOUND
------------------------------------
region2d decomposes the legal (beta, t) region into open 2-cells
    C = { beta_lo < beta < beta_hi,  s_{j-1}(beta) < t < s_j(beta) }
(s_j continuous sections of the boundary curves; s_j = +inf for the top cell) and
decides each by an exact sample.  A certificate's stable_cells entry
(lo, hi, b0, t_lo, t_hi) is the rational hull of one stable 2-cell: lo/hi are the
nearest endpoints of the isolating intervals of beta_lo/beta_hi, t_lo/t_hi those of
the t-roots s_{j-1}(b0), s_j(b0) at the sample b0 (t_hi = None iff s_j = +inf).

(S0) sup over the stable set = max over stable 2-cells of sup over the cell.
     The stable set is open (Lemma O, notes/lemma_F_reduction.md section 7), so every
     stable point has stable points of a 2-cell arbitrarily close (2-cells are dense),
     and the verdict is constant on a 2-cell; q is continuous.  Undecided 2-cells
     (stable_fp None) are reported separately; with any present the result is only a
     lower bound.
(S1) UPPER bound for a cell: q < beta - 1 < beta_hi - 1 on C.  beta_hi is the
     smallest beta_critical root >= hi (no critical root lies in the open hull
     (lo, hi), checked).  beta_hi = None (cell unbounded in beta) gives +inf.
(S2) If t_hi = None the bound is attained as a supremum: for every beta in the
     beta-cell, (beta, t) is in C for all large t, and q -> beta - 1 as t -> inf;
     letting beta -> beta_hi gives sup_C q = beta_hi - 1 EXACTLY.
(S3) If t_hi is finite: LOWER bound (b0-1) tau/(tau+1), tau = s_j(b0) the t-root
     at the sample (points (b0, t), t -> tau-, lie in C).  This is exact algebraic.
(S4) Optional sharpening (needs the graph; --graph or key in the record): the exact
     sup of Q(beta) = (beta-1) s_j/(s_j+1) over the beta-cell, by
       * a decision procedure "sup_C q > ell ?" for rational ell: the section s_j can
         cross the level curve tau_ell(beta) = ell/(beta-1-ell) only at real roots of
         G_f(beta, ell), G_f(beta, q) := (beta-1-q)^d f(beta, q/(beta-1-q)), f the
         irreducible curve carrying s_j (fixed on the cell by CAD delineability); on
         each root-free subinterval one rational sample decides; a YES comes with an
         exact certified stable point (decide_point and ctln.analyze_support) whose
         q exceeds ell;
       * a finite candidate set Z containing sup_C q: if Res_t(f, h) with
         h = t(t+1) f_t - (beta-1) f_beta has no root in the beta-cell, Q is strictly
         monotone there (Q' = (q_beta f_t - q_t f_beta)/f_t and f_t != 0 on the
         cell), so sup = an endpoint limit, which is a root of
         Res_beta(G_f(beta, q), m(beta)) for m the minimal polynomial of beta_hi
         (resp. beta_lo), or beta_hi - 1 (section escaping to infinity);
       * bisection on ell until exactly one candidate of Z lies in (lo, hi].
Aggregation: exact maximum of the per-cell sups (and of the upper bounds) over all
cells of all certificates.

Comparisons of algebraic numbers: isolating intervals refined with sympy's
dup_refine_real_root (every refinement re-checked by an exact sign change with
our own Horner evaluation; fallback exact bisection) until disjoint; equal
irreducible polynomials with overlapping intervals are identified by an exact
root count of the hull.

Usage:
  python3 qstar.py --n6                      # the three n=6 HITS files
  python3 qstar.py --fnfl                    # the finite-eps n=7 candidate
  python3 qstar.py --n7-limit                # lower bound from limit certificates
  python3 qstar.py --region2d-jsonl FILE...  # any *_HITS.jsonl
  python3 qstar.py --selftest
  add --sharpen to compute exact sups of t-bounded cells (needs the graph)
"""
import argparse
import json
import os
import sys
from fractions import Fraction as F

REPO = os.environ.get("CLTN_REPO", "/home/user/CLTN")
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.rootisolation import dup_refine_real_root

X = sp.Symbol("x")
Bs, Ts, Qs = sp.symbols("beta t q")

# ======================================================================
# exact univariate helpers (coefficients hi -> lo, Python ints / Fractions)
# ======================================================================


def horner(p, x):
    v = 0
    for c in p:
        v = v * x + c
    return v


def sgn(v):
    return (v > 0) - (v < 0)


def _F(x):
    if isinstance(x, F):
        return x
    if isinstance(x, int):
        return F(x)
    if isinstance(x, str):
        return F(x)
    return F(int(x.numerator), int(x.denominator))


def prim_int(poly_expr_or_poly, var=X):
    """Integer primitive coefficient list (hi->lo, positive leading coeff)."""
    P = poly_expr_or_poly if isinstance(poly_expr_or_poly, sp.Poly) else sp.Poly(poly_expr_or_poly, var)
    cs = [sp.Rational(c) for c in P.all_coeffs()]
    den = sp.ilcm(*[c.q for c in cs]) if cs else 1
    ints = [int(c * den) for c in cs]
    g = 0
    for c in ints:
        g = sp.igcd(g, c)
    ints = [c // g for c in ints] if g else ints
    if ints and ints[0] < 0:
        ints = [-c for c in ints]
    return ints


def irreducible_factors(p):
    """Irreducible factors over ZZ (as primitive int lists) of positive degree."""
    _, facs = sp.factor_list(sp.Poly(p, X))
    out = []
    for f, _m in facs:
        fi = prim_int(f)
        if len(fi) >= 2:
            out.append(fi)
    return out


def count_roots(p, a, b):
    """exact number of distinct real roots of p in the CLOSED interval [a, b] (Sturm, sympy)."""
    return int(sp.Poly(p, X).count_roots(sp.Rational(a.numerator, a.denominator),
                                         sp.Rational(b.numerator, b.denominator)))


# ======================================================================
# real algebraic numbers
# ======================================================================

class Alg:
    """A real algebraic number: either rational (self.rat) or the unique root of an
    irreducible primitive integer polynomial self.p (deg >= 2) in the OPEN interval
    (self.a, self.b) with a sign change.  All fields exact."""

    def __init__(self, p=None, a=None, b=None, rat=None):
        if rat is not None:
            self.rat, self.p, self.a, self.b = _F(rat), [1, -_F(rat)], _F(rat), _F(rat)
            return
        self.rat = None
        self.p, self.a, self.b = list(p), _F(a), _F(b)
        assert len(self.p) >= 3, "Alg: degree >= 2 expected"
        assert self.a < self.b
        assert sgn(horner(self.p, self.a)) * sgn(horner(self.p, self.b)) < 0, "Alg: no sign change"

    # ---- construction from (any integer poly, closed interval) ----
    @staticmethod
    def from_interval(poly, a, b):
        """The unique real root of `poly` in the closed interval [a, b]; `poly` need not
        be irreducible.  Picks the irreducible factor that vanishes there and checks
        uniqueness exactly."""
        a, b = _F(a), _F(b)
        poly = [int(c) for c in poly]
        if a == b:
            assert horner(poly, a) == 0, "degenerate interval is not a root"
            return Alg(rat=a)
        hits = []
        for f in irreducible_factors(poly):
            k = count_roots(f, a, b)
            if k:
                hits.append((f, k))
        assert len(hits) == 1 and hits[0][1] == 1, f"interval [{a},{b}] does not isolate exactly one root: {hits}"
        f = hits[0][0]
        if len(f) == 2:
            return Alg(rat=F(-f[1], f[0]))
        return Alg(f, a, b)       # irreducible deg >= 2: no rational roots, endpoints not roots

    def is_rat(self):
        return self.rat is not None

    def lo(self):
        return self.rat if self.is_rat() else self.a

    def hi(self):
        return self.rat if self.is_rat() else self.b

    def refine(self):
        if self.is_rat():
            return
        a, b = self.a, self.b
        try:
            na, nb = dup_refine_real_root(self.p, a, b, ZZ, steps=1)
            na, nb = _F(na), _F(nb)
            ok = (a <= na < nb <= b and (nb - na) < (b - a)
                  and sgn(horner(self.p, na)) * sgn(horner(self.p, nb)) < 0)
        except Exception:
            ok = False
        if not ok:                                   # exact bisection fallback
            m = (a + b) / 2
            sm = sgn(horner(self.p, m))
            assert sm != 0
            if sm == sgn(horner(self.p, a)):
                na, nb = m, b
            else:
                na, nb = a, m
        self.a, self.b = na, nb

    def refine_to(self, width):
        while not self.is_rat() and self.b - self.a > width:
            self.refine()

    def cmp_rat(self, r):
        r = _F(r)
        if self.is_rat():
            return sgn(self.rat - r)
        assert horner(self.p, r) != 0
        while self.a <= r <= self.b:
            self.refine()
        return 1 if self.a > r else -1

    def __repr__(self):
        if self.is_rat():
            return f"Alg({self.rat})"
        return f"Alg(root of {self.p} in ({self.a}, {self.b}))"

    def to_json(self, digits=15):
        if self.is_rat():
            return {"rational": str(self.rat), "approx_decimal_NOT_A_CERTIFICATE": dec(self.rat, digits)}
        c = Alg(self.p, self.a, self.b)
        c.refine_to(F(1, 10 ** (digits + 2)))
        return {"minpoly_hi_to_lo": self.p, "isolating_interval": [str(self.a), str(self.b)],
                "degree": len(self.p) - 1,
                "approx_decimal_NOT_A_CERTIFICATE": dec((c.a + c.b) / 2, digits)}


def dec(r, digits):
    r = _F(r)
    s = "-" if r < 0 else ""
    r = abs(r)
    ip = r.numerator // r.denominator
    frac = r - ip
    ds = (frac * 10 ** digits).numerator * 1 // (frac * 10 ** digits).denominator
    return f"{s}{ip}.{str(ds).zfill(digits)}"


def cmp(x, y):
    """exact comparison of two Alg; returns -1, 0, 1."""
    if x.is_rat():
        return -y.cmp_rat(x.rat)
    if y.is_rat():
        return x.cmp_rat(y.rat)
    if x.p == y.p and not (x.b < y.a or y.b < x.a):
        lo, hi = min(x.a, y.a), max(x.b, y.b)
        if count_roots(x.p, lo, hi) == 1:
            return 0
    guard = 0
    while not (x.b < y.a or y.b < x.a):
        if x.p == y.p:
            lo, hi = min(x.a, y.a), max(x.b, y.b)
            if count_roots(x.p, lo, hi) == 1:
                return 0
        x.refine(); y.refine()
        guard += 1
        assert guard < 5000, "comparison did not separate"
    return -1 if x.b < y.a else 1


def minpoly_int(z):
    if z.is_rat():
        return [z.rat.denominator, -z.rat.numerator]
    return z.p


def alg_max(items, key=lambda z: z):
    best = None
    for it in items:
        if best is None or cmp(key(it), key(best)) > 0:
            best = it
    return best


def transform(x, num, den, image_lo, image_hi):
    """y = M(x) for a monotone rational map M on [x.lo, x.hi]; the minimal polynomial of
    y is found from  P(y) := den(y)^d * p(num(y)/den(y))  where x = num(y)/den(y) is the
    inverse map (num, den sympy expressions in X).  image_lo/hi: M at the interval ends."""
    if x.is_rat():
        raise ValueError("use rational arithmetic")
    d = len(x.p) - 1
    expr = sum(sp.Integer(c) * num ** (d - k) * den ** k for k, c in enumerate(x.p))
    P = prim_int(sp.expand(expr))
    a, b = (image_lo, image_hi) if image_lo < image_hi else (image_hi, image_lo)
    return Alg.from_interval(P, a, b)


def shift(x, c):
    """x - c for rational c."""
    c = _F(c)
    if x.is_rat():
        return Alg(rat=x.rat - c)
    cq = sp.Rational(c.numerator, c.denominator)
    return transform(x, X + cq, sp.Integer(1), x.a - c, x.b - c)


def q_of_tau(tau, c):
    """c * tau/(tau+1) for rational c > 0 and tau > 0."""
    c = _F(c)
    if tau.is_rat():
        return Alg(rat=c * tau.rat / (tau.rat + 1))
    assert tau.a > -1
    cq = sp.Rational(c.numerator, c.denominator)
    # tau = y/(c - y)
    return transform(tau, X, cq - X, c * tau.a / (tau.a + 1), c * tau.b / (tau.b + 1))


# ======================================================================
# per-certificate cell analysis (no graph needed)
# ======================================================================

def _crit_list(entries):
    return [Alg.from_interval(e["minpoly_hi_to_lo"], *e["interval"]) for e in entries]


def analyse_region2d(cert, label):
    """Per stable 2-cell: exact upper bound U = beta_hi - 1, exact lower bound L, flags."""
    crit = _crit_list(cert.get("beta_critical", []))
    rows = cert.get("cells", [])
    out = {"label": label, "cells": [], "undecided_cells": [], "verdict": cert.get("verdict")}
    for row in rows:
        for tc in row["t_cells"]:
            if tc.get("stable_fp") is None:
                out["undecided_cells"].append({"beta_cell": row["beta_cell"], "beta_sample": row["beta_sample"],
                                               "t_lo": tc["t_lo"], "t_hi": tc["t_hi"]})
    for sc in cert.get("stable_cells", []):
        lo, hi, b0, tlo, thi = sc
        lo_r, b0_r = F(lo), F(b0)
        hi_r = None if hi is None else F(hi)
        # locate the row and the t-cell; certificate self-consistency
        row = [r for r in rows if r["beta_sample"] == b0 and list(r["beta_cell"]) == [lo, hi]]
        assert len(row) == 1, f"{label}: no unique row for stable cell {sc}"
        row = row[0]
        jj = [k for k, tc in enumerate(row["t_cells"]) if tc["t_lo"] == tlo and tc["t_hi"] == thi]
        assert len(jj) == 1 and row["t_cells"][jj[0]]["stable_fp"] is True, f"{label}: t-cell mismatch {sc}"
        j = jj[0]
        assert lo_r < b0_r and (hi_r is None or b0_r < hi_r)
        # no critical root strictly inside the open hull (lo, hi)
        for r in crit:
            inside = r.cmp_rat(lo_r) > 0 and (hi_r is None or r.cmp_rat(hi_r) < 0)
            assert not inside, f"{label}: critical root {r} inside the hull of stable cell {sc}"
        if hi_r is None:
            beta_hi = None
        else:
            above = [r for r in crit if r.cmp_rat(hi_r) >= 0]
            assert above, f"{label}: no critical root >= {hi} for cell {sc}"
            beta_hi = above[0]
            for r in above[1:]:
                if cmp(r, beta_hi) < 0:
                    beta_hi = r
        below = [r for r in crit if r.cmp_rat(lo_r) <= 0]
        beta_lo = Alg(rat=2)
        for r in below:
            if cmp(r, beta_lo) > 0:
                beta_lo = r
        U = None if beta_hi is None else shift(beta_hi, 1)
        rec = {"beta_cell_hull": [lo, hi], "beta_sample": b0, "t_hull": [tlo, thi], "t_index": j,
               "beta_lo": beta_lo, "beta_hi": beta_hi, "U": U, "t_bounded": thi is not None}
        if thi is None:
            rec["L"] = U
            rec["exact"] = U is not None
            rec["L_reason"] = "t unbounded above on the cell: sup = beta_hi - 1 (S2)"
        else:
            tcrit = _crit_list(row["t_critical"])
            thi_r = F(thi)
            cand = [r for r in tcrit if r.cmp_rat(thi_r) >= 0]
            assert cand, f"{label}: t_hi root not found"
            tau = cand[0]
            for r in cand[1:]:
                if cmp(r, tau) < 0:
                    tau = r
            rec["tau_sample"] = tau
            rec["L"] = q_of_tau(tau, b0_r - 1)
            rec["exact"] = False
            rec["L_reason"] = "sample slice: (b0-1) tau/(tau+1), tau = top t-root at b0 (S3)"
        out["cells"].append(rec)
    return out


def aggregate(analyses):
    cells = [(a["label"], c) for a in analyses for c in a["cells"]]
    undecided = [(a["label"], u) for a in analyses for u in a["undecided_cells"]]
    if not cells:
        return {"n_stable_cells": 0}
    unb = [lc for lc in cells if lc[1]["U"] is None]
    res = {"n_stable_cells": len(cells), "n_certificates": len(analyses),
           "undecided_cells": [{"label": l, **u} for l, u in undecided]}
    if unb:
        res["upper"] = "+inf (a stable cell is unbounded in beta)"
        res["upper_attained_by"] = [l for l, _ in unb]
    else:
        lu, cu = cells[0]
        for l, c in cells[1:]:
            if cmp(c["U"], cu["U"]) > 0:
                lu, cu = l, c
        res["upper"] = cu["U"].to_json()
        res["upper_attained_by"] = sorted({l for l, c in cells if cmp(c["U"], cu["U"]) == 0})
        res["upper_beta_hi"] = cu["beta_hi"].to_json()
    ll, cl = cells[0]
    for l, c in cells[1:]:
        if cmp(c["L"], cl["L"]) > 0:
            ll, cl = l, c
    res["lower"] = cl["L"].to_json()
    res["lower_attained_by"] = sorted({l for l, c in cells if cmp(c["L"], cl["L"]) == 0})
    if not unb:
        eq = cmp(cu["U"], cl["L"]) == 0
        res["exact"] = eq and not undecided
        if eq:
            res["q_star"] = cu["U"].to_json()
    return res


# ======================================================================
# limit certificates (t = infinity): a LOWER bound only
# ======================================================================

def analyse_limit(rec, label):
    """Limit-stable on (beta_lo, beta_hi) at t = inf  =>  for each such beta, stable for
    all sufficiently large t (p = phi + t omega, omega Hurwitz of degree m-1 with
    positive leading coefficient m, phi monic of degree m: one root ~ -m t, the rest
    near the roots of omega; positivity/off-support conditions are t-free signs once
    sign(D) = sign(s)), so sup over the finite region >= beta_hi - 1.  Not an upper
    bound: finite-eps regions may extend beyond the limit set."""
    cert = rec["cert"]
    out = []
    for iv in cert.get("stable_intervals", []):
        h = iv["hi"]
        if h["kind"] == "inf":
            out.append(None)
            continue
        if h["kind"] == "root":
            b = Alg.from_interval(h["beta_minpoly_hi_to_lo"], *h["beta_interval"])
        else:
            raise ValueError(f"{label}: unknown endpoint kind {h['kind']}")
        out.append(shift(b, 1))
    return out


# ======================================================================
# sharpening of t-bounded cells (needs the graph)
# ======================================================================

class CellCtx:
    def __init__(self, A, sigma, cert, cell, label, verbose=False):
        import fastregion, region2d
        self.fr, self.r2 = fastregion, region2d
        self.A, self.sigma, self.label, self.verbose = A, sorted(sigma), label, verbose
        self.bd = fastregion.beta_data(A, self.sigma)
        curves, _, _ = region2d.projection(self.bd, region2d.curves_lc(self.bd))
        self.curves = curves
        self.cell = cell
        self.j = cell["t_index"]
        row = [r for r in cert["cells"] if r["beta_sample"] == cell["beta_sample"]][0]
        self.n_tcells = len(row["t_cells"])
        self.b0 = F(cell["beta_sample"])
        # reconstruct the slice at b0 and check it reproduces the certificate
        tiso, cells = self.slice(self.b0)
        got = [{"interval": [str(_F(a)), str(_F(b))], "minpoly_hi_to_lo": [int(c) for c in f]} for (a, b), _, f in tiso]
        assert got == row["t_critical"], f"{label}: recomputed t-slice at b0 differs from the certificate"
        assert len(cells) == self.n_tcells
        # identify the curve carrying the top section s_j
        tau = cell["tau_sample"]
        owners = []
        for k, (nm, f) in enumerate(self.curves):
            u = self._tpoly(f, self.b0)
            if u is not None and sp.rem(sp.Poly(u, X), sp.Poly(minpoly_int(tau), X)).is_zero:
                owners.append(k)
        assert len(owners) == 1, f"{label}: top section owner ambiguous {owners}"
        self.k_top = owners[0]
        self.f_top = self.curves[owners[0]][1]
        self.top_name = self.curves[owners[0]][0]
        d = self.f_top.degree(Ts)
        fe = self.f_top.as_expr()
        coeffs = [sp.expand(fe).coeff(Ts, k) for k in range(d + 1)]
        self.G = sp.Poly(sp.expand(sum(coeffs[k] * Qs ** k * (Bs - 1 - Qs) ** (d - k) for k in range(d + 1))), Bs, Qs)
        self.witnesses = []

    def _tpoly(self, f, b):
        u = sp.Poly(f.as_expr().subs(Bs, sp.Rational(b.numerator, b.denominator)), Ts)
        if u.degree() < 1:
            return None
        return prim_int(u.as_expr().subs(Ts, X))

    def slice(self, b):
        tpolys = [p for p in (self._tpoly(f, b) for _, f in self.curves) if p is not None]
        t_min = 1 / (b - 2)
        tiso = self.r2._isolate(tpolys, t_min, None, True) if tpolys else []
        tiso = [iv for iv in tiso if not (_F(iv[0][0]) == _F(iv[0][1]) == t_min)]
        cells = self.fr._cells(tiso, t_min)
        return tiso, cells

    def section(self, b):
        """(lower boundary, upper boundary) of the j-th t-cell at rational b, as Alg
        (upper None = +inf)."""
        tiso, cells = self.slice(b)
        assert len(cells) == self.n_tcells, f"{self.label}: t-cell count changed inside the beta-cell (delineability broken?)"
        lo_i, hi_i, _, _, _ = cells[self.j]
        t_min = 1 / (b - 2)
        low = Alg(rat=t_min) if lo_i == -1 else Alg.from_interval(tiso[lo_i][2], _F(tiso[lo_i][0][0]), _F(tiso[lo_i][0][1]))
        up = None if hi_i is None else Alg.from_interval(tiso[hi_i][2], _F(tiso[hi_i][0][0]), _F(tiso[hi_i][0][1]))
        return low, up

    def certify_point(self, b, t):
        import ctln
        assert t * (b - 2) > 1, "not legal"
        v = self.fr.decide_point(self.bd, t, b)
        eps = 1 / (t + 1)
        delta = (b - 1) * eps
        r = ctln.analyze_support(self.A, self.sigma, eps, delta)
        ok = v["positive"] and v["off_ok"] and v["stable"] is True and r["positive"] and r["off_support_ok"] and r["stable"] is True
        return ok, {"beta": str(b), "t": str(t), "eps": str(eps), "delta": str(delta), "q": str((b - 1) * t / (t + 1)),
                    "legal": t * (b - 2) > 1, "decide_point": [v["positive"], v["off_ok"], v["stable"]],
                    "analyze_support": [r["positive"], r["off_support_ok"], r["stable"], r.get("stability_detail")]}

    def _between(self, u, v):
        """a rational strictly between Alg u < Alg v."""
        while u.hi() >= v.lo():
            u.refine(); v.refine()
            if u.is_rat() and v.is_rat():
                break
        return (u.hi() + v.lo()) / 2

    def _roots_in_cell(self, P_int):
        lo_alg, hi_alg = self.cell["beta_lo"], self.cell["beta_hi"]
        a, b = lo_alg.lo(), hi_alg.hi()
        roots = []
        for fac in irreducible_factors(P_int):
            for (ra, rb), _mult in sp.Poly(fac, X).intervals(inf=sp.Rational(a.numerator, a.denominator),
                                                       sup=sp.Rational(b.numerator, b.denominator)):
                z = Alg.from_interval(fac, _F(ra), _F(rb))
                if cmp(z, lo_alg) > 0 and cmp(z, hi_alg) < 0:
                    roots.append(z)
        return roots

    def level_test(self, ell):
        """True (+ certified witness with q > ell) iff some point of the cell has q > ell."""
        ell = _F(ell)
        el = sp.Rational(ell.numerator, ell.denominator)
        P = prim_int(sp.expand(self.G.as_expr().subs(Qs, el)).subs(Bs, X))
        assert any(P), "G(beta, ell) vanishes identically"
        brk = [self.cell["beta_lo"], self.cell["beta_hi"]] + self._roots_in_cell(P)
        one = Alg(rat=1 + ell)
        if cmp(one, self.cell["beta_lo"]) > 0 and cmp(one, self.cell["beta_hi"]) < 0:
            brk.append(one)
        # sort exactly, dedupe
        srt = []
        for z in brk:
            k = 0
            while k < len(srt) and cmp(srt[k], z) < 0:
                k += 1
            if k < len(srt) and cmp(srt[k], z) == 0:
                continue
            srt.insert(k, z)
        for u, v in zip(srt, srt[1:]):
            s = self._between(u, v)
            if s - 1 - ell <= 0:
                continue
            tau_ell = ell / (s - 1 - ell)
            low, up = self.section(s)
            up_gt = True if up is None else up.cmp_rat(tau_ell) > 0
            if up_gt:
                # witness t strictly inside the t-cell and above tau_ell
                if up is None:
                    t = max(tau_ell, low.hi()) + 1
                else:
                    guard = 0
                    while max(tau_ell, low.hi()) >= up.lo():
                        low.refine(); up.refine()
                        guard += 1
                        assert guard < 5000
                    t = (max(tau_ell, low.hi()) + up.lo()) / 2
                ok, w = self.certify_point(s, t)
                assert ok, f"{self.label}: level witness failed exact certification: {w}"
                assert F(w["q"]) > ell
                return True, w
        return False, None

    def candidates(self):
        """finite set Z (list of Alg q-values) containing sup over the cell, valid when
        Q has no interior critical point; also returns the interior-critical report."""
        f = self.f_top
        fe = f.as_expr()
        h = sp.expand(Ts * (Ts + 1) * sp.diff(fe, Ts) - (Bs - 1) * sp.diff(fe, Bs))
        R = sp.Poly(sp.resultant(fe, h, Ts), Bs)
        assert not R.is_zero, "Res_t(f, h) vanishes identically"
        Rint = prim_int(R.as_expr().subs(Bs, X)) if R.degree() >= 1 else None
        interior = self._roots_in_cell(Rint) if Rint else []
        Z = []
        for end in ("beta_lo", "beta_hi"):
            e = self.cell[end]
            if e.is_rat():
                Pq = prim_int(sp.expand(self.G.as_expr().subs(Bs, sp.Rational(e.rat.numerator, e.rat.denominator))).subs(Qs, X))
            else:
                mb = sum(sp.Integer(c) * Bs ** (len(e.p) - 1 - k) for k, c in enumerate(e.p))
                Pq = prim_int(sp.Poly(sp.resultant(self.G.as_expr(), mb, Bs), Qs).as_expr().subs(Qs, X))
            Z.append((end, "root of G(beta_end, q)", Pq))
            Z.append((end, "beta_end - 1", shift(e, 1)))
        return Z, interior, (R.degree() if R.degree() >= 0 else None)

    def sharpen(self, max_iter=400, min_width=F(1, 10 ** 12)):
        cell = self.cell
        # starting bracket: lo from the sample slice (certified), hi from U
        ok, w0 = self.level_test(cell["L"].lo() - F(1, 10 ** 6))
        assert ok
        lo = F(w0["q"])
        self.witnesses.append(w0)
        U = Alg(cell["U"].p, cell["U"].a, cell["U"].b) if not cell["U"].is_rat() else cell["U"]
        hi = U.hi()
        Z, interior, degR = self.candidates()
        info = {"top_curve": self.top_name, "Res_t(f,h)_degree": degR,
                "interior_critical_betas_in_cell": [z.to_json(8) for z in interior]}
        if interior:
            info["status"] = "interior critical points present: candidate set incomplete; bounds only"
        # build candidate Alg list lazily within [lo, hi]
        def cand_in(lo_, hi_):
            got = []
            for end, kind, obj in Z:
                if isinstance(obj, Alg):
                    zs = [obj]
                else:
                    zs = []
                    for fac in irreducible_factors(obj):
                        for (ra, rb), _mult in sp.Poly(fac, X).intervals(inf=sp.Rational(lo_.numerator, lo_.denominator),
                                                                   sup=sp.Rational(hi_.numerator, hi_.denominator)):
                            zs.append(Alg.from_interval(fac, _F(ra), _F(rb)))
                for z in zs:
                    if z.cmp_rat(lo_) > 0 and z.cmp_rat(hi_) <= 0:
                        if not any(cmp(z, y[0]) == 0 for y in got):
                            got.append((z, end, kind))
            return got
        import time as _time
        t_start = _time.time()
        def log(msg):
            if self.verbose:
                print(f"  [{self.label} {(_time.time() - t_start):7.1f}s] {msg}", file=sys.stderr, flush=True)
        # phase 1: bisection on dyadic levels (lo := mid on YES keeps denominators small;
        # sup > mid is what a YES proves) until exactly one candidate is left in (lo, hi]
        it = 0
        while True:
            cs = cand_in(lo, hi)
            log(f"step {it}: bracket width ~2^{(hi - lo).denominator.bit_length() - (hi - lo).numerator.bit_length()}, candidates {len(cs)}")
            if not interior:
                assert len(cs) >= 1, f"{self.label}: no candidate in ({lo}, {hi}] -- candidate set unsound?"
                if len(cs) == 1:
                    break
            if it >= max_iter or (interior and hi - lo < min_width):
                break
            mid = (lo + hi) / 2
            ok, w = self.level_test(mid)
            if ok:
                lo = mid
                self.witnesses.append(w)
            else:
                hi = mid
            it += 1
        info["bisection_steps"] = it
        info["unique_candidate_bracket"] = [str(lo), str(hi)]
        if not interior and len(cs) == 1:
            z, end, kind = cs[0]
            # phase 2: two-sided confirmation at the ends of a min_width isolating interval of z
            zz = Alg(z.p, z.a, z.b) if not z.is_rat() else z
            if not zz.is_rat():
                zz.refine_to(min_width)
                a2, b2 = zz.a, zz.b
            else:
                a2, b2 = zz.rat - min_width, zz.rat
            okA, wA = self.level_test(a2)
            log(f"confirm below: {okA}")
            okB, _ = self.level_test(b2)
            log(f"confirm above: {okB}")
            assert okA is True and okB is False, f"{self.label}: two-sided confirmation failed ({okA}, {okB})"
            cs2 = cand_in(a2, b2)
            assert len(cs2) == 1 and cmp(cs2[0][0], z) == 0
            self.witnesses.append(wA)
            info["status"] = "exact"
            info["sup"] = zz
            info["attained_as"] = f"limit at {end}: {kind}"
            info["confirmation"] = {"level_below": str(a2), "sup_gt_level_below": okA, "certified_witness_above_level_below": wA,
                                    "level_above": str(b2), "sup_gt_level_above": okB}
        else:
            info["bracket"] = [str(lo), str(hi)]
        info["seconds"] = round(_time.time() - t_start, 1)
        return info


# ======================================================================
# drivers
# ======================================================================

def load_n6_hits():
    out = []
    for s in range(3):
        p = os.path.join(REPO, "certificates", f"region2d_n6_shard{s}of3_HITS.jsonl")
        with open(p) as fh:
            for line in fh:
                if line.strip():
                    out.append(json.loads(line))
    return out


def graph_from_record(rec, n):
    if "A" in rec:
        return [[int(v) for v in row] for row in rec["A"]]
    from graphs import from_key
    from enum_cache import key_from_int
    return from_key(key_from_int(rec["key"], n), n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n6", action="store_true")
    ap.add_argument("--fnfl", action="store_true")
    ap.add_argument("--n7-limit", action="store_true")
    ap.add_argument("--region2d-jsonl", nargs="*", default=[])
    ap.add_argument("--n", type=int, default=None, help="node count for key decoding (jsonl inputs)")
    ap.add_argument("--sharpen", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    report = {}

    def run_region2d(records, tag, n):
        analyses, graphs_, cert_of = [], {}, {}
        for rec in records:
            assert rec.get("complete") is False, "complete (clique) graph in a HITS input"
            lab = f"idx{rec['idx']}" if "idx" in rec else str(rec.get("d6", rec.get("key")))
            analyses.append(analyse_region2d(rec["cert"], lab))
            cert_of[lab] = rec["cert"]
            if n is not None and ("A" in rec or "key" in rec):
                A = graph_from_record(rec, n)
                from graphs import is_complete_bidirected
                assert not is_complete_bidirected(A)
                graphs_[lab] = A
        agg = aggregate(analyses)
        percell = []
        for a in analyses:
            for c in a["cells"]:
                e = {"label": a["label"], "beta_cell_hull": c["beta_cell_hull"], "beta_sample": c["beta_sample"],
                     "t_hull": c["t_hull"], "upper_beta_hi_minus_1": c["U"].to_json() if c["U"] else "+inf",
                     "lower": c["L"].to_json() if c["L"] else None, "exact_sup": c["exact"], "lower_reason": c["L_reason"]}
                if args.sharpen and c["t_bounded"] and a["label"] in graphs_:
                    ctx = CellCtx(graphs_[a["label"]], range(n), cert_of[a["label"]], c, a["label"], verbose=True)
                    info = ctx.sharpen()
                    if "sup" in info:
                        info["sup"] = info["sup"].to_json()
                    e["sharpened"] = info
                percell.append(e)
        report[tag] = {"aggregate": agg, "cells": percell}

    if args.n6:
        run_region2d(load_n6_hits(), "n6_full_support", 6)
    for p in args.region2d_jsonl:
        with open(p) as fh:
            recs = [json.loads(l) for l in fh if l.strip()]
        run_region2d(recs, os.path.basename(p), args.n)
    if args.fnfl:
        d = json.load(open(os.path.join(REPO, "certificates", "COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json")))
        g = d["graph"]
        A = [[int(ch) for ch in row] for row in g["adjacency_rows"]]
        from graphs import is_complete_bidirected
        assert not is_complete_bidirected(A)
        an = analyse_region2d(g["region2d"], "FNFL &FNFL~[udo?")
        agg = aggregate([an])
        cells = []
        for c in an["cells"]:
            e = {"beta_cell_hull": c["beta_cell_hull"], "beta_sample": c["beta_sample"], "t_hull": c["t_hull"],
                 "upper_beta_hi_minus_1": c["U"].to_json(), "lower_sample_slice": c["L"].to_json(), "exact_sup": c["exact"]}
            if args.sharpen and c["t_bounded"]:
                ctx = CellCtx(A, range(7), g["region2d"], c, "FNFL", verbose=True)
                info = ctx.sharpen()
                if "sup" in info:
                    info["sup"] = info["sup"].to_json()
                e["sharpened"] = info
            cells.append(e)
        report["n7_finite_FNFL_CANDIDATE"] = {"status": "CANDIDATE pending human hand-verification", "aggregate": agg, "cells": cells}
    if args.n7_limit:
        p = os.path.join(REPO, "certificates", "COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl")
        best, best_lab, n_iv, unb = None, [], 0, []
        vals = []
        with open(p) as fh:
            for line in fh:
                if not line.strip():
                    continue
                rec = json.loads(line)
                assert rec.get("complete") is False
                for v in analyse_limit(rec, rec["d6"]):
                    n_iv += 1
                    if v is None:
                        unb.append(rec["d6"])
                        continue
                    vals.append((rec["d6"], v))
        for lab, v in vals:
            if best is None or cmp(v, best) > 0:
                best = v
        best_lab = sorted({lab for lab, v in vals if cmp(v, best) == 0})
        report["n7_limit_CANDIDATES_lower_bound"] = {
            "status": "CANDIDATES pending human hand-verification; LOWER bound on q_7* only",
            "n_graphs_intervals": n_iv, "unbounded": unb,
            "max_beta_c_minus_1": best.to_json() if best else None, "attained_by": best_lab}
    s = json.dumps(report, indent=1, default=str)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(s + "\n")
    print(s)


if __name__ == "__main__":
    main()
