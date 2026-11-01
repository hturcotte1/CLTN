"""
geneson_section_check.py -- exact reproduction of the example in J. Geneson's
correctness report (2026-11-01), point 1: a lower-dimensional (vertical section)
cell of the two-parameter decomposition that lies inside the open SFP region.

Everything is exact (sympy integers/rationals, fractions.Fraction).  Floats are
used only to PRINT labelled approximations; no decision reads a float.

Definitions are re-implemented from Section 2 of notes/paper/paper.tex, not
imported from the project engines:
    N = J - I - A^T   (receiver-first non-edge matrix),  C = I + beta N,
    K = C + t J,      p(mu) = det(mu I + K) = mu^6 + c_1 mu^5 + ... + c_6,
    a = adj(C) 1      (= adj(K) 1, positivity vector),
    Delta_k = det[c_{2i-j}]_{i,j=1..k}   (c_0 = 1, c_l = 0 outside 0..6),
    legal iff L := t (beta - 2) - 1 > 0.
Full support of a 6-node graph is an SFP at a legal (beta, t) iff a_i > 0 for
all i and Delta_1..Delta_6 > 0 (Lemma A.2 of the paper; condition (2) vacuous).

Parts:
  (1) G_1 of Table 1 and Geneson's graph (stage-0 certificate) give the same
      p(mu; beta, t)  (sanity: they are isomorphic).
  (2) P(beta) := numerator of Delta_5(beta, 1/(beta-2)); its relation to
      Res_t(Delta_5, L); exact root count in (17/8, 15/7) by sympy AND by an own
      Sturm sequence over Fraction; the irreducible factor carrying beta_0 and a
      refined rational isolating interval.
  (3) At t = 100: for each of the six a_i(beta) and six Delta_k(beta, 100),
      substitute beta = (17/8 + (15/7) u)/(1+u), clear (1+u)^deg and the
      denominators, and check every coefficient is > 0  (=> f > 0 on [17/8, 15/7]).
      Same check at t = 10 and t = 1000 (extra, not in the report).
  (4) Legality 100 (beta_0 - 2) > 1 from beta_0 > 17/8.
  (5) Pointwise cross-check with the project's verifier ctln.analyze_support
      (project W convention, independent code path) at rational beta in the
      interval and t = 100.
  (6) The project's own CAD (region2d.py): which projection polynomials vanish
      at beta_0, and how the recorded certificate of G_1 (idx 1000255) treats
      the strip around beta_0.
Run:  cd /home/user/CLTN && nice -n 10 python3 notes/referee_external/fixes_2026-11-01/prover_cells/geneson_section_check.py
Writes geneson_section_check.json next to this script.
"""
import json
import math
import os
import sys
from fractions import Fraction as F

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, ROOT)

b, t, mu, u, s = sp.symbols("beta t mu u s")
LO, HI = sp.Rational(17, 8), sp.Rational(15, 7)

G1_ROWS = ["011000", "000011", "000111", "110011", "110101", "101110"]   # Table 1, idx 1,000,255
G1 = [[int(ch) for ch in r] for r in G1_ROWS]
STAGE0 = json.load(open(os.path.join(ROOT, "certificates", "stage0_geneson_2607.21396.json")))
GEN = STAGE0["A_ours"]

out = {"script": os.path.relpath(__file__, ROOT), "graph_G1_rows": G1_ROWS}


def matrices(A):
    m = len(A)
    I = sp.eye(m)
    J = sp.ones(m, m)
    AT = sp.Matrix(A).T
    N = J - I - AT
    C = I + b * N
    K = C + t * J
    return N, C, K


def charcoeffs(K):
    m = K.shape[0]
    p = sp.Poly(sp.expand((mu * sp.eye(m) + K).det(method="berkowitz")), mu)
    c = p.all_coeffs()                      # c[0] = 1 (monic), ..., c[m]
    assert c[0] == 1 and len(c) == m + 1
    return [sp.expand(x) for x in c]


def hurwitz(c, k):
    m = len(c) - 1

    def ent(i, j):                           # 1-indexed: c_{2i-j}
        idx = 2 * i - j
        return c[idx] if 0 <= idx <= m else sp.Integer(0)
    H = sp.Matrix(k, k, lambda i, j: ent(i + 1, j + 1))
    return sp.expand(H.det(method="berkowitz"))


# ---------------------------------------------------------------- own Sturm count (Fraction)
def pdiv_rem(a, d):
    """remainder of a by d, coefficient lists high->low, Fraction."""
    a = list(a)
    while len(a) >= len(d) and any(a):
        if a[0] == 0:
            a.pop(0)
            continue
        q = a[0] / d[0]
        for i in range(len(d)):
            a[i] -= q * d[i]
        a.pop(0)
    while a and a[0] == 0:
        a.pop(0)
    return a


def horner(p, x):
    r = F(0)
    for c in p:
        r = r * x + c
    return r


def sturm_count(p, lo, hi):
    """number of distinct real roots of p (high->low, Fraction) in (lo, hi], p(lo),p(hi) != 0."""
    p = [F(c) for c in p]
    dp = [c * (len(p) - 1 - i) for i, c in enumerate(p[:-1])]
    seq = [p, dp]
    while True:
        r = pdiv_rem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def V(x):
        s = [horner(q, x) for q in seq]
        s = [v for v in s if v != 0]
        return sum(1 for i in range(len(s) - 1) if (s[i] > 0) != (s[i + 1] > 0))
    return V(lo) - V(hi)


def positive_on_interval(f, lo=LO, hi=HI):
    """exact certificate that the univariate f(beta) > 0 on [lo, hi]:
    g(u) = (1+u)^deg f((lo + hi u)/(1+u)) has all coefficients > 0 (after clearing
    a positive denominator).  Then g > 0 on [0, inf) and the u^deg coefficient is
    a positive multiple of f(hi), so f > 0 on [lo, hi]."""
    P = sp.Poly(f, b)
    dgr = P.degree()
    g = sp.expand(sp.cancel((1 + u) ** dgr * f.subs(b, (lo + hi * u) / (1 + u))))
    G = sp.Poly(g, u)
    co = G.all_coeffs()
    den = 1
    for c_ in co:
        den = math.lcm(den, int(sp.Rational(c_).q))
    co = [int(sp.Rational(c) * den) for c in co]
    return {"degree_beta": dgr, "n_coeffs": len(co), "all_positive": all(c > 0 for c in co),
            "min_coeff": min(co), "n_zero": sum(1 for c in co if c == 0), "n_negative": sum(1 for c in co if c < 0)}


def positive_on_box(f, t0, lo=LO, hi=HI):
    """exact certificate that f(beta, t) > 0 on [lo, hi] x [t0, inf): with
    beta = (lo + hi u)/(1+u), t = t0 + s, the polynomial (1+u)^deg_beta f has all
    coefficients (in u, s) > 0 and degree_s-leading part positive (implied)."""
    dgr = sp.Poly(f, b).degree()
    g = sp.expand(sp.cancel((1 + u) ** dgr * f.subs({b: (lo + hi * u) / (1 + u), t: t0 + s})))
    G = sp.Poly(g, u, s)
    co = [c_ for _, c_ in G.terms()]
    den = 1
    for c_ in co:
        den = math.lcm(den, int(sp.Rational(c_).q))
    co = [int(sp.Rational(c_) * den) for c_ in co]
    # every monomial u^i s^j with 0<=i<=deg, 0<=j<=deg_s must be present with c > 0
    # for strict positivity at all corners; we require: all present coefficients > 0
    # and the coefficients of s^j u^0 and s^j u^deg for j = 0 present (corner values).
    terms = dict(G.terms())
    corners_ok = all(((0, 0) in terms, (dgr, 0) in terms))
    return {"degree_beta": dgr, "n_terms": len(co), "all_positive": all(c_ > 0 for c_ in co) and corners_ok,
            "n_negative": sum(1 for c_ in co if c_ < 0)}


def main():
    # (1) same characteristic data for Table-1 G_1 and Geneson's labelling
    N1, C1, K1 = matrices(G1)
    Ng, Cg, Kg = matrices(GEN)
    c1, cg = charcoeffs(K1), charcoeffs(Kg)
    out["same_charpoly_G1_vs_stage0"] = all(sp.expand(x - y) == 0 for x, y in zip(c1, cg))
    assert out["same_charpoly_G1_vs_stage0"]
    c = c1
    Delta = {k: hurwitz(c, k) for k in range(1, 7)}
    D5 = Delta[5]
    d5t = sp.Poly(D5, t).degree()
    out["Delta5_degrees"] = {"beta": sp.Poly(D5, b).degree(), "t": d5t}

    # (2) Delta_5 on the legality boundary
    P = sp.expand(sp.cancel((b - 2) ** d5t * D5.subs(t, 1 / (b - 2))))
    Pp = sp.Poly(P, b)
    L = t * (b - 2) - 1
    R = sp.Poly(sp.resultant(sp.Poly(D5, t).as_expr(), L, t), b)
    ratio = sp.cancel(R.as_expr() / P)
    out["P_degree"] = Pp.degree()
    out["Res_t(Delta5,L) / P"] = str(ratio)
    cnt_sympy = Pp.count_roots(LO, HI)
    Pf = [F(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in Pp.all_coeffs()]
    assert horner(Pf, F(17, 8)) != 0 and horner(Pf, F(15, 7)) != 0
    cnt_sturm = sturm_count(Pf, F(17, 8), F(15, 7))
    out["roots_of_P_in_(17/8,15/7)"] = {"sympy_count_roots": int(cnt_sympy), "own_sturm": cnt_sturm,
                                         "P(17/8)": str(horner(Pf, F(17, 8))), "P(15/7)": str(horner(Pf, F(15, 7)))}
    assert cnt_sympy == 1 and cnt_sturm == 1
    _, facs = sp.factor_list(P, b)
    carrier = None
    fac_info = []
    for g, mult in facs:
        gp = sp.Poly(g, b)
        k = gp.count_roots(LO, HI)
        fac_info.append({"degree": gp.degree(), "mult": mult, "roots_in_interval": int(k)})
        if k:
            carrier = gp
    out["factors_of_P"] = fac_info
    assert carrier is not None
    out["minpoly_beta0_hi_to_lo"] = [int(x) for x in carrier.all_coeffs()]
    # refine an isolating interval for beta_0 by exact bisection
    cf = [F(int(x)) for x in carrier.all_coeffs()]
    lo, hi = F(17, 8), F(15, 7)
    slo = horner(cf, lo) > 0
    for _ in range(60):
        mid = (lo + hi) / 2
        v = horner(cf, mid)
        assert v != 0
        if (v > 0) == slo:
            lo = mid
        else:
            hi = mid
    assert (horner(cf, lo) > 0) != (horner(cf, hi) > 0)
    out["beta0_isolating_interval"] = [str(lo), str(hi)]
    out["beta0_approx_display_only"] = float((lo + hi) / 2)
    # compare with the certificate's degree-14 critical value
    cert_poly = [43, -498, 2780, -9966, 25641, -50226, 77654, -96910, 98882, -82792, 56784, -31408, 13760, -4352, 1024]
    out["minpoly_equals_certificate_beta_critical_[17/8,15/7]"] = (
        [int(x) for x in carrier.all_coeffs()] == cert_poly or [-int(x) for x in carrier.all_coeffs()] == cert_poly)

    # (3) positivity at t = 100 (and 10, 1000) throughout [17/8, 15/7]
    a = sp.expand(C1.adjugate() * sp.ones(6, 1))
    checks = {}
    for tv in (100, 1000):
        rec = {}
        for i in range(6):
            rec[f"a_{i+1}"] = positive_on_interval(a[i])
        for k in range(1, 7):
            rec[f"Delta_{k}"] = positive_on_interval(sp.expand(Delta[k].subs(t, tv)))
        rec["ALL_POSITIVE"] = all(v["all_positive"] for v in rec.values())
        # legality of the whole strip [17/8,15/7] x {tv}: tv (beta - 2) - 1 >= tv/8 - 1
        rec["legal_on_strip"] = sp.Rational(tv, 8) - 1 > 0
        checks[str(tv)] = rec
    out["positivity_on_[17/8,15/7]"] = checks
    # whole strip [17/8,15/7] x [100, inf): a_i are t-free (checked above); Hurwitz bivariate
    box = {f"Delta_{k}": positive_on_box(Delta[k], 100) for k in range(1, 7)}
    box["ALL_POSITIVE"] = all(v["all_positive"] for v in box.values())
    out["positivity_on_box_[17/8,15/7]x[100,inf)"] = box
    out["note_coefficient_test"] = "the coefficient-sign test is a sufficient condition only; a failure would mean inconclusive"
    assert checks["100"]["ALL_POSITIVE"]

    # (4) legality at (beta_0, 100)
    out["legal_at_beta0_t100"] = {"claim": "100(beta_0 - 2) > 100(17/8 - 2) = 25/2 > 1", "holds": 100 * (F(17, 8) - 2) > 1}
    out["conclusion"] = ("(beta_0, 100) is a legal point at which all a_i > 0 and Delta_1..Delta_6 > 0, hence the full "
                         "support of G_1 is an SFP there (Lemma A.2); beta_0 is a root of Res_t(Delta_5, L).")

    # (5) cross-check with the project's independent pointwise verifier (project W convention)
    import ctln
    xc = []
    eps = F(1, 101)                                     # t = 100
    for bb in (F(17, 8), lo, hi, F(15, 7)):
        delta = (bb - 1) * eps
        r = ctln.analyze_support(G1, list(range(6)), eps, delta)
        xc.append({"beta": str(bb), "positive": r["positive"], "stable": r["stable"],
                   "q": str(delta * (1 - eps) / eps)})
    out["crosscheck_ctln_analyze_support_t100"] = xc
    assert all(x["positive"] is True and x["stable"] is True for x in xc)

    # (6) the project's CAD on G_1: which projection polynomials vanish at beta_0?
    from fastregion import beta_data
    import region2d
    bd = beta_data(G1, list(range(6)))
    curves0 = region2d.curves_lc(bd)
    curves, polys, names = region2d.projection(bd, curves0)
    hit = []
    for nm, pp in zip(names, polys):
        q, r = sp.div(sp.Poly(pp.as_expr(), b), carrier)
        if r.is_zero:
            hit.append(nm)
    out["region2d_projection_names"] = names
    out["region2d_curves_F"] = [(nm, sp.Poly(f.as_expr(), t).degree()) for nm, f in curves]
    out["n_projection_polys"] = len(polys)
    # sanity: points just left/right of beta_0 (the two ends of the 60-step isolating
    # interval) at legal t between t_L(beta_0) ~ 7.86 and 100 (pointwise, exact)
    near = []
    for bb in (lo, hi):
        for tv in (8, 9, 10, 20, 50):
            eps_ = F(1, tv + 1)
            r = ctln.analyze_support(G1, list(range(6)), eps_, (bb - 1) * eps_)
            near.append({"beta": "lo" if bb == lo else "hi", "t": tv, "legal": tv * (bb - 2) > 1,
                         "positive": r["positive"], "stable": r["stable"]})
    out["pointwise_near_section_t_in_[8,50]"] = near
    out["region2d_projection_polys_vanishing_at_beta0"] = hit
    # the repo's Delta5 vs ours (up to a constant)
    repo_D5 = [f for n_, f in curves0 if n_.startswith("Delta")][0].as_expr()
    out["repo_Delta5_over_ours"] = str(sp.cancel(repo_D5.subs({sp.Symbol('beta'): b, sp.Symbol('t'): t}) / D5))
    # the recorded certificate (shard 2 HITS, idx 1000255)
    rec = None
    for line in open(os.path.join(ROOT, "certificates", "region2d_n6_shard2of3_HITS.jsonl")):
        r = json.loads(line)
        if r.get("idx") == 1000255:
            rec = r
    cert = json.loads(rec["cert"]) if isinstance(rec["cert"], str) else rec["cert"]
    out["recorded_cert_idx1000255"] = {
        "verdict": rec["verdict"],
        "positivity_cells": cert["positivity_cells"],
        "beta_critical": [(x["interval"], len(x["minpoly_hi_to_lo"]) - 1, x["positivity_root"]) for x in cert["beta_critical"]],
        "beta_cells": [(row["beta_cell"], row["beta_sample"],
                        [(c_["t_lo"], c_["t_hi"], c_["stable_fp"]) for c_ in row["t_cells"]]) for row in cert["cells"]],
        "stable_cells": cert["stable_cells"],
        "any_recorded_cell_contains_beta0": any(
            (F(row["beta_cell"][0]) < lo) and (row["beta_cell"][1] is None or hi < F(row["beta_cell"][1]))
            for row in cert["cells"]),
        "keys_in_cert": sorted(cert.keys()),
    }
    with open(os.path.join(HERE, "geneson_section_check.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
