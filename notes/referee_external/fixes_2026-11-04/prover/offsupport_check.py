"""
offsupport_check.py -- exact reproduction of the example in J. Geneson's second
correctness report (2026-11-04), very minor point 1: the opening boundary list
of the two-parameter paragraph in Section 3.1 omits the off-support curves
o_k = 0 (k not in sigma) for proper supports.

All decisions exact (sympy integers/rationals, fractions.Fraction). Floats are
only printed as labelled approximations.

Definitions re-implemented from Section 2 of notes/paper/paper.tex:
    N = J - I - A_sigma^T   (N_ab = 1 iff a != b and sigma_b -/-> sigma_a)
    C = I + beta N,  K = C + t J,  a = adj(C) 1,  d = det C,  s = 1^T a,
    D = det K = d + t s,  p(mu) = det(mu I + K),
    o_k = d - beta * sum_b N_kb a_b,  N_kb = 1 iff sigma_b -/-> k   (Lemma 2.1(v)),
    L = t (beta - 2) - 1.
The off-support drive is ALSO derived directly from Definition 2.2(2) with the
weight matrix (2.1) and checked to equal theta * o_k / D.

Graph: G_1 of Table 1 (rows 011000 000011 000111 110011 110101 101110) on
vertices 1..6, plus vertex 7 receiving arcs from exactly 2,3,5,6 and sending
none. The stage-0 labelling of Geneson's graph is also run for comparison.

Run: cd /home/user/CLTN && nice -n 10 python3 notes/referee_external/fixes_2026-11-04/prover/offsupport_check.py
Writes offsupport_check.json next to this script.
"""
import json
import os
import sys
from fractions import Fraction as F

import sympy as sp

sys.set_int_max_str_digits(0)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, ROOT)
import ctln  # noqa: E402  (cross-check only)

b, t, mu, th = sp.symbols("beta t mu theta", positive=True)
LO, HI = sp.Rational(52, 25), sp.Rational(209, 100)
T0 = sp.Integer(100)

G1_ROWS = ["011000", "000011", "000111", "110011", "110101", "101110"]
G1 = [[int(c) for c in r] for r in G1_ROWS]
STAGE0 = json.load(open(os.path.join(ROOT, "certificates", "stage0_geneson_2607.21396.json")))["A_ours"]
RECV_FROM_1BASED = [2, 3, 5, 6]


def extend(A6):
    """7-node digraph: vertex 7 (index 6) receives from RECV_FROM, sends none.
    Convention A[i][j] = 1 iff arc i -> j."""
    A = [row[:] + [0] for row in A6] + [[0] * 7]
    for v in RECV_FROM_1BASED:
        A[v - 1][6] = 1
    return A


def hurwitz(coeffs):
    """coeffs = [c0, c1, ..., cm], c0 > 0. Returns [Delta_1..Delta_m]."""
    m = len(coeffs) - 1
    c = lambda l: coeffs[l] if 0 <= l <= m else 0
    return [sp.Matrix(k, k, lambda i, j: c(2 * (i + 1) - (j + 1))).det(method="berkowitz")
            for k in range(1, m + 1)]


def own_sturm_count(poly_expr, lo, hi):
    """Number of distinct real roots in (lo, hi] by an own Sturm chain over Fraction."""
    P = sp.Poly(sp.expand(poly_expr), b)
    seq = [P, P.diff(b)]
    while not seq[-1].is_zero and seq[-1].degree() > 0:
        r = -seq[-2].rem(seq[-1])
        if r.is_zero:
            break
        seq.append(r)
    def var(x):
        vals = [(lambda r: F(int(r.p), int(r.q)))(sp.Rational(q.eval(x))) for q in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))
    return var(lo) - var(hi)


def analyse(A6, label):
    A7 = extend(A6)
    sigma = list(range(6))
    m = 6
    As = sp.Matrix(6, 6, lambda i, j: A6[i][j])
    J = sp.ones(6, 6)
    N = J - sp.eye(6) - As.T
    C = sp.eye(6) + b * N
    K = C + t * J
    adjC = C.adjugate()
    a = [sp.expand(e) for e in (adjC * sp.ones(6, 1))]
    d = sp.expand(C.det(method="berkowitz"))
    s = sp.expand(sum(a))
    D = sp.expand(K.det(method="berkowitz"))
    assert sp.expand(D - (d + t * s)) == 0
    # o_7 from Lemma 2.1(v)
    k = 6
    Nk = [0 if A7[bb][k] == 1 else 1 for bb in sigma]       # N_kb = 1 iff sigma_b -/-> k
    o7 = sp.expand(d - b * sum(Nk[bb] * a[bb] for bb in sigma))
    # direct derivation from Definition 2.2 with W of (2.1), eps = 1/(1+t), delta = (beta-1)/(1+t)
    eps = 1 / (1 + t)
    delta = (b - 1) / (1 + t)
    W = sp.Matrix(7, 7, lambda i, j: 0 if i == j else (-1 + eps if A7[j][i] == 1 else -1 - delta))
    IWs = sp.eye(6) - W[:6, :6]
    # check closed form (2.2) and IminusW = eps K
    closed = -delta * sp.eye(6) + (1 + delta) * J - (eps + delta) * As.T
    assert sp.simplify(IWs - closed) == sp.zeros(6, 6)
    assert sp.simplify(IWs - eps * K) == sp.zeros(6, 6)
    x = (th / eps) * sp.Matrix(a) / D
    assert sp.simplify(IWs * x - th * sp.ones(6, 1)) == sp.zeros(6, 1)
    drive = th + sum(W[k, j] * x[j] for j in sigma)
    assert sp.simplify(drive - th * o7 / D) == 0
    # characteristic polynomial of K side: p(mu) = det(mu I + K); eigenvalues of -I+W_sigma are -eps*spec(K)
    p = sp.Poly(sp.expand((mu * sp.eye(6) + K).det(method="berkowitz")), mu)
    coeffs = [sp.expand(c) for c in p.all_coeffs()]
    assert coeffs[0] == 1
    coeffs100 = [sp.expand(c.subs(t, T0)) for c in coeffs]
    H = [sp.expand(h) for h in hurwitz(coeffs100)]
    L = sp.expand(T0 * (b - 2) - 1)
    D100 = sp.expand(D.subs(t, T0))

    fam = {f"a_{i+1}": a[i] for i in range(6)}
    fam["D(t=100)"] = D100
    for i, h in enumerate(H):
        fam[f"Delta_{i+1}(t=100)"] = h
    fam["L(t=100)"] = L
    fam["o_7"] = o7
    rows = {}
    for name, f in fam.items():
        P = sp.Poly(f, b)
        cnt_closed = P.count_roots(LO, HI)          # sympy exact, closed interval
        cnt_own = own_sturm_count(f, LO, HI)        # own Sturm, (LO, HI]
        vlo, vhi = f.subs(b, LO), f.subs(b, HI)
        rows[name] = {"degree": P.degree(), "roots_in_closed_segment_sympy": int(cnt_closed),
                      "roots_in_(lo,hi]_own_sturm": int(cnt_own),
                      "value_at_52/25": str(vlo), "value_at_209/100": str(vhi),
                      "sign_lo": int(sp.sign(vlo)), "sign_hi": int(sp.sign(vhi))}
    target = -b**5 + 3*b**4 - 3*b**3 + 3*b**2 - 2*b + 1
    res = {
        "label": label,
        "A7": ["".join(map(str, r)) for r in A7],
        "N_7b_(b=1..6)": Nk,
        "o_7": str(o7),
        "o_7_equals_report_polynomial": sp.expand(o7 - target) == 0,
        "o_7_over_report_polynomial": str(sp.cancel(o7 / target)),
        "d": str(d), "a": [str(e) for e in a],
        "family": rows,
    }
    # endpoint verdicts from Lemma 2.3 (exact)
    def verdict(beta0):
        pos = all(e.subs(b, beta0) > 0 for e in a)
        off = o7.subs(b, beta0) < 0
        hur = all(h.subs(b, beta0) > 0 for h in H)
        return {"positivity": bool(pos), "off_support_o7<0": bool(off),
                "hurwitz_all_positive": bool(hur), "legal": bool(L.subs(b, beta0) > 0),
                "SFP": bool(pos and off and hur)}
    res["lemma_verdict_52/25"] = verdict(LO)
    res["lemma_verdict_209/100"] = verdict(HI)
    # cross-check with ctln.analyze_support (project convention, independent code path)
    xc = {}
    for nm, beta0 in [("52/25", F(52, 25)), ("209/100", F(209, 100))]:
        e = F(1, 101)
        dl = (beta0 - 1) / 101
        rec = ctln.analyze_support(A7, sigma, e, dl)
        # outside drive recomputed from rec x for the report
        Wc = ctln.weight_matrix(A7, e, dl)
        xs = [F(v) for v in rec["x"]]
        drv = 1 + sum(Wc[6][j] * xs[j] for j in range(6))
        Dr = sp.Rational(D.subs({t: 100, b: sp.Rational(beta0.numerator, beta0.denominator)}))
        orr = sp.Rational(o7.subs(b, sp.Rational(beta0.numerator, beta0.denominator)))
        D_val, o_val = F(int(Dr.p), int(Dr.q)), F(int(orr.p), int(orr.q))
        xc[nm] = {"eps": str(e), "delta": str(dl), "q": rec["q"],
                  "positive": rec["positive"], "off_support_ok": rec["off_support_ok"],
                  "stable": rec["stable"], "BOUNDARY": rec["BOUNDARY"],
                  "stable_fixed_point": rec["stable_fixed_point"],
                  "outside_drive": str(drv), "drive_equals_o7_over_D": drv == o_val / D_val}
    res["ctln_cross_check"] = xc
    return res


def main():
    out = {"script": os.path.relpath(__file__, ROOT),
           "segment": [str(LO), str(HI)], "t": 100,
           "G1_table1": analyse(G1, "G_1 of Table 1 (canonical labelling)"),
           "stage0_labelling": analyse(STAGE0, "Geneson graph, stage-0 certificate labelling")}
    g = out["G1_table1"]
    listed = [k for k in g["family"] if k != "o_7"]
    out["summary"] = {
        "o7_matches_report": g["o_7_equals_report_polynomial"],
        "o7(52/25)": g["family"]["o_7"]["value_at_52/25"],
        "o7(209/100)": g["family"]["o_7"]["value_at_209/100"],
        "listed_family_root_free_and_positive_on_segment": all(
            g["family"][k]["roots_in_closed_segment_sympy"] == 0 and g["family"][k]["roots_in_(lo,hi]_own_sturm"] == 0
            and g["family"][k]["sign_lo"] == 1 and g["family"][k]["sign_hi"] == 1 for k in listed),
        "o7_roots_in_segment": g["family"]["o_7"]["roots_in_closed_segment_sympy"],
        "SFP_at_52/25": g["lemma_verdict_52/25"]["SFP"],
        "SFP_at_209/100": g["lemma_verdict_209/100"]["SFP"],
        "ctln_SFP_at_52/25": g["ctln_cross_check"]["52/25"]["stable_fixed_point"],
        "ctln_SFP_at_209/100": g["ctln_cross_check"]["209/100"]["stable_fixed_point"],
    }
    json.dump(out, open(os.path.join(HERE, "offsupport_check.json"), "w"), indent=1)
    print(json.dumps(out["summary"], indent=1))
    for lab in ("G1_table1", "stage0_labelling"):
        r = out[lab]
        print(lab, "o_7 =", r["o_7"], "| matches report:", r["o_7_equals_report_polynomial"],
              "| ratio:", r["o_7_over_report_polynomial"])
        for k, v in r["family"].items():
            print(f"  {k:16s} deg {v['degree']:2d} roots[seg] {v['roots_in_closed_segment_sympy']} own {v['roots_in_(lo,hi]_own_sturm']}"
                  f"  signs {v['sign_lo']:+d},{v['sign_hi']:+d}")
        print("  verdicts:", r["lemma_verdict_52/25"], r["lemma_verdict_209/100"])
        print("  ctln:", r["ctln_cross_check"])


if __name__ == "__main__":
    main()
