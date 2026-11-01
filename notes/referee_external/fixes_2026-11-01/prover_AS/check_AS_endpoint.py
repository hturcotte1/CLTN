#!/usr/bin/env python3
"""
Exact check of the nu_inf = 1/2 case of the first-order expansion (paper thm:AS / lem:Bmech)
on E_0 = &FNFL~[udo?.

All decisions are exact: sympy polynomials over QQ, arithmetic in QQ[z]/(z^2 - z + 1) = QQ(sqrt(-3)),
fractions.Fraction, exact Hurwitz minors, exact real-root counting. Floats are printed for display only.

Usage:  python3 check_AS_endpoint.py        (exit code 0 iff every check passes)
Writes: check_AS_endpoint_out.json next to this script.
"""
import json, os, sys
from fractions import Fraction as F
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, REPO)
import ctln  # read-only use of the repository verifier (exact)

z, x, y = sp.symbols("z x y")
OUT = {"graph": "&FNFL~[udo?", "checks": {}}
FAIL = []


def check(name, ok, detail=None):
    OUT["checks"][name] = {"ok": bool(ok), "detail": detail}
    print(("PASS " if ok else "FAIL ") + name + ("" if detail is None else ": " + str(detail)))
    if not ok:
        FAIL.append(name)


def decode_d6(s):
    assert s[0] == "&"
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits += [(v >> (5 - k)) & 1 for k in range(6)]
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]


def hurwitz_minors(coeffs):
    """coeffs a0..am (descending), exact; returns the m leading principal Hurwitz minors."""
    m = len(coeffs) - 1
    a = lambda k: coeffs[k] if 0 <= k <= m else 0
    H = sp.Matrix(m, m, lambda i, j: a(2 * (i + 1) - (j + 1)))
    return [sp.Rational(H[:k, :k].det()) for k in range(1, m + 1)]


# ---------------------------------------------------------------- 1. the graph
ROWS = ["0011110", "0011100", "1101111", "1110111", "0011011", "0100101", "1100000"]  # paper, Sec. 5
A = decode_d6("&FNFL~[udo?")
check("digraph6 decodes to the paper's adjacency rows", ["".join(map(str, r)) for r in A] == ROWS)
n = len(A)
check("E_0 is not a clique (complete bidirected)", not all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j))
Am = sp.Matrix(A)
I = sp.eye(n)
one = sp.ones(n, 1)

# ---------------------------------------------------------------- 2. chi, psi
Mz = z * I + Am
chi = sp.Poly(Mz.det(method="berkowitz"), z, domain="QQ")
psi = sp.Poly((one.T * Mz.adjugate(method="berkowitz") * one)[0, 0], z, domain="QQ")
OUT["chi"] = str(chi.as_expr().expand())
OUT["psi"] = str(psi.as_expr().expand())
OUT["psi_factored"] = str(sp.factor(psi.as_expr()))
# transpose invariance (A vs A^T conventions do not matter for chi, psi)
MzT = z * I + Am.T
check("chi, psi invariant under A -> A^T",
      sp.Poly(MzT.det(method="berkowitz"), z, domain="QQ") == chi and
      sp.Poly((one.T * MzT.adjugate(method="berkowitz") * one)[0, 0], z, domain="QQ") == psi)
# cross-check against the paper's compressed-adjacency identity: psi(z) = 7 chi_perp(-z)
lam = sp.symbols("lam")
seven_chi_perp = sp.expand((lam + 1) * (lam**2 + lam + 1) * (7 * lam**3 + 14 * lam**2 + 14 * lam + 5))
check("psi(z) == 7 chi_perp(-z) with the paper's 7 chi_perp factorisation",
      sp.expand(psi.as_expr() - seven_chi_perp.subs(lam, -z)) == 0)

m = sp.Poly(z**2 - z + 1, z, domain="QQ")          # roots z0, z0bar = 1/2 -+ i sqrt(3)/2
qm, rm = sp.div(psi, m)
check("z^2 - z + 1 divides psi", rm.is_zero)
check("z^2 - z + 1 divides psi exactly once (z0 simple)", not sp.div(qm, m)[1].is_zero)
check("z^2 - z + 1 is irreducible over QQ", m.is_irreducible)
g = sp.gcd(chi, m)
check("chi(z0) != 0 (gcd(chi, z^2-z+1) = 1)", g.degree() == 0)
OUT["psi_cofactor"] = str(qm.as_expr())

# every other root of psi has Re > 1/2: cofactor C(z), substitute z = 1/2 - x; C(1/2 - x) Hurwitz in x
C = qm
Cx = sp.Poly(sp.expand(C.as_expr().subs(z, sp.Rational(1, 2) - x)), x)
cx = Cx.all_coeffs()
if cx[0] < 0:
    cx = [-c for c in cx]
hm = hurwitz_minors(cx)
OUT["cofactor_shift_coeffs"] = [str(c) for c in cx]
OUT["cofactor_shift_hurwitz_minors"] = [str(h) for h in hm]
check("all coefficients of +-C(1/2 - x) positive", all(c > 0 for c in cx), [str(c) for c in cx])
check("all Hurwitz minors of +-C(1/2 - x) strictly positive (other roots Re z > 1/2)", all(h > 0 for h in hm),
      [str(h) for h in hm])

# hence alpha = -1/2, nu_inf = 1/2 exactly (roots of psi are -spec(A|perp))
check("roots of z^2-z+1 have real part exactly 1/2", sp.Rational(1, 2) == sp.Rational(-m.all_coeffs()[1], 2 * m.all_coeffs()[0]))

# ---------------------------------------------------------------- 3. R = Re(chi(z0)/psi'(z0)), exactly
# Work in K = QQ[z]/(z^2 - z + 1). An element u + v z with z = z0 = 1/2 - i sqrt(3)/2 has
# Re = u + v/2, Im = -v sqrt(3)/2.   (z0bar gives the conjugate.)
dpsi = psi.diff(z)
chi_r = sp.rem(chi, m)
dpsi_r = sp.rem(dpsi, m)
check("psi'(z0) != 0 in K", not dpsi_r.is_zero)
inv = sp.invert(dpsi_r.as_expr(), m.as_expr(), z)            # exact inverse in K
quot = sp.Poly(sp.rem(sp.Poly(sp.expand(chi_r.as_expr() * inv), z), m), z, domain="QQ")
coef = quot.all_coeffs()
coef = [sp.Rational(0)] * (2 - len(coef)) + coef
v, u = coef[0], coef[1]
Re = u + v / 2
Im_over_sqrt3 = -v / 2
OUT["chi_mod_m"] = str(chi_r.as_expr()); OUT["dpsi_mod_m"] = str(dpsi_r.as_expr())
OUT["chi_over_dpsi_at_z0"] = {"as_u_plus_v_z0": [str(u), str(v)], "Re": str(Re), "Im/sqrt(3)": str(Im_over_sqrt3)}
check("chi(z0)/psi'(z0) = 1/4 - (sqrt(3)/12) i", Re == sp.Rational(1, 4) and Im_over_sqrt3 == -sp.Rational(1, 12),
      f"Re={Re}, Im={Im_over_sqrt3}*sqrt(3)")
# independent cross-check with sympy radicals
z0 = sp.Rational(1, 2) - sp.sqrt(3) / 2 * sp.I
val = sp.nsimplify(sp.radsimp(sp.expand(chi.as_expr().subs(z, z0)) / sp.expand(dpsi.as_expr().subs(z, z0))))
check("radical cross-check of chi(z0)/psi'(z0)", sp.simplify(val - (sp.Rational(1, 4) - sp.sqrt(3) / 12 * sp.I)) == 0, str(val))
# the same via 1/s^'(z0), s^ = psi/chi
s = psi.as_expr() / chi.as_expr()
val2 = sp.radsimp(sp.simplify(1 / sp.diff(s, z).subs(z, z0)))
check("1/s^'(z0) equals chi(z0)/psi'(z0)", sp.simplify(val2 - val) == 0, str(val2))
R = F(int(sp.fraction(Re)[0]), int(sp.fraction(Re)[1]))
check("R = 1/4 > 1/8 (exact Fraction comparison)", R > F(1, 8), str(R))
OUT["R"] = str(R)

# ---------------------------------------------------------------- 4. legality boundary expansion
# nu_P2(c) = c - sqrt(c^2 - c) = 1/(1 + sqrt(1 - x)),  x = 1/c in (0, 1).
ser = sp.series(1 / (1 + sp.sqrt(1 - x)), x, 0, 5).removeO()
OUT["nu_P2_series_in_x=1/c"] = str(sp.expand(ser))
check("nu_P2 = 1/2 + x/8 + x^2/16 + 5x^3/128 + 7x^4/256 + O(x^5)",
      sp.expand(ser - (sp.Rational(1, 2) + x / 8 + x**2 / 16 + 5 * x**3 / 128 + 7 * x**4 / 256)) == 0)
check("identity c - sqrt(c^2-c) == 1/(1+sqrt(1-1/c)) (c>1)",
      sp.simplify((sp.Symbol('c', positive=True) - sp.sqrt(sp.Symbol('c', positive=True)**2 - sp.Symbol('c', positive=True)))
                  - 1 / (1 + sp.sqrt(1 - 1 / sp.Symbol('c', positive=True)))) == 0)
# rigorous two-sided bound on (0, 1/2]:  1/2 + x/8 < nu_P2 < 1/2 + x/8 + x^2/8.
# With s = sqrt(1-x) > 0:  1/(1+s) > L  <=>  s < 1/L - 1,  and  1/(1+s) < U  <=>  s > 1/U - 1.
L = sp.Rational(1, 2) + x / 8
U = sp.Rational(1, 2) + x / 8 + x**2 / 8
# lower: need 1/L - 1 > 0 and (1-x) < (1/L - 1)^2 on (0,1); numerator polynomial:
low = sp.Poly(sp.numer(sp.together((1 / L - 1)**2 - (1 - x))), x)
upp = sp.Poly(sp.numer(sp.together((1 - x) - (1 / U - 1)**2)), x)
den_low = sp.Poly(sp.denom(sp.together((1 / L - 1)**2 - (1 - x))), x)
den_upp = sp.Poly(sp.denom(sp.together((1 - x) - (1 / U - 1)**2)), x)
OUT["lower_bound_poly"] = str(low.as_expr()); OUT["upper_bound_poly"] = str(upp.as_expr())
half = sp.Rational(1, 2)
def positive_on_open_interval(p, a, b):
    # no root in [a,b] except possibly x = 0 handled by factoring; and positive at midpoint
    q = p
    while q.eval(0) == 0:
        q = sp.Poly(sp.quo(q, sp.Poly(x, x)), x)
    return q.count_roots(a, b) == 0 and q.eval((a + b) / 2) > 0, str(q.as_expr())
okl, ql = positive_on_open_interval(low * den_low, sp.Rational(0), sp.Rational(1))
def no_root_and_positive(expr, a, b):
    # expr rational in x; numerator*denominator has no root in [a,b] and is positive at a point inside
    nd = sp.Poly(sp.numer(sp.together(expr)), x) * sp.Poly(sp.denom(sp.together(expr)), x)
    return nd.count_roots(a, b) == 0 and nd.eval((a + b) / 2) > 0
check("1/2 + 1/(8c) < c - sqrt(c^2-c) for all c > 1 (exact root count on (0,1))",
      okl and no_root_and_positive(1 / L - 1, sp.Rational(0), sp.Rational(1)), ql)
oku, qu = positive_on_open_interval(upp * den_upp, sp.Rational(0), half)
check("c - sqrt(c^2-c) < 1/2 + 1/(8c) + 1/(8c^2) for all c >= 2 (exact root count on (0,1/2])",
      oku and no_root_and_positive(1 / U - 1, sp.Rational(0), half), qu)

# ---------------------------------------------------------------- 5. finite-c sanity check (exact)
# At c = 10^k choose gamma = 1/2 + 3/(16c) (midway between the first-order legality and stability coefficients
# 1/8 and R = 1/4) -- a rational point; then beta = 1/(1-gamma), t = beta(c-1), eps = 1/(1+t),
# delta = (beta-1)/(1+t) and run the repository verifier ctln.analyze_support (exact) on the full support.
fin = []
for c in [F(10), F(100), F(1000), F(10000)]:
    gam = F(1, 2) + F(3, 16) / c
    beta = 1 / (1 - gam)
    t = beta * (c - 1)
    eps = 1 / (1 + t)
    delta = (beta - 1) / (1 + t)
    rec = ctln.analyze_support(A, list(range(n)), eps, delta)
    leg = ctln.legal(eps, delta)
    entry = {"c": str(c), "gamma": str(gam), "beta": str(beta), "t": str(t), "eps": str(eps), "delta": str(delta),
             "q": rec["q"], "legal": leg, "positive": rec["positive"], "stable": rec["stable"],
             "stability_detail": rec["stability_detail"], "stable_fixed_point": rec["stable_fixed_point"],
             "q_display_only": float(F(rec["q"]))}
    fin.append(entry)
    check(f"c={c}: legal, positive, Hurwitz-stable full support (SFP) at gamma=1/2+3/(16c)",
          leg and rec["stable_fixed_point"] is True and rec["stable"] is True,
          f"eps={eps}, delta={delta}, {rec['stability_detail']}")
OUT["finite_c_witnesses"] = fin

# Upper side: at gamma' = 1/2 + 5/(16c) (first-order coefficient 5/16 > R = 1/4) the point is legal, and the
# full support should be Hurwitz-UNSTABLE for large c (stability of K does not depend on positivity; at c = 10,
# beta = 32/15 lies beyond the positivity cell (2, beta_1'), beta_1' ~ 2.083, so 'positive' is False there). Together with the witnesses above this brackets
# c (nu_A(c) - 1/2) in (3/16, 5/16) at each tested c, exactly (a three-valued result; None would be reported).
upper = []
for c in [F(10), F(100), F(1000), F(10000)]:
    gam = F(1, 2) + F(5, 16) / c
    beta = 1 / (1 - gam); t = beta * (c - 1); eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
    rec = ctln.analyze_support(A, list(range(n)), eps, delta)
    upper.append({"c": str(c), "gamma": str(gam), "eps": str(eps), "delta": str(delta), "legal": ctln.legal(eps, delta),
                  "positive": rec["positive"], "stable": rec["stable"], "stability_detail": rec["stability_detail"]})
    check(f"c={c}: legal, Hurwitz-UNSTABLE (decided, not undecided) at gamma=1/2+5/(16c)",
          ctln.legal(eps, delta) and rec["stable"] is False,
          f"positive={rec['positive']}, {rec['stability_detail']}")
OUT["finite_c_upper_side"] = upper

# Consistency with the recorded E_0 window (paper Sec. 5: at beta = 53/26 the stable range is 26 < t < t_1, t_1 in (53,54)):
# the first-order curve beta = 2 + 8R/t of thm:AS at nu_inf = 1/2 predicts t_1 ~ 8R/(beta - 2), and legality is t > 1/(beta-2).
b = F(53, 26)
OUT["first_order_prediction_at_beta_53_26"] = {"8R/(beta-2)": str(8 * R / (b - 2)), "1/(beta-2)": str(1 / (b - 2))}
check("first-order prediction 8R/(beta-2) at beta=53/26 equals 52 (paper: t_1 in (53,54); legality t > 26)",
      8 * R / (b - 2) == 52 and 1 / (b - 2) == 26)

OUT["all_pass"] = not FAIL
OUT["failed"] = FAIL
with open(os.path.join(HERE, "check_AS_endpoint_out.json"), "w") as f:
    json.dump(OUT, f, indent=1)
print("ALL PASS" if not FAIL else f"FAILURES: {FAIL}")
sys.exit(0 if not FAIL else 1)
