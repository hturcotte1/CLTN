"""Independent exact check of the 4 November report example (Section 3.1 boundary list).

Built from CLAUDE.md section 2 only (A[i][j]=1 iff i->j; W[i][j] = -1+eps if A[j][i] else -1-delta),
no project module in the decision path. Exact arithmetic throughout (Fraction, sympy over QQ).
"""
import sys, json
from fractions import Fraction as F
import sympy as sp

b, t, mu = sp.symbols("beta t mu")

# G_1 of Table 1 (paper.tex line 578), rows = who i sends to; vertex 7 receives from 2,3,5,6 (1-based)
rows = ["011000", "000011", "000111", "110011", "110101", "101110"]
n = 7
A = [[int(c) for c in r] + [0] for r in rows] + [[0] * 7]
for j in (2, 3, 5, 6):
    A[j - 1][6] = 1
sigma = list(range(6)); m = 6
assert all(A[6][j] == 0 for j in range(7))

# (beta,t) form: N_ab = 1 iff a != b and sigma_b -/-> sigma_a ; C = I + beta N ; K = C + t J
N = sp.Matrix(m, m, lambda i, j: 1 if (i != j and A[j][i] == 0) else 0)
C = sp.eye(m) + b * N
K = C + t * sp.ones(m, m)
adjC = C.adjugate()
a = [sp.expand(sum(adjC[i, j] for j in range(m))) for i in range(m)]
d = sp.expand(C.det()); s = sp.expand(sum(a))
D = sp.expand(K.det())
assert sp.expand(D - (d + t * s)) == 0
# o_7 = d - beta * sum_{b: sigma_b -/-> 7} a_b
o7 = sp.expand(d - b * sum(a[j] for j in range(m) if A[j][6] == 0))
claimed = -b**5 + 3*b**4 - 3*b**3 + 3*b**2 - 2*b + 1
nonsend = [j + 1 for j in range(m) if A[j][6] == 0]
out = {"nonsenders_to_7_1based": nonsend,
       "o7": str(o7), "o7_equals_claimed": sp.expand(o7 - claimed) == 0,
       "o7_equals_d_minus_beta_a1_plus_a4": sp.expand(o7 - (d - b * (a[0] + a[3]))) == 0}

lo, hi = sp.Rational(52, 25), sp.Rational(209, 100)
out["o7_at_lo"] = str(o7.subs(b, lo)); out["o7_at_hi"] = str(o7.subs(b, hi))
assert o7.subs(b, lo) == sp.Rational(417793, 9765625)
assert o7.subs(b, hi) == sp.Rational(-1006161749, 10**10)
out["o7_roots_in_closed_segment"] = sp.Poly(o7, b).count_roots(lo, hi)

# polynomials of the old list at t = 100
T = 100
pt = sp.Poly(sp.expand((mu * sp.eye(m) + K.subs(t, T)).det()), mu)
c = [sp.expand(x) for x in pt.all_coeffs()]          # c[0] = 1 (monic), c[k] coeff of mu^(6-k)
assert c[0] == 1
def coef(k):
    return c[k] if 0 <= k <= m else 0
def hurwitz(k):
    H = sp.Matrix(k, k, lambda i, j: coef(2 * (j + 1) - (i + 1)))
    return sp.expand(H.det())
Delta = [hurwitz(k) for k in range(1, m + 1)]
L = T * (b - 2) - 1
fam = {f"a{i+1}": a[i] for i in range(m)}
fam["D"] = sp.expand(D.subs(t, T))
for k in range(m):
    fam[f"Delta{k+1}"] = Delta[k]
fam["L"] = sp.expand(L)

# independent Sturm count over Fraction on [lo, hi]
def to_fr(p):
    return [F(int(sp.numer(x)), int(sp.denom(x))) for x in sp.Poly(p, b).all_coeffs()]
def ev(p, x):
    v = F(0)
    for co in p:
        v = v * x + co
    return v
def prem(f, g):
    f = f[:]
    while len(f) >= len(g) and any(f):
        q = f[0] / g[0]
        for i in range(len(g)):
            f[i] -= q * g[i]
        f.pop(0)
    while f and f[0] == 0:
        f.pop(0)
    return f
def deriv(p):
    k = len(p) - 1
    return [co * (k - i) for i, co in enumerate(p[:-1])]
def sturm_count(p, x0, x1):
    seq = [p, deriv(p)]
    while len(seq[-1]) > 1:
        r = prem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-x for x in r])
    def V(x):
        vals = [ev(q, x) for q in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for u, w in zip(vals, vals[1:]) if (u > 0) != (w > 0))
    # (x0, x1]; endpoints checked separately
    return V(x0) - V(x1)

Lf, Hf = F(52, 25), F(209, 100)
res = {}
allok = True
for name, p in fam.items():
    P = sp.Poly(p, b)
    sq = sp.Poly(sp.sqf_part(sp.Poly(p, b)).as_expr(), b) if sp.Poly(p, b).degree() > 0 else sp.Poly(p, b)
    fr = to_fr(sq.as_expr())
    vlo, vhi = ev(to_fr(p), Lf), ev(to_fr(p), Hf)
    n_sym = P.count_roots(lo, hi)
    n_st = sturm_count(fr, Lf, Hf) + (1 if ev(fr, Lf) == 0 else 0) if P.degree() > 0 else 0
    res[name] = {"degree": P.degree(), "sign_lo": (vlo > 0) - (vlo < 0), "sign_hi": (vhi > 0) - (vhi < 0),
                 "roots_closed_sympy": int(n_sym), "roots_closed_sturm": int(n_st)}
    allok &= vlo > 0 and vhi > 0 and n_sym == 0 and n_st == 0
out["listed_family_t100"] = res
out["listed_family_all_positive_no_roots"] = bool(allok)
# sturm sanity: o7 must have a root there
fo = to_fr(sp.sqf_part(sp.Poly(o7, b)).as_expr())
out["o7_roots_sturm"] = sturm_count(fo, Lf, Hf)

# direct check from CLAUDE.md conventions at (eps, delta) = (1/(1+t), (beta-1)/(1+t)), theta = 1
def direct(beta):
    eps = F(1, T + 1); delta = (beta - 1) / (T + 1)
    W = [[F(0) if i == j else (F(-1) + eps if A[j][i] else F(-1) - delta) for j in range(n)] for i in range(n)]
    M = sp.Matrix(m, m, lambda i, j: sp.Rational((1 if i == j else 0) - W[i][j].numerator, 1) if False else 0)
    M = sp.Matrix(m, m, lambda i, j: sp.Rational(((1 if i == j else 0) - W[i][j]).numerator, ((1 if i == j else 0) - W[i][j]).denominator))
    x = M.LUsolve(sp.ones(m, 1))
    pos = all(v > 0 for v in x)
    drive = 1 + sum(sp.Rational(W[6][j].numerator, W[6][j].denominator) * x[j] for j in range(m))
    J = -sp.eye(m) + sp.Matrix(m, m, lambda i, j: sp.Rational(W[i][j].numerator, W[i][j].denominator))
    cp = sp.Poly((mu * sp.eye(m) - J).det(), mu).all_coeffs()  # roots = eigenvalues of J
    cf = lambda k: cp[k] if 0 <= k <= m else 0
    Hs = [sp.Matrix(k, k, lambda i, j: cf(2 * (j + 1) - (i + 1))).det() for k in range(1, m + 1)]
    stab = "stable" if all(h > 0 for h in Hs) else ("undecided" if any(h == 0 for h in Hs) else "not_stable")
    # drive should equal theta * o7/D / ... check proportionality sign with D>0
    return {"x_positive": pos, "drive": str(drive), "drive_negative": drive < 0, "stability": stab,
            "SFP": bool(pos and drive < 0 and stab == "stable")}
out["direct_lo"] = direct(Lf); out["direct_hi"] = direct(Hf)
assert out["direct_lo"]["SFP"] is False and out["direct_hi"]["SFP"] is True
print(json.dumps(out, indent=1, default=str))
