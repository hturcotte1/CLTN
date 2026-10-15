"""Referee ref_Q: scale law (5.3), many-class trace ceiling (5.5), family ceilings (display values bracketed exactly),
and the limit polynomial == det(zI - D + bbar 1^T) identity. Exact."""
import random, json
import sympy as sp
from fractions import Fraction as F
import fam, exact, sym
z, r = sp.symbols('z r')
out = {}
# --- scale law on random real(=rational/integer) weights
random.seed(1)
ok_phi = ok_u = 0; trials = 0
for _ in range(25):
    m = random.randint(3, 5)
    D = sp.Matrix(m, m, lambda i, j: 0 if i == j else random.randint(0, 1))
    mu = [random.randint(1, 4) for _ in range(m)]; N = random.randint(2, 3)
    def Phi(muv):
        M = sp.diag(*muv); return sp.expand((sp.Matrix([muv]) * (z * sp.eye(m) - D * M).adjugate() * sp.ones(m, 1))[0])
    def u(muv, rr):
        M = sp.diag(*muv); A = D.T * M + rr * sp.eye(m)
        if A.det() == 0: return None
        return A.LUsolve(sp.ones(m, 1))
    P1 = Phi([N * v for v in mu]); P0 = Phi(mu)
    ok_phi += sp.expand(P1 - N ** m * P0.subs(z, z / N)) == 0
    rr = sp.Rational(random.randint(1, 9), random.randint(10, 30))
    u1 = u([N * v for v in mu], rr); u0 = u(mu, rr / N)
    if u1 is not None and u0 is not None:
        ok_u += all(sp.simplify(a - c / N) == 0 for a, c in zip(u1, u0)) and \
                sp.simplify(sum(N * mu[i] * u1[i] for i in range(m)) - sum(mu[i] * u0[i] for i in range(m))) == 0
    trials += 1
out['scale_law'] = dict(trials=trials, Phi_identity_ok=ok_phi, u_identity_ok=ok_u)
# --- limit polynomial == p' * det(zI - D + bbar 1^T) and trace bound
for nm in ('Q', 'Q2', 'Y'):
    f = sym.FAMS[nm]; B = f['B']; m = len(B)
    Phi, zpsi, Lk, chi = sym.family_Phi(B, f['classes'])
    pk = sum(c for _, c in f['classes'])
    D = sp.Matrix(sym.Dmat(B))
    bbar = sum((sp.Rational(c, pk) * sp.Matrix([1 - s for s in S]) for S, c in f['classes']), sp.zeros(m, 1))
    Lmat = D - bbar * sp.ones(1, m)
    same = sp.expand(Lk - pk * (z * sp.eye(m) - Lmat).det()) == 0
    out[nm] = dict(limit_eq_det=same, trace_L=str(Lmat.trace()), trace_bound_minRe_le=str(Lmat.trace() / m))
# --- ceilings: bracket -1/minRe of the limit polynomial with exact Routh tests (display values from notes)
def stable_at(poly, beta):
    co = [F(int(c.p), int(c.q)) for c in reversed(sp.Poly(poly, z).all_coeffs())]
    return exact.routh(exact.poly_shift_neg(co, 1 / beta))
brk = {}
brk['Q quartic, beta 2.8441 / 2.8442'] = (stable_at(z**4 + z**3 + 3*z**2 + z + 1, F(28441, 10000)), stable_at(z**4 + z**3 + 3*z**2 + z + 1, F(28442, 10000)))
brk['Q2 (2z^2+z+2), beta 3999/1000 / 4'] = (stable_at(2*z**2 + z + 2, F(3999, 1000)), stable_at(2*z**2 + z + 2, F(4)))
Y6 = 5*z**6 + 5*z**5 + 12*z**4 + 7*z**3 + 7*z**2 + 2*z + 1
lo, hi = F(4), F(6)
for _ in range(40):
    mid = (lo + hi) / 2
    if stable_at(Y6, mid): lo = mid
    else: hi = mid
brk['Y sextic ceiling bracket'] = (str(float(lo)), str(float(hi)), 'display only')
brk['Y sextic, beta 24/5 stable'] = stable_at(Y6, F(24, 5))
out['ceilings'] = brk
print(json.dumps(out, indent=1, default=str))
json.dump(out, open('misc_results.json', 'w'), indent=1, default=str)
