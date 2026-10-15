"""Referee: is beta_c(k) strictly decreasing for k >= 8, and -> 2?  r = 1/beta; r_c(k) = largest root in (0,1/2) of D3(k, .).
Exact resultants/root counting; floats only for display."""
import sympy as sp
z, k, x, r = sp.symbols('z k x r')
txt = open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/referee_theory_2/minors.txt').read().splitlines()
D3 = sp.sympify(txt[0].split('=', 1)[1]); a4 = sp.sympify(txt[1].split('=', 1)[1])
c3 = sp.Poly(D3, k).all_coeffs()[0]
print('leading coeff in k (c3(r)) =', sp.factor(c3))
Dk = sp.diff(D3, k); Dr = sp.diff(D3, r)
res_k = sp.Poly(sp.resultant(D3, Dk, r), k)
res_r = sp.Poly(sp.resultant(D3, Dr, r), k)
print('Res_r(D3, dD3/dk): degree', res_k.degree(), ' real roots >= 8:', res_k.count_roots(8, None) if False else None)
def roots_ge(P, lo):
    return [rt for rt in P.real_roots() if rt >= lo]
rk = roots_ge(res_k, 4); rr = roots_ge(res_r, 4)
print('real roots k >= 4 of Res(D3, D3_k):', [sp.N(t, 8) for t in rk])
print('real roots k >= 4 of Res(D3, D3_r) (discriminant):', [sp.N(t, 8) for t in rr])
# per-k values of r_c (display) for a range
import math
for kk in [4, 5, 6, 8, 9, 20, 100, 1000, 10**5]:
    P = sp.Poly(D3.subs(k, kk), r)
    rts = [t for t in P.real_roots() if 0 < t < sp.Rational(1, 2)]
    rc = max(rts, key=lambda t: sp.N(t, 30))
    print(kk, 'beta_c ~', sp.N(1/rc, 10), ' k*(beta_c-2) ~', sp.N(kk*(1/rc - 2), 6))
# a4 never binds: a4 = k r^2 (r^2 - r + 1) + p0(r) with p0 > 0 on [0, 1/2]
p0 = sp.Poly(6*r**4 - 11*r**3 + 16*r**2 - 10*r + 3, r)
assert sp.expand(a4 - (k*r**2*(r**2 - r + 1) + p0.as_expr())) == 0
assert p0.count_roots(0, sp.Rational(1, 2)) == 0 and p0.eval(0) > 0
assert sp.Poly(64*(k + 6)**3, k) == sp.Poly(sp.Poly(D3, r).all_coeffs()[0], k)
assert len(rk) == 1 and rk[0] < 5 and len(rr) == 0
print('a4 > 0 for all k >= 0, r in [0,1/2]; lead_r(D3) = 64(k+6)^3 > 0; unique critical k in [4, inf) is < 5, no double roots: beta_c strictly monotone on k >= 5 (decreasing, as beta_c(5) > beta_c(6)), -> 2 since c3(r) < 0 on (0, 1/2)')
