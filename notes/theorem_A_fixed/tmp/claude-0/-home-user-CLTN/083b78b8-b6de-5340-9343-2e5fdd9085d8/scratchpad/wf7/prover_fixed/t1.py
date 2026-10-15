import fixlib as fl, sympy as sp, time
from fractions import Fraction as F
B = fl.from_d6('&F^S{^ReL|?'); SQ=(1,1,1,1,0,1,1)
t0=time.time()
Psi = fl.Psi_expr(B,[SQ],[1])
print(sp.factor(Psi.coeff(fl.k)), time.time()-t0)
def k0_at(Psi, beta, t, kmin=0):
    beta=sp.Rational(beta); t=sp.Rational(t)
    co = fl.shifted_coeffs(Psi.subs(fl.a_s, t/beta), 1/beta)
    mins = fl.hurwitz_minors(co)
    k0 = kmin
    for mn in mins:
        kk = fl.smallest_k0(mn, kmin)
        if kk is None: return None
        k0 = max(k0, kk)
    return k0
for beta in [F(21,10), F(11,5), F(23,10), F(12,5), F(5,2), F(13,5), F(27,10), F(14,5)]:
    row=[]
    for t in [F(1,1)/(beta-2)*F(11,10), 5, 10, 20, 50, 100, 1000, 10**6]:
        if t*(beta-2) <= 1: row.append('-'); continue
        row.append(k0_at(Psi, beta, t))
    print(beta, row, time.time()-t0)
