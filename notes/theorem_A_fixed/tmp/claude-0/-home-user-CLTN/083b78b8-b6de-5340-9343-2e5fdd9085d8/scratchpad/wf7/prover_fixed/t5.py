import fixlib as fl, sympy as sp
A = fl.from_d6('&F^U|m^QTq?'); S=(1,1,0,1,1,1,1)
g, j0 = fl.reduced_psi(A,[S],[1])
print('g =', sp.collect(g, fl.z))
gk = g.subs(fl.a_s, sp.Rational(200,21))
co, mins = fl.region_minors(gk, sp.Rational(10,21))
print('h coeffs (x):', [sp.factor(c) for c in co])
for j, mn in enumerate(mins):
    P = sp.Poly(mn, fl.k); den = sp.lcm([sp.fraction(c)[1] for c in P.all_coeffs()])
    print('Delta_%d * %s =' % (j+1, den), sp.Poly(sp.expand(mn*den), fl.k).as_expr())
# region version at beta1=43/20 with a = 400/129 + u
co2, mins2 = fl.region_minors(g, sp.Rational(20,43))
for j, mn in enumerate(mins2):
    q = sp.Poly(sp.expand(mn.subs(fl.a_s, sp.Rational(400,129)+fl.u_s)), fl.u_s, fl.k)
    print('region minor', j+1, 'total degree', q.total_degree(), 'n terms', len(q.terms()), 'min coeff', min(q.coeffs()))
