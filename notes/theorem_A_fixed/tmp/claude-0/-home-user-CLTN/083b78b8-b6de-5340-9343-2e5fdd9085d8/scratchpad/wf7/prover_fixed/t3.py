import fixlib as fl, sympy as sp, time
A = fl.from_d6('&F^U|m^QTq?'); S=(1,1,0,1,1,1,1)
g, j0 = fl.reduced_psi(A,[S],[1]); print('j0', j0, sp.factor(g))
for beta1 in [sp.Rational(21,10), sp.Rational(43,20), sp.Rational(216,100)]:
    r1 = 1/beta1; a1 = fl.legal_a_bound(r1)
    co, mins = fl.region_minors(g, r1)
    print(beta1, 'a1', a1, [fl.sign_pattern_certificate(mn, a1, 0)[:2] for mn in mins])
