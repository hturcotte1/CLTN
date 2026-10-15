import fixlib as fl, sympy as sp
from fractions import Fraction as F
from theorem_B import FAMILIES
for nm, pts in (('Q2', [(F(13,4),1), (F(7,2),1), (F(16,5),1)]), ('Y', [(F(13,4),1), (F(7,2),1), (F(4),1)])):
    base, classes, rhos = FAMILIES[nm]; A = fl.from_d6(base)
    for beta, t in pts:
        k0, mins = fl.k_threshold(A, classes, rhos, beta, t, kmin=1)
        eps, delta = fl.eps_delta(beta, t)
        print(nm, beta, t, 'eps', eps, 'delta', delta, 'q', fl.ctln.q_of(eps,delta), 'k0', k0, 'n0', len(A) + sum(rhos)*k0 if k0 else None, flush=True)
