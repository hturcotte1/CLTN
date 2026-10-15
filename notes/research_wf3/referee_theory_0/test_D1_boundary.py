import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/n7')
import sympy as sp, random, ctln, fastregion as fr
from curve_test import from_d6
sys.argv = sys.argv[:1]; exec(open("test_D1_D2.py").read().split("stats = {")[0])
y = sp.Symbol('y')
def verdicts(A, beta):
    m = len(A); As = sp.Matrix(A); D = sp.ones(m, m) - sp.eye(m) - As
    bd = fr.beta_data(A, list(range(m)))
    v_repo = ctln.hurwitz_stable([fr.peval(p, beta) for p in fr.stability_coeffs(bd, fr.INF)])[0]
    v_D = shifted_hurwitz(psi(D), sp.Rational(beta.denominator, beta.numerator))
    return v_repo, v_D
# FNFL: alpha = -1/2 exactly, so at beta = 2 boundary
A = from_d6('&FNFL~[udo?')
for b in [F(2), F(199,100), F(201,100), F(3)]:
    print('FNFL', b, verdicts(A, b))
# Geneson graph and K_n limit check
g = open('/home/user/CLTN/n7/geneson.d6').read().split()[0]
A = from_d6(g)
for b in [F(2), F(21,10), F(11,5), F(5,2)]:
    print('Geneson', b, verdicts(A, b))
# random near-complete graphs: count True verdict agreement
random.seed(7); cnt = {True:0, False:0, None:0}
for _ in range(200):
    n = random.randint(3, 7)
    A = [[0 if i==j else (1 if random.random() < 0.9 else 0) for j in range(n)] for i in range(n)]
    for b in [F(2), F(5,2), F(3), F(9,4)]:
        v = verdicts(A, b); assert v[0] == v[1], (A, b, v); cnt[v[0]] += 1
print(cnt)
