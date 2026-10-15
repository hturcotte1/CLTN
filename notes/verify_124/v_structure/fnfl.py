import json, sys
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf6/v_structure')
from exact import *
import sympy as sp
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf6/v_structure/'
A = d6_decode('&FNFL~[udo?'); n = 7
AJ = [[A[i][j] + 1 for j in range(n)] for i in range(n)]
mchi = poly_sub(charpoly(A), charpoly(AJ))
chi = [c / n for c in mchi]
res = {'d6': '&FNFL~[udo?', 'rows': [''.join(map(str, r)) for r in A], 'm_chi_perp': [str(c) for c in mchi],
       'chi_perp': [str(c) for c in chi]}
# own polynomial long division by l^2 + l + 1
def pdiv(p, d):
    p = list(p); qd = []
    while len(p) >= len(d):
        c = p[0] / d[0]; qd.append(c)
        for k in range(len(d)): p[k] -= c * d[k]
        p.pop(0)
    return qd, p
quo, rem = pdiv(chi, [F(1), F(1), F(1)])
res['quotient_by_l2+l+1'] = [str(c) for c in quo]; res['remainder'] = [str(c) for c in rem]
assert all(r == 0 for r in rem), rem
# l^2+l+1 shifted by -1/2: roots exactly -1/2 +- i sqrt(3)/2
sh = poly_shift([1, 1, 1], F(-1, 2)); res['l2+l+1_at_x-1/2'] = [str(c) for c in sh]
assert sh == [1, 0, F(3, 4)]
# quotient: every root real part < -1/2  <=> quo(x - 1/2) Hurwitz
st, minors = hurwitz(poly_shift(quo, F(-1, 2)))
res['quotient_shifted_hurwitz'] = st; res['quotient_shifted_minors'] = [str(m) for m in minors]
# l^2+l+1 does not divide the quotient again
q2, r2 = pdiv(quo, [F(1), F(1), F(1)]); res['l2+l+1_divides_quotient'] = all(r == 0 for r in r2)
# sympy factorization (assistant, cross-check only)
l = sp.symbols('l')
P = sp.Poly([sp.Rational(c.numerator, c.denominator) for c in mchi], l)
res['sympy_factor_list_of_m_chi_perp'] = str(sp.factor_list(P.as_expr()))
# also full chi_perp spectral abscissa: alpha = -1/2 exactly; full shifted polynomial at shift 1/2 must have zero minor
st2, m2 = hurwitz(poly_shift(chi, F(-1, 2)))
res['chi_perp_shift_1/2_hurwitz'] = st2
# limit stability needs alpha < -gamma with gamma = 1 - 1/beta > 1/2 for legal beta > 2: impossible since alpha = -1/2
for b in ['2', '201/100', '53/26', '3', '100']:
    beta = F(b); gamma = 1 - 1 / beta
    s, _ = hurwitz(poly_shift(chi, -gamma))
    res.setdefault('limit_shifted_hurwitz_at_beta', {})[b] = s
# finite witnesses from the FNFL certificate
def W_of(A, eps, delta):
    return [[F(0) if i == j else ((-1 + eps) if A[j][i] == 1 else (-1 - delta)) for j in range(n)] for i in range(n)]
wit = []
for eps, delta in [(F(1, 49), F(33, 1568)), (F(1, 36), F(17, 588))]:
    q = delta * (1 - eps) / eps
    W = W_of(A, eps, delta)
    IW = [[F(int(i == j)) - W[i][j] for j in range(n)] for i in range(n)]
    x = solve(IW, [1] * n)
    st3, _ = hurwitz(charpoly([[W[i][j] - int(i == j) for j in range(n)] for i in range(n)]))
    t = (1 - eps) / eps; beta = 1 + delta / eps
    wit.append({'eps': str(eps), 'delta': str(delta), 'q': str(q), 'beta': str(beta), 't': str(t),
                'legal': 0 < eps < 1 and delta > 0 and q > 1, 'positive': x is not None and all(v > 0 for v in x),
                'stable': st3, 'x': [str(v) for v in x]})
res['finite_witnesses'] = wit
print(json.dumps(res, indent=1))
json.dump(res, open(OUT + 'fnfl_out.json', 'w'), indent=1)
