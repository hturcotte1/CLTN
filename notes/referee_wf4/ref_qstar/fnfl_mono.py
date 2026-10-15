"""Interior critical points of Q(beta) = (beta-1) s(beta)/(s(beta)+1) along f1 = 0 in (2, beta_pos).
Own resultant + Sturm count. Also: does lc_t(f1) vanish at beta_pos? Exact."""
import pickle, time
import sympy as sp
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_qstar/'
b, t, x = sp.symbols('b t x')
d = pickle.load(open(OUT + 'fnfl_f.pkl', 'rb'))
f1 = sp.Poly(sp.sympify(d['f1']), b, t)
bpl, bph = sp.Rational(d['bpl']), sp.Rational(d['bph'])
m = sp.Poly(b**5 - 3*b**4 + 3*b**3 - 3*b**2 + 2*b - 1, b)
lc = sp.Poly(f1.as_expr().coeff(t, f1.degree(t)), b)
print('lc_t(f1) =', sp.factor(lc.as_expr()), ' roots in [2, bph]:', lc.count_roots(2, bph))
disc = sp.Poly(sp.discriminant(f1.as_expr(), t), b)
print('disc_t(f1) roots in [2,bph] (closed):', sp.Poly(disc, b).count_roots(2, bph))
fe = f1.as_expr()
h = sp.expand(t * (t + 1) * sp.diff(fe, t) - (b - 1) * sp.diff(fe, b))
t0 = time.time()
R = sp.Poly(sp.resultant(fe, h, t), b)
print('Res_t(f1,h) degree', R.degree(), 'time', time.time() - t0)
print('roots of Res in closed [2, beta_pos_hi]:', R.count_roots(2, bph), ' in [bpl,bph]:', R.count_roots(bpl, bph))
# also check the legality line: f1 meeting t = 1/(b-2) inside (2, beta_pos)?  (would split the cell)
L = sp.Poly(sp.numer(sp.together(fe.subs(t, 1 / (b - 2)))), b)
print('f1 on legality line: roots in (2,bph]:', L.count_roots(sp.Rational(2) + sp.Rational(1, 10**30), bph))
