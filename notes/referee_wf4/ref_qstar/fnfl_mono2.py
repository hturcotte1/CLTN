"""Interior critical points of Q(beta) = (beta-1)s/(s+1) along f1 = 0 inside (2, beta_pos): own
resultant, factor, exact isolation.  Also f1 against the legality line."""
import pickle, time
import sympy as sp
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_qstar/'
b, t = sp.symbols('b t')
d = pickle.load(open(OUT + 'fnfl_f.pkl', 'rb'))
f1 = sp.Poly(sp.sympify(d['f1']), b, t)
bpl, bph = sp.Rational(d['bpl']), sp.Rational(d['bph'])
fe = f1.as_expr()
h = sp.expand(t * (t + 1) * sp.diff(fe, t) - (b - 1) * sp.diff(fe, b))
R = sp.Poly(sp.resultant(fe, h, t), b)
print('Res degree', R.degree(), flush=True)
t0 = time.time()
c, fl = sp.factor_list(R)
print('factor time', time.time() - t0, [(sp.Poly(f, b).degree(), mm) for f, mm in fl], flush=True)
tot = 0
for f, mm in fl:
    P = sp.Poly(f, b)
    ivs = P.intervals(inf=2, sup=bph)
    tot += len(ivs)
    print('deg', P.degree(), 'roots in [2,bph]:', [(str(i[0][0]), str(i[0][1])) for i in ivs], flush=True)
print('TOTAL roots of Res_t(f1,h) in [2, bph]:', tot)
L = sp.Poly(sp.numer(sp.together(fe.subs(t, 1 / (b - 2)))), b)
cL, fL = sp.factor_list(L)
print('f1 on legality line, roots in [2,bph]:', [(sp.Poly(f, b).degree(), [(str(i[0][0]), str(i[0][1])) for i in sp.Poly(f, b).intervals(inf=2, sup=bph)]) for f, mm in fL])
