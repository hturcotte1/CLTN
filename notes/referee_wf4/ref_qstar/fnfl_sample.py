"""Exact sampling of the FNFL stable region around the claimed supremum with ctln.analyze_support.
Floats only in printed labels."""
import sys, pickle
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
from fractions import Fraction as F
import ctln
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf4/ref_qstar/'
b, t, x = sp.symbols('b t x')
d = pickle.load(open(OUT + 'fnfl_f.pkl', 'rb'))
f1 = sp.Poly(sp.sympify(d['f1']), b, t)
bpl, bph = F(d['bpl']), F(d['bph'])
rows = ['0011110', '0011100', '1101111', '1110111', '0011011', '0100101', '1100000']
A = [[int(c) for c in r] for r in rows]
QA, QB = F(685003, 657781), F(2170535, 2084278)   # claimed sup in (QA, QB)
assert QA < QB
def R(fr): return sp.Rational(fr.numerator, fr.denominator)
def Fr(r): return F(int(r.p), int(r.q))
def verdict(beta, tt):
    assert tt * (beta - 2) > 1
    eps = 1 / (tt + 1); delta = (beta - 1) * eps
    r = ctln.analyze_support(A, list(range(7)), eps, delta)
    return (beta - 1) * tt / (tt + 1), (r['positive'] is True and r['off_support_ok'] is True and r['stable'] is True), r['stable'], r['positive']
def section(beta):
    """roots of f1(beta, .) above t_min, as rational isolating intervals"""
    u = sp.Poly(f1.as_expr().subs(b, R(beta)), t)
    tmin = 1 / (beta - 2)
    return [(Fr(i[0][0]), Fr(i[0][1])) for i in u.intervals(eps=sp.Rational(1, 10**30), inf=R(tmin))]
print('--- slice structure check: just below / just above the f1 section')
for beta in [F(201, 100), F(53, 26), F(206, 100), F(208, 100), bpl - F(1, 10**6)]:
    secs = section(beta)
    s_lo, s_hi = secs[0]
    tmin = 1 / (beta - 2)
    below = verdict(beta, s_lo * (1 - F(1, 10**9)))
    above = verdict(beta, s_hi * (1 + F(1, 10**9)))
    mid = verdict(beta, (tmin + s_lo) / 2)
    print(f'beta~{float(beta):.6f} tmin~{float(tmin):.3f} sections~{[float(a) for a, _ in secs]} mid:{mid[1]} justbelow:{below[1]} justabove:{above[1]} (stable={above[2]}, pos={above[3]})')
print('--- BELOW: approach the corner (beta_pos, s(beta_pos))')
for h in [F(1, 10**4), F(1, 10**6), F(1, 10**9), F(1, 10**12)]:
    beta = bpl - h
    s_lo, s_hi = section(beta)[0]
    tt = s_lo * (1 - F(1, 10**15))
    q, ok, st, pos = verdict(beta, tt)
    print(f'h=1e-{len(str(h.denominator))-1}: stable positive fp: {ok}; q < QA ? {q < QA}; QA - q ~ {float(QA - q):.3e}')
print('--- ABOVE: levels ell > QB; points with q >= ell need t >= tau_ell(beta) = ell/(beta-1-ell)')
bad = 0; tested = 0; und = 0
for ell in [QB, QB + F(1, 10**9), QB + F(1, 10**6), QB + F(1, 10**3), F(106, 100)]:
    for k in range(1, 400):
        beta = 2 + (bpl - 2) * F(k, 400)
        if beta - 1 - ell <= 0: continue
        tau = ell / (beta - 1 - ell)
        for tt in [tau, tau * F(1000001, 1000000), tau * F(11, 10), tau * 2, tau * 10]:
            if tt * (beta - 2) <= 1: continue
            q, ok, st, pos = verdict(beta, tt)
            assert q >= ell
            tested += 1; bad += ok; und += (st is None)
    for beta in [bph + F(1, 10**9), bph + F(1, 10**4), F(21, 10)]:
        for tt in [F(13), F(25), F(30), F(100), F(10**6)]:
            if tt * (beta - 2) <= 1 or (beta - 1) * tt / (tt + 1) < ell: continue
            q, ok, st, pos = verdict(beta, tt)
            tested += 1; bad += ok; und += (st is None)
print(f'ABOVE: tested {tested} legal points with q >= ell > QB: stable-positive {bad}, undecided {und}')
