"""Exact isolation of the roots of chi_perp (charpoly of A compressed to 1-perp):
7 chi_perp(l) = chi_A(l) - chi_{A+J}(l). Prints exact isolating rectangles for every root and
decides alpha(A|perp) vs a rational threshold rigorously (rectangle strictly on one side)."""
import sys, json
import sympy as sp
from sympy.polys.rootoftools import ComplexRootOf
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf/lemmaM_refute')
from scan import from_d6

def chi_perp(A):
    n = len(A); l = sp.symbols('l')
    M = sp.Matrix(A); J = sp.ones(n, n)
    p = sp.expand((l * sp.eye(n) - M).det() - (l * sp.eye(n) - M - J).det())
    return sp.Poly(p, l)

def alpha_vs(A, thr, prec=40):
    P = chi_perp(A)
    thr = sp.Rational(thr)
    out = []
    best = None
    for f, mult in sp.factor_list(P)[1]:
        for k in range(f.degree()):
            r = ComplexRootOf(f, k)
            if r.is_Rational:
                side = '>' if r > thr else ('<' if r < thr else '=')
                out.append(dict(factor=str(f.as_expr()), k=k, real=True, re_lo=str(r), re_hi=str(r), side=side, approx=str(r)))
                continue
            # refine until rectangle is strictly on one side of thr or small
            for _ in range(60):
                iv = r._get_interval()
                if r.is_real:
                    lo, hi = iv.a, iv.b
                else:
                    lo, hi = iv.ax, iv.bx
                if lo > thr or hi < thr: break
                r.eval_rational(dx=sp.Rational(1, 10 ** (_ // 3 + 2)), dy=sp.Rational(1, 10 ** (_ // 3 + 2)))
            iv = r._get_interval()
            lo, hi = (iv.a, iv.b) if r.is_real else (iv.ax, iv.bx)
            side = '>' if lo > thr else ('<' if hi < thr else '?')
            out.append(dict(factor=str(f.as_expr()), k=k, real=bool(r.is_real), re_lo=str(lo), re_hi=str(hi), side=side,
                            approx=str(sp.N(r, 15))))
    return P, out

if __name__ == '__main__':
    A = from_d6(sys.argv[1]); thr = sys.argv[2]
    P, out = alpha_vs(A, thr)
    print('7*chi_perp =', P.as_expr())
    for o in out: print(json.dumps(o))
    sides = [o['side'] for o in out]
    print('VERDICT alpha', '>' if '>' in sides else ('<' if all(s == '<' for s in sides) else '?'), thr)
