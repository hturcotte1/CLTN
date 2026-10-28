"""Exact count of roots of an integer/rational polynomial by real part sign (Re < 0, Re = 0, Re > 0), with multiplicity,
via sympy exact complex root isolation (rational rectangles) plus an exact axis test gcd(Re p(iy), Im p(iy))."""
from fractions import Fraction as F
import sympy as sp
y = sp.Symbol('y'); w = sp.Symbol('w')

def root_signs(coeffs_hi_to_lo):
    p = sp.Poly([sp.Rational(c.numerator, c.denominator) if isinstance(c, F) else sp.Rational(c) for c in coeffs_hi_to_lo], w)
    neg = zero = pos = 0
    for fac, mult in p.sqf_list()[1]:
        if fac.degree() == 0: continue
        # axis roots of this squarefree factor: real y with fac(i y) = 0
        cl = fac.all_coeffs()[::-1]  # lo to hi: c_k w^k
        rec = {}; imc = {}
        for k, c in enumerate(cl):
            # (i y)^k = i^k y^k ; i^k real for k even: (-1)^(k/2); imaginary for k odd: i (-1)^((k-1)/2)
            if k % 2 == 0: rec[k] = c * (-1) ** (k // 2)
            else: imc[k] = c * (-1) ** ((k - 1) // 2)
        re = sp.Poly(sum(c * y**k for k, c in rec.items()), y) if rec else sp.Poly(0, y)
        im = sp.Poly(sum(c * y**k for k, c in imc.items()), y) if imc else sp.Poly(0, y)
        g = sp.gcd(re, im)
        n_axis = 0
        if g.degree() > 0:
            n_axis = sp.Poly(g, y).count_roots()  # real roots of the gcd: each is an axis root i*y0 of fac (fac squarefree -> distinct)
        # isolate all roots of fac in rectangles / intervals; refine until no rectangle straddles Re = 0 except the axis ones
        eps = sp.Rational(1, 2**20)
        for attempt in range(6):
            real_iv, cplx_iv = fac.intervals(all=True, eps=eps)
            n_neg = n_pos = n_strad = 0
            for (a, b), m in real_iv:
                if b < 0: n_neg += 1
                elif a > 0: n_pos += 1
                else:
                    if fac.eval(0) == 0: pass  # the axis root 0 (counted in n_axis)
                    else: n_strad += 1
            for (c1, c2), m in cplx_iv:
                ax, ay, bx, by = sp.re(c1), sp.im(c1), sp.re(c2), sp.im(c2)
                if bx < 0: n_neg += 1
                elif ax > 0: n_pos += 1
                else: n_strad += 1
            # straddling boxes must be exactly the axis roots
            if n_strad + (1 if fac.eval(0) == 0 else 0) == n_axis:
                break
            eps = eps / 2**20
        else:
            raise RuntimeError('could not separate roots from the axis')
        neg += mult * n_neg; pos += mult * n_pos; zero += mult * n_axis
    assert neg + zero + pos == p.degree(), (neg, zero, pos, p.degree())
    return neg, zero, pos

if __name__ == '__main__':
    # calibration
    print(root_signs([1, 0, 1]))            # w^2 + 1 -> (0, 2, 0)
    print(root_signs([1, 2, 1]))            # (w+1)^2 -> (2, 0, 0)
    print(root_signs([1, 0, -1]))           # w^2 - 1 -> (1, 0, 1)
    print(root_signs([1, 0, 1, 1]))         # w^3 + w + 1: one real negative root, two complex with Re > 0
    print(root_signs([1, 1, 1, 1]))         # (w+1)(w^2+1) -> (1,2,0)
    print(root_signs([F(1), F(3, 2), F(1, 2)]))  # (w+1)(w+1/2)
