"""Referee ref_Q: independent exact re-check of the four stored witnesses (no ctln, no theory scripts):
 (i) W built from CLAUDE.md section 2 verbatim; I - W == eps*(I + beta N + t J) entrywise;
 (ii) legality, x = (I-W)^{-1} 1 > 0 (Fraction Gauss);
 (iii) finite-(eps,delta) stability: charpoly of the integer matrix s*K (Bareiss + interpolation), Routh array on p(-y);
 (iv) exact limit test from the explicit adjacency: chi of A|perp in the basis e_i - e_n, Routh on chi(-gamma + y);
      plus psi_{D'} == z^(p-1) Phi(z,k) and n*chi_{D'|perp} == psi_{D'}."""
import sys, json, time
from fractions import Fraction as F
import sympy as sp
import fam, exact, sym

def W_of(A, eps, delta):
    n = len(A)
    return [[F(0) if i == j else ((-1 + eps) if A[j][i] == 1 else (-1 - delta)) for j in range(n)] for i in range(n)]

def solve(M, rhs):
    n = len(M); a = [list(map(F, r)) + [F(v)] for r, v in zip(M, rhs)]
    for c in range(n):
        p = next(r for r in range(c, n) if a[r][c] != 0); a[c], a[p] = a[p], a[c]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c] / a[c][c]; a[r] = [u - f * v for u, v in zip(a[r], a[c])]
    return [a[i][n] / a[i][i] for i in range(n)]

def perp_matrix(A):
    """X with V X = P A V, V = [e_i - e_n], P = I - J/n; returns n*X as an integer matrix."""
    n = len(A)
    AV = [[A[i][j] - A[i][n - 1] for j in range(n - 1)] for i in range(n)]
    col = [sum(AV[i][j] for i in range(n)) for j in range(n - 1)]
    # (P A V)[i][j] = AV[i][j] - col[j]/n ; times n:
    return [[n * AV[i][j] - col[j] for j in range(n - 1)] for i in range(n - 1)]

def run(name, kk, A, beta, Phi, p):
    n = len(A); t = F(10000); eps = 1 / (t + 1); delta = (beta - 1) * eps
    rec = dict(family=name, k=kk, n=n, beta=str(beta), eps=str(eps), delta=str(delta))
    q = delta * (1 - eps) / eps; rec['q'] = str(q); rec['legal'] = q > 1 and t * (beta - 2) > 1
    W = W_of(A, eps, delta)
    N = [[1 if (i != j and not A[j][i]) else 0 for j in range(n)] for i in range(n)]
    K = [[(1 if i == j else 0) + beta * N[i][j] + t for j in range(n)] for i in range(n)]
    rec['IminusW_eq_epsK'] = all(((1 if i == j else 0) - W[i][j]) == eps * K[i][j] for i in range(n) for j in range(n))
    x = solve([[(1 if i == j else 0) - W[i][j] for j in range(n)] for i in range(n)], [1] * n)
    rec['x_positive'] = all(v > 0 for v in x); rec['x_min'] = str(min(x))
    s = beta.denominator
    sK = [[int(s * v) for v in r] for r in K]
    cp = exact.charpoly_int(sK)            # det(zI - sK), low->high
    # eigenvalues of -I+W are -eps*eig(K); stable iff sK positive stable iff roots of cp(-y) in Re y < 0
    pm = [c * (-1) ** i for i, c in enumerate(cp)][::-1]   # cp(-y) high->low (up to sign)
    rec['finite_stable_routh'] = exact.routh(pm)
    # limit test, directly from A
    nX = perp_matrix(A); cpp = exact.charpoly_int(nX)      # det(w I - nX), w = n*lambda
    gamma = 1 - 1 / beta
    lam = sp.symbols('lam'); y = sp.symbols('y')
    chi = sum(sp.Integer(c) * (n * lam) ** i for i, c in enumerate(cpp))          # = n^(n-1) chi_{A|perp}(lam)
    shifted = sp.Poly(sp.expand(chi.subs(lam, -sp.Rational(gamma.numerator, gamma.denominator) + y)), y)
    rec['limit_alpha_lt_minus_gamma_routh'] = exact.routh([F(int(c.p), int(c.q)) for c in shifted.all_coeffs()])
    # D-picture consistency: chi_{A|perp}(lam) = (-1)^(n-1) chi_{D|perp}(-1-lam); n chi_{D|perp} = psi_{D'} = z^(p-1) Phi
    z = sym.z
    psi = sp.Poly(list(reversed(exact.psi_int(sym.Dmat(A)))), z).as_expr()
    rec['psi_eq_zp1_Phi'] = sp.expand(psi - z ** (p - 1) * Phi.subs(sym.k, kk)) == 0
    chiA_from_psi = sp.expand((-1) ** (n - 1) * psi.subs(z, -1 - lam) / n)
    chiA_direct = sp.expand(chi / sp.Integer(n) ** (n - 1))
    rec['chiAperp_eq_from_psi'] = sp.expand(chiA_from_psi - chiA_direct) == 0
    # Routh on the reduced family polynomial at the same beta
    r = 1 / beta
    ph = sp.Poly(Phi.subs(sym.k, kk), z)
    co = [F(int(c)) for c in reversed(ph.all_coeffs())]
    rec['Phi_routh_at_beta'] = exact.routh(exact.poly_shift_neg(co, r))
    return rec

if __name__ == '__main__':
    out = []
    Phis = {nm: sym.family_Phi(sym.FAMS[nm]['B'], sym.FAMS[nm]['classes'])[0] for nm in ('Q', 'Q2', 'Y')}
    cases = [('Q2', 3, F(3), 4), ('Q', 25, F(5, 2), 1), ('Q2', 10, F(7, 2), 4), ('Y', 13, F(3), 5)]
    for nm, kk, beta, ppk in cases:
        t0 = time.time()
        rec = run(nm, kk, sym.FAMS[nm]['build'](kk), beta, Phis[nm], ppk * kk)
        rec['seconds'] = round(time.time() - t0)
        print(json.dumps(rec), flush=True); out.append(rec)
    json.dump(out, open('inst_results.json', 'w'), indent=1)
