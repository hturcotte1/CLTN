"""Exact first-order (Lemma AS, notes/lemma_M_refuted.md section 4) data at the cube-root boundary.
A-picture: chi(z) = det(zI + A), psi(z) = 1^T adj(zI + A) 1 = det(zI + A + J) - det(zI + A).
nu_A(c) = min Re spec(cJ - A);  roots of psi are the limits (c -> inf) of the moving eigenvalues.
If z0 is a simple root of psi, not a root of chi, then the eigenvalue branch is z0 + 1/(c s'(z0)) + O(c^-2),
s = psi/chi, and 1/s'(z0) = chi(z0)/psi'(z0).  R := Re(chi(z0)/psi'(z0)).
Legal finite-t stability near (beta, t) -> (2, inf) needs nu_A -> 1/2 from above faster than nu_P2(c) = 1/2 + 1/(8c) + ...,
i.e. R > 1/8 at every root z0 with Re z0 = 1/2 (Theorem AS)."""
import sympy as sp
import tlib
z = sp.symbols('z')
Z0 = sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2      # = 1 + omega,  omega = e^{2 pi i/3}
def chi_psi(A):
    m = len(A)
    negA = [[-x for x in r] for r in A]
    c1 = tlib.chi(negA)                     # det(zI + A)
    AJ = [[-(A[i][j] + 1) for j in range(m)] for i in range(m)]
    c2 = tlib.chi(AJ)                       # det(zI + A + J)
    chi = sp.Poly(c1, z); psi = sp.Poly(c2, z) - chi
    return chi, psi
def R_at(A, z0=Z0):
    chi, psi = chi_psi(A)
    pv = sp.nsimplify(sp.expand(psi.as_expr().subs(z, z0)))
    if sp.simplify(pv) != 0:
        return None, 'z0 not a root of psi'
    cv = sp.expand(chi.as_expr().subs(z, z0))
    if sp.simplify(cv) == 0:
        return None, 'z0 is also a root of chi (fixed eigenvalue)'
    d = sp.expand(sp.diff(psi.as_expr(), z).subs(z, z0))
    if sp.simplify(d) == 0:
        return None, 'z0 multiple root of psi'
    val = sp.radsimp(cv / d)
    return sp.nsimplify(sp.re(sp.expand(val))), val
if __name__ == '__main__':
    A = tlib.from_d6('&FNFL~[udo?')
    R, v = R_at(A); print('FNFL: R* =', R, '=', float(R) if R is not None else None, ' > 1/8:', R is not None and R > sp.Rational(1, 8))
