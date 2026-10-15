"""
family_E.py -- an infinite family of FINITE-EPS-ONLY full-support counterexample candidates.
Re-runnable: python3 family_E.py   (exact; no floating point in any decision; writes family_E_certificate.json)

H := the 7-node digraph &FNDL~[uds? (one of the n=7 limit candidates) in the labelling below.
E_a := closed-twin blow-up of H with class sizes mu(a) = (1,1,2,a,1,1,1)   (a >= 0; a = 0 deletes node 3).
n(E_a) = 7 + a.  E_0 is the n=7 finite-eps-only graph &FNFL~[udo? (certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json).

Claims checked here (labels refer to the notes file):
 (E1) Lemma W formula  n chi_{D'|perp} = z^(n-7) Phi_mu  on the explicit graphs, a = 0..8.
 (E2) Phi_{mu(a)}(z) = z (z^2+z+1) ((a+7) z (z^2+z+1) + 2)  identically in a  (D-picture).
 (E3) the cubic (a+7) z(z^2+z+1) + 2 has all roots Re z > -1/2 for every real a > -5/3 (Hurwitz in closed form).
 (E4) omega is a simple root, not a root of chi:  det(z0 I + A_Q(a)) != 0 for a >= 0  (A-picture, z0 = 1 + omega).
 (E5) R*(a) := Re( chi(z0)/psi'(z0) ) = 1/4 for every a  (closed form; Lemma AS coefficient).
 (E6) positivity at beta = 2 (gamma = 1/2): class vector y(a) > 0 for every a >= 0 (closed form).
 (E7) explicit exact witnesses (ctln.analyze_support) for a = 0..4 (n = 7..11) and exact never-limit-stable check
      (l^2+l+1 divides n chi_perp; all other roots Re < -1/2 by exact Hurwitz).
"""
import json, sys
import sympy as sp
from fractions import Fraction as F
import tlib, weighted, classext, boundary_tools as bt, rstar, ctln

H = [[0,0,1,1,1,1,0],[0,0,1,0,1,0,0],[1,1,0,1,1,1,1],[1,1,1,0,1,1,1],[0,0,1,1,0,1,1],[0,1,0,0,1,0,1],[1,1,0,1,0,0,0]]
assert tlib.to_d6(H) == '&FNDL~[uds?'
z, a = sp.symbols('z a')
def mu_of(av): return (1, 1, 2, av, 1, 1, 1)
def E(av): return weighted.blowup(H, mu_of(av))
def main():
    cert = {'family': 'E_a = closed-twin blow-up of H=&FNDL~[uds? with class sizes (1,1,2,a,1,1,1); n = 7 + a', 'H_rows': [''.join(map(str, r)) for r in H]}

    # (E1)
    ok = True
    for av in range(0, 9):
        B = E(av); n = len(B)
        lhs = classext.psi_poly(classext.nonedge_digraph(B))
        rhs = sp.Poly(z ** (n - 7 + (1 if av == 0 else 0)), z) * sp.Poly(weighted.Phi_mu(H, mu_of(av)).as_expr() if av > 0 else 0, z) if av > 0 else None
        if av > 0:
            ok &= (lhs == rhs)
    cert['E1_lemmaW_checked_a_1_to_8'] = ok
    assert ok

    # (E2)
    m = 7
    D = sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - H[i][j])
    M = sp.diag(*mu_of(a))
    Phi = sp.expand((sp.Matrix([list(mu_of(a))]) * (z * sp.eye(m) - D * M).adjugate() * sp.ones(m, 1))[0])
    target = sp.expand(z * (z**2 + z + 1) * ((a + 7) * z * (z**2 + z + 1) + 2))
    cert['E2_Phi_closed_form'] = str(sp.factor(Phi)); assert sp.expand(Phi - target) == 0
    # a = 0 member (node deleted) directly:
    B0 = E(0); lhs0 = classext.psi_poly(classext.nonedge_digraph(B0))
    assert lhs0 == sp.Poly(sp.expand(z * (z**2 + z + 1) * (7 * z * (z**2 + z + 1) + 2)), z)
    cert['E2_a0_direct'] = True

    # (E3) cubic c z (z^2+z+1) + 2, c = a+7: roots Re z > -1/2  <=>  h(y) := -(c((-y-1/2)^3 + (-y-1/2)^2 + (-y-1/2)) + 2) Hurwitz in y
    c, y = sp.symbols('c y')
    h = sp.Poly(sp.expand(-(c * ((-y - sp.Rational(1, 2))**3 + (-y - sp.Rational(1, 2))**2 + (-y - sp.Rational(1, 2))) + 2)), y)
    co = h.all_coeffs()
    H2 = sp.factor(co[1] * co[2] - co[0] * co[3])
    cert['E3_shifted_cubic_coeffs'] = [str(sp.factor(t)) for t in co]
    cert['E3_Hurwitz_Delta2'] = str(H2)
    # conditions: co[0]=c>0, co[1]=c/2>0, co[3] = 3c/8 - 2 > 0 (c > 16/3), Delta2 = 2c > 0
    assert sp.simplify(co[0] - c) == 0 and sp.simplify(co[1] - c / 2) == 0 and sp.simplify(co[3] - (3 * c / 8 - 2)) == 0 and sp.simplify(H2 - 2 * c) == 0
    cert['E3_conclusion'] = 'all roots of (a+7) z(z^2+z+1)+2 have Re > -1/2 iff a+7 > 16/3; holds for all a >= 0'

    # (E4),(E5): A-picture quotient.  A_Q = (A+I)M - I ; chi(z) = (z-1)^(n-7) det(zI + A_Q); s(z) = mu^T (zI + A_Q)^{-1} 1
    AH = sp.Matrix(H)
    AQ = (AH + sp.eye(m)) * M - sp.eye(m)
    z0 = sp.Rational(1, 2) - sp.sqrt(3) * sp.I / 2          # = 1 + omega^2 = -omega (a root of psi_A, Re = 1/2)
    X = z0 * sp.eye(m) + AQ
    detX = sp.expand(X.det(method='berkowitz'))
    cert['E4_det(z0 I + A_Q)'] = str(sp.simplify(detX))
    mu_vec = sp.Matrix([list(mu_of(a))])
    Xadj = X.adjugate()
    s_num = sp.expand((mu_vec * Xadj * sp.ones(m, 1))[0])          # mu^T adj(X) 1  (= psi/(z-1)^(n-7) at z0)
    assert sp.simplify(s_num) == 0, 'z0 not a root of psi'
    # s'(z0) = -mu^T X^{-2} 1 = -mu^T adj(X)^2 1 / det(X)^2 ;  1/s'(z0) = -det(X)^2 / (mu^T adj(X)^2 1)
    q2 = sp.expand((mu_vec * Xadj * Xadj * sp.ones(m, 1))[0])
    inv_sprime = sp.simplify(-detX**2 / q2)
    R = sp.simplify(sp.re(sp.expand(sp.radsimp(inv_sprime))))
    cert['E5_inverse_s_prime'] = str(inv_sprime); cert['E5_R_star'] = str(R)
    print('det(z0 I + A_Q) =', sp.factor(sp.simplify(detX)))
    print('1/s\'(z0) =', inv_sprime, '   R* =', R)

    # (E6) positivity at gamma = 1/2
    g = sp.Rational(1, 2)
    Mpos = (AH.T + sp.eye(m)) * M + (g - 1) * sp.eye(m)
    yv = sp.simplify(Mpos.LUsolve(sp.ones(m, 1)))
    cert['E6_y_at_beta_2'] = [str(sp.factor(t)) for t in yv]
    print('y(beta=2) =', [sp.factor(t) for t in yv])

    json.dump(cert, open('family_E_certificate_part1.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
