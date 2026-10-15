"""Integer weight vectors mu (node blow-up sizes) of a base graph H for which the cube root omega is a zero of
Phi_mu (weighted.py), i.e. alpha(A'|perp) has the eigenvalue with real part exactly -1/2 (D-picture zeta = omega)."""
import sympy as sp, tlib, sys, itertools
z = sp.symbols('z')
def omega_UV(H):
    m = len(H)
    mus = sp.symbols('m0:%d' % m)
    D = sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - H[i][j])
    X = z * sp.eye(m) - D * sp.diag(*mus)
    Phi = sp.expand((sp.Matrix([mus]) * X.adjugate() * sp.ones(m, 1))[0])
    rem = sp.rem(sp.Poly(Phi, z), sp.Poly(z**2 + z + 1, z))
    co = rem.all_coeffs()
    V = co[0] if len(co) == 2 else sp.Integer(0); U = co[-1]
    return mus, sp.expand(U), sp.expand(V)
if __name__ == '__main__':
    A = tlib.from_d6('&FNFL~[udo?')
    keep = [0,1,3,4,5,6]
    H = [[A[i][j] for j in keep] for i in keep]
    mus, U, V = omega_UV(H)
    print('U =', sp.factor(U)); print('V =', sp.factor(V))
    fU = sp.lambdify(mus, U, 'math'); fV = sp.lambdify(mus, V, 'math')   # integer inputs -> exact python ints
    sols = []
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for mu in itertools.product(range(1, B + 1), repeat=6):
        if fU(*mu) == 0 and fV(*mu) == 0: sols.append(mu)
    print(len(sols), sols[:50])
