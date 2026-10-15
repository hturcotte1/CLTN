"""Weighted (node blow-up) version.  Blow node a of base graph H (m nodes) into a bidirected clique class of
size mu_a (closed twins).  Lemma W (proof in notes): with D = J - I - A_H, M = diag(mu), n = sum mu,
    n chi_{D'|perp}(z) = z^(n-m) * Phi_mu(z),   Phi_mu(z) := mu^T adj(zI - D M) 1,
and positivity: w' is class-constant, w'|class a = y_a with ((A^T + I) M + (gamma - 1) I) y = 1.
Checked exactly against explicit blow-ups below."""
import sympy as sp, itertools
import tlib, classext
z = sp.symbols('z')
def blowup(H, mu):
    m = len(H); idx = []
    for a in range(m): idx += [a] * mu[a]
    n = len(idx)
    return [[0 if i == j else (1 if idx[i] == idx[j] else H[idx[i]][idx[j]]) for j in range(n)] for i in range(n)]
def Phi_mu(H, mu):
    m = len(H)
    D = sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - H[i][j])
    M = sp.diag(*mu)
    X = z * sp.eye(m) - D * M
    return sp.Poly(sp.expand((sp.Matrix([mu]) * X.adjugate() * sp.ones(m, 1))[0]), z)
def check(H, mus):
    for mu in mus:
        B = blowup(H, mu); n = len(B); m = len(H)
        lhs = sp.Poly(tlib.chi_perp_times_m([[0 if i == j else 1 - B[i][j] for j in range(n)] for i in range(n)]) if False else [0], z)
        Dp = classext.nonedge_digraph(B)
        lhs = classext.psi_poly(Dp)
        rhs = sp.Poly(z ** (n - m), z) * Phi_mu(H, mu)
        if lhs != rhs: return False, mu
    return True, None
if __name__ == '__main__':
    import load_hits
    H6 = load_hits.n6_hits()
    for H in H6.values():
        print(check(H, [(1,1,1,1,1,1), (2,1,1,1,1,1), (1,3,1,2,1,1), (2,2,1,1,1,3)]))
