"""
family_prover.py -- automatic EXACT proofs of infinite clique-class families G[S, all; k].

For a base digraph G (m nodes) and an in-set S (T = all: the class sends to every node of G), Lemma X1 gives
    n chi_{D'|perp}(z) = z^(k-1) Phi_k(z),   Phi_k = z psi_D + k chi_L  (L = D - b 1^T, b = 1 - S),
and Lemma X2 gives  Pos(G[S,all;k], beta)  <=>  w_G(beta) > 0  and  sigma_S(beta) = sum_{i in S} w_G,i(beta) > 1  (k-free).
alpha(G[S,all;k]|perp) < -1/2  <=>  every root of Phi_k has Re > -1/2  <=>  h_k(x) := (-1)^m Phi_k(-1/2 - x) Hurwitz.
The Hurwitz minors of h_k are polynomials in k; if each is positive on [k0, inf) (exact real-root counting), then
alpha < -1/2 for every integer k >= k0.  With positivity strictly at beta = 2 (exact), continuity gives positivity on
(2, 2 + eta), so for every k >= k0 there is beta in (2, min(beta_c(k), 2 + eta)) where the full support is an SFP for all
large t (Lemma I.3 combined statement).  This proves an infinite family.  No floating point in any decision.
"""
import sys, json
import sympy as sp
from fractions import Fraction as F
import tlib, classext, ctln
z, k, x = sp.symbols('z k x')

def hurwitz_minors(coeffs):
    d = len(coeffs) - 1
    a = lambda i: coeffs[i] if 0 <= i <= d else 0
    H = sp.Matrix(d, d, lambda i, j: a(2 * (i + 1) - (j + 1)))
    return [sp.expand(H[:j, :j].det(method='berkowitz')) for j in range(1, d + 1)]

def positive_on(p, k0):
    """exact: polynomial p(k) > 0 for all real k >= k0."""
    P = sp.Poly(p, k)
    if P.is_zero: return False
    if P.degree() == 0: return P.LC() > 0
    if P.eval(k0) <= 0: return False
    # roots in [k0, inf): count real roots > k0 via sturm on (k0, big) and leading coeff sign
    n_roots = len([r for r in P.real_roots() if r >= k0])
    return n_roots == 0 and P.LC() > 0

def prove_family(A, S, k0=1):
    m = len(A); T = (1,) * m
    Q, P = classext.Phi(A, S, T)
    Phi = sp.expand(Q.as_expr() + k * P.as_expr())
    h = sp.Poly(sp.expand((-1) ** m * Phi.subs(z, -sp.Rational(1, 2) - x) * 2 ** m), x)
    co = [sp.expand(c) for c in h.all_coeffs()]
    if not positive_on(co[0], k0): return {'ok': False, 'why': 'leading coefficient'}
    mins = hurwitz_minors(co)
    for j, mn in enumerate(mins):
        if not positive_on(mn, k0):
            return {'ok': False, 'why': 'Hurwitz minor %d not positive on [%d, inf)' % (j + 1, k0)}
    # positivity at beta = 2 (gamma = 1/2), exact
    w = tlib.w_vec(A, F(2))
    if w is None or not all(v > 0 for v in w): return {'ok': False, 'why': 'w_G(2) not > 0'}
    sig = sum(w[i] for i in range(m) if S[i])
    if not sig > 1: return {'ok': False, 'why': 'sigma_S(2) <= 1'}
    Pl = sp.factor(P.as_expr())
    return {'ok': True, 'Phi_k': str(sp.factor(Phi)), 'chi_L': str(Pl), 'z_psi_D': str(sp.factor(Q.as_expr())),
            'hurwitz_minors_in_k': [str(sp.factor(mn)) for mn in mins], 'sigma_S_at_2': str(sig), 'w_G_at_2': [str(v) for v in w]}

if __name__ == '__main__':
    pass
