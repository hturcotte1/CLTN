"""
multiclass.py -- several SINK clique classes (each sends to every other node; classes mutually bidirected).
Class j: size rho_j * k, in-set S_j (i -> class iff S_j[i] = 1), b_j = 1 - S_j.
Lemma X1' (proof in notes): n chi_{D'|perp} = z^(p-1) Phi,  Phi = z psi_D + p chi_D + 1^T adj(zI - D) B,
p = sum rho_j k,  B = sum rho_j k b_j.   Lemma X2': positivity <=> w_G > 0 and sigma_{S_j} > 1 for every j (k-free).
"""
import sympy as sp
import tlib, classext, family_prover as fp
z, k, x = sp.symbols('z k x')
def build(A, classes, kv):
    """classes: list of (S, rho).  Explicit graph with class sizes rho*kv."""
    m = len(A); sizes = [rho * kv for S, rho in classes]; n = m + sum(sizes)
    B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
    idx = m
    owner = {}
    for c, (S, rho) in enumerate(classes):
        for _ in range(rho * kv): owner[idx] = c; idx += 1
    for v in range(m, n):
        c = owner[v]; S = classes[c][0]
        for i in range(m): B[i][v] = S[i]; B[v][i] = 1
        for w in range(m, n):
            if w != v: B[v][w] = 1
    return B
def Phi(A, classes):
    m = len(A)
    D = classext.nonedge_digraph(A)
    Dm = sp.Matrix(D)
    R_adj = (z * sp.eye(m) - Dm).adjugate()
    psi = classext.psi_poly(D).as_expr(); chiD = classext.chi_poly(D).as_expr()
    p = sum(rho for S, rho in classes) * k
    Bv = sp.zeros(m, 1)
    for S, rho in classes:
        for i in range(m):
            if not S[i]: Bv[i] += rho * k
    return sp.expand(z * psi + p * chiD + (sp.ones(1, m) * R_adj * Bv)[0])
def check(A, classes, ks=(1, 2, 3)):
    P = Phi(A, classes)
    for kv in ks:
        G = build(A, classes, kv); n = len(G)
        lhs = classext.psi_poly(classext.nonedge_digraph(G))
        p = sum(rho for S, rho in classes) * kv
        rhs = sp.Poly(sp.expand(z ** (p - 1) * P.subs(k, kv)), z)
        if lhs != rhs: return False
    return True
def k0_for(A, classes, bstar):
    """smallest integer k0 such that every Hurwitz minor of (-1)^m Phi_k(-1/bstar - x) is positive for all k >= k0 (exact);
    None if some minor has non-positive leading coefficient."""
    m = len(A); P = Phi(A, classes)
    h = sp.Poly(sp.expand((-1) ** m * P.subs(z, -1 / bstar - x)), x)
    co = [sp.expand(c) for c in h.all_coeffs()]
    mins = [co[0]] + fp.hurwitz_minors(co)
    k0 = 1
    for mn in mins:
        Pk = sp.Poly(mn, k)
        if Pk.is_zero or Pk.LC() <= 0: return None
        if Pk.degree() > 0:
            rr = [r for r in Pk.real_roots()]
            if rr: k0 = max(k0, int(sp.floor(max(rr))) + 1)
    # confirm no root >= k0 and value positive at k0
    for mn in mins:
        Pk = sp.Poly(mn, k)
        assert (Pk.degree() == 0 and Pk.LC() > 0) or (Pk.count_roots(k0, None) == 0 and Pk.eval(k0) > 0)
    return k0
