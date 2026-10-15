"""Independent re-proof of catalogue entries (families G[S,all;k]).
Phi derived from Lemma W push-through: Phi(z,k) = mu^T adj(zI - D_+ M) 1 = det(X + 1 mu^T) - det(X), X = zI - D_+ M.
Validated against direct compression chi_{D'|perp} = det(zI - P D')/z (P = I - J/n) on explicit graphs.
Positivity of Hurwitz minors on [k0, inf): own Sturm sequences over Q (exact), plus Descartes on p(k0+u).
No floating point in any decision."""
import sys, json, random, time
from fractions import Fraction as F
import sympy as sp
z, k, x, u = sp.symbols('z k x u')
CAT = '/home/user/CLTN/notes/research_wf3/theory/'

def d6_decode(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits += [(v >> (5 - i)) & 1 for i in range(6)]
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

def blow(A, S, kk):
    m = len(A); n = m + kk
    B = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j] = A[i][j]
        for v in range(m, n): B[i][v] = S[i]; B[v][i] = 1   # T = all
    for a in range(m, n):
        for b in range(m, n):
            if a != b: B[a][b] = 1
    return B

def Phi_W(A, S):
    m = len(A)
    Ap = blow(A, S, 1)              # G + v
    Dp = sp.Matrix(m + 1, m + 1, lambda i, j: 0 if i == j else 1 - Ap[i][j])
    mu = [1] * m + [k]
    M = sp.diag(*mu)
    X = z * sp.eye(m + 1) - Dp * M
    one_muT = sp.Matrix(m + 1, m + 1, lambda i, j: mu[j])
    return sp.expand(X.add(one_muT).det(method='berkowitz') - X.det(method='berkowitz'))

def direct_nchi_perp_D(B):
    """n * chi_{D'|perp}(z) for explicit graph B, via det(zI - P D) = z chi_perp, P = I - J/n."""
    n = len(B)
    D = sp.Matrix(n, n, lambda i, j: 0 if i == j else 1 - B[i][j])
    P = sp.eye(n) - sp.ones(n, n) / n
    from sympy.polys.matrices import DomainMatrix
    dm = DomainMatrix.from_Matrix(P * D).convert_to(sp.QQ)
    cp = dm.charpoly()      # coefficients hi->lo
    poly = sp.Poly([sp.Rational(c) for c in cp], z)
    q, r = sp.div(poly, sp.Poly(z, z))
    assert r.is_zero
    return sp.expand(n * q.as_expr())

def hurwitz_minors(co):
    d = len(co) - 1
    a = lambda i: co[i] if 0 <= i <= d else 0
    H = sp.Matrix(d, d, lambda i, j: a(2 * (i + 1) - (j + 1)))
    return [sp.factor(H[:j, :j].det(method='bareiss')) for j in range(1, d + 1)]

# ---- exact Sturm counting over Q
def pdiv_rem(a, b):
    a = a[:]
    while len(a) >= len(b) and any(a):
        c = a[0] / b[0]; 
        for i in range(len(b)): a[i] -= c * b[i]
        a.pop(0)
    while a and a[0] == 0: a.pop(0)
    return a

def ev(p, t):
    r = F(0)
    for c in p: r = r * t + c
    return r

def sturm_count_above(p, a0):
    """number of distinct real roots of p in (a0, inf), p a list of Fractions hi->lo, p(a0) != 0."""
    seq = [p]
    dp = [c * (len(p) - 1 - i) for i, c in enumerate(p[:-1])]
    seq.append(dp)
    while len(seq[-1]) > 1 or (len(seq[-1]) == 1 and False):
        r = pdiv_rem(seq[-2], seq[-1])
        if not r: break
        seq.append([-c for c in r])
        if len(r) == 1: break
    def changes(vals):
        vals = [v for v in vals if v != 0]
        return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))
    at_a = [ev(s, a0) for s in seq]
    at_inf = [s[0] for s in seq]
    return changes(at_a) - changes(at_inf)

def positive_on(pexpr, k0):
    P = sp.Poly(pexpr, k)
    if P.is_zero: return False, 'zero'
    co = [F(int(c.p), int(c.q)) for c in [sp.Rational(c) for c in P.all_coeffs()]]
    if len(co) == 1: return co[0] > 0, 'const'
    if co[0] <= 0: return False, 'LC<=0'
    v0 = ev(co, F(k0))
    if v0 <= 0: return False, 'p(k0)<=0'
    sh = sp.Poly(sp.expand(pexpr.subs(k, k0 + u)), u).all_coeffs()
    descartes = all(c >= 0 for c in sh)
    nroots = sturm_count_above(co, F(k0))
    # cross-check with sympy count_roots
    nr2 = P.count_roots(k0, None)
    assert nr2 == nroots, (nr2, nroots)
    return nroots == 0, 'descartes=%s sturm_roots_above=%d' % (descartes, nroots)

def solve_q(Mx, b):
    m = len(Mx); M = [list(r) + [b[i]] for i, r in enumerate(Mx)]
    for c in range(m):
        p = next((r for r in range(c, m) if M[r][c] != 0), None)
        if p is None: return None
        M[c], M[p] = M[p], M[c]; pv = M[c][c]; M[c] = [v / pv for v in M[c]]
        for r in range(m):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [a - f * cc for a, cc in zip(M[r], M[c])]
    return [M[i][m] for i in range(m)]

def w_G(A, beta):
    g = 1 - 1 / F(beta); m = len(A)
    return solve_q([[F(A[j][i]) + (g if i == j else 0) for j in range(m)] for i in range(m)], [F(1)] * m)

def limit_stable(Phi, m, kk, beta):
    """all roots of Phi_kk have Re > -1/beta  (exact Hurwitz on (-1)^m Phi(-1/beta - x)), three-valued."""
    e = sp.expand(Phi.subs(k, kk).subs(z, -1 / sp.Rational(beta) - x) * (-1) ** m)
    co = sp.Poly(e, x).all_coeffs()
    if co[0] <= 0: return None
    if any(c == 0 for c in co): return None
    if any(c < 0 for c in co): return False
    for mn in hurwitz_minors(co):
        if mn == 0: return None
        if mn < 0: return False
    return True

def reprove(entry, k0=1):
    name = entry['base']; d6 = name.split(':')[-1]
    A = d6_decode(d6); S = list(entry['S']); m = len(A)
    out = {'base': name, 'S': S}
    assert not all(A[i][j] for i in range(m) for j in range(m) if i != j)   # base not complete
    Phi = Phi_W(A, S)
    # compare with catalogue's Phi_k
    out['Phi_matches_catalogue'] = sp.expand(Phi - sp.sympify(entry['detail']['Phi_k'])) == 0
    # validate Phi against direct compression on explicit graphs k = 1, 2, 3
    for kk in (1, 2, 3):
        B = blow(A, S, kk)
        lhs = direct_nchi_perp_D(B)
        rhs = sp.expand(z ** (kk - 1) * Phi.subs(k, kk))
        assert sp.expand(lhs - rhs) == 0, ('W1 fails', kk)
    out['Phi_validated_direct_k123'] = True
    h = sp.Poly(sp.expand((-1) ** m * 2 ** m * Phi.subs(z, -sp.Rational(1, 2) - x)), x)
    co = [sp.expand(c) for c in h.all_coeffs()]
    assert h.degree() == m
    checks = []
    ok, why = positive_on(co[0], k0); checks.append(('a0', ok, why))
    mins = hurwitz_minors(co)
    for j, mn in enumerate(mins):
        ok, why = positive_on(mn, k0); checks.append(('Delta%d' % (j + 1), ok, why))
    out['minors'] = [str(mn) for mn in mins]
    out['minor_checks'] = checks
    out['alpha_lt_-1/2_all_k>=%d' % k0] = all(c[1] for c in checks)
    # catalogue's minors equal ours (up to the same normalisation)
    cm = entry['detail']['hurwitz_minors_in_k']
    out['minors_match_catalogue'] = all(sp.expand(sp.sympify(a) - b) == 0 for a, b in zip(cm, mins)) and len(cm) == len(mins)
    w = w_G(A, 2)
    sig = sum(w[i] for i in range(m) if S[i])
    out['w_G(2)'] = [str(v) for v in w]; out['sigma_S(2)'] = str(sig)
    out['positive_at_2'] = w is not None and all(v > 0 for v in w) and sig > 1
    out['PROVED'] = out['alpha_lt_-1/2_all_k>=%d' % k0] and out['positive_at_2']
    return out, A, S, Phi

if __name__ == '__main__':
    cat6 = json.load(open(CAT + 'families_catalogue_n6.json'))
    cat7 = json.load(open(CAT + 'families_catalogue_n7.json'))
    allc = [('n6', i, e) for i, e in enumerate(cat6)] + [('n7', i, e) for i, e in enumerate(cat7)]
    rng = random.Random(20261015)
    pick = rng.sample(allc, 10)
    # also add the unique n7 entry not proved for k >= 0, as a control
    extra = [('n7', i, e) for i, e in enumerate(cat7) if not e['proved_for_all_k>=0']]
    res = []
    for tag, i, e in pick + extra:
        t0 = time.time()
        r, A, S, Phi = reprove(e)
        r['catalogue'] = tag; r['index'] = i; r['catalogue_says_k>=1'] = e['proved_for_all_k>=1']
        if (tag, i, e) in extra:
            r0, _, _, _ = reprove(e, k0=0); r['k>=0_recheck'] = r0['alpha_lt_-1/2_all_k>=0']
            r['k>=0_detail'] = [c for c in r0['minor_checks'] if not c[1]]
        r['secs'] = round(time.time() - t0, 1)
        print(tag, i, e['base'], e['S'], 'PROVED' if r['PROVED'] else 'FAILED', r['Phi_matches_catalogue'], r['minors_match_catalogue'], r['secs'], flush=True)
        res.append(r)
    json.dump(res, open('reprove_results.json', 'w'), indent=1)
