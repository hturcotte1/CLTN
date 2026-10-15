"""Independent re-verification of every exact instance record in the theory certificate files.
Own digraph6 decoder (cross-checked against n7/curve_test.from_d6), own integer Berkowitz charpoly, own Hurwitz
leading-minor computation (exact Fractions, Gaussian elimination without pivoting; a zero pivot = undecided = FAIL).
Finite point: ctln.analyze_support (repo verifier) AND, for n <= 20, my own charpoly/Hurwitz of -I+W as a cross-check.
Limit claim at beta: every root of psi_A(z) = 1^T adj(zI - A) 1 has Re < -gamma, gamma = 1 - 1/beta.
Boundary (finite-eps-only) claim: (z^2+z+1) divides psi_A exactly once, cofactor has every root with Re < -1/2.
No floating point anywhere.  Usage: python3 verify_instances.py [maxn]"""
import sys, json, os
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/n7')
import ctln
from curve_test import from_d6 as repo_from_d6
TH = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/theory/'
MAXN = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

def my_d6(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63
    assert s[1] != '~' and 0 < n
    if n >= 63: print('  NOTE: non-standard digraph6 header chr(%d) for n=%d (standard form needs ~ + 18 bits)' % (ord(s[1]), n))
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits += [(v >> (5 - i)) & 1 for i in range(6)]
    assert len(bits) >= n * n and not any(bits[n * n:])
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

def berkowitz(M):
    """charpoly det(zI - M), hi->lo, exact (ints or Fractions)."""
    n = len(M)
    C = [1]
    for k in range(n):        # add row/col k
        # leading principal submatrix of size k, column R = M[0..k-1][k], row S = M[k][0..k-1], a = M[k][k]
        a = M[k][k]
        R = [M[i][k] for i in range(k)]; S = [M[k][j] for j in range(k)]
        # Toeplitz column: 1, -a, -S R, -S A R, -S A^2 R, ...
        col = [1, -a]
        v = R[:]
        for p in range(k):
            col.append(-sum(S[j] * v[j] for j in range(k)))
            v = [sum(M[i][j] * v[j] for j in range(k)) for i in range(k)]
        # new C = T * C, T lower-triangular Toeplitz (k+2) x (k+1)
        newC = []
        for i in range(k + 2):
            newC.append(sum(col[i - j] * C[j] for j in range(max(0, i - (k + 1)), min(i, k) + 1) if 0 <= i - j < len(col)))
        C = newC
    return C

def polysub(a, b):
    L = max(len(a), len(b)); a = [0] * (L - len(a)) + a; b = [0] * (L - len(b)) + b
    r = [x - y for x, y in zip(a, b)]
    while r and r[0] == 0: r = r[1:]
    return r

def psi_A(A):
    n = len(A)
    AmJ = [[A[i][j] - 1 for j in range(n)] for i in range(n)]
    return polysub(berkowitz(AmJ), berkowitz(A))   # det(zI-A+J) - det(zI-A) = 1^T adj(zI-A) 1

def shift(p, c):
    """coefficients of p(x + c), hi->lo."""
    p = [F(v) for v in p]; n = len(p) - 1
    out = p[:]
    for i in range(n):
        for j in range(1, n - i + 1):
            out[j] += c * out[j - 1]
    return out

def hurwitz_all_positive(p):
    """True iff leading coeff > 0 and all leading Hurwitz minors > 0 (exact). False/None otherwise (None: zero minor)."""
    p = [F(v) for v in p]
    if p[0] < 0: p = [-v for v in p]
    d = len(p) - 1
    if d == 0: return True
    a = lambda i: p[i] if 0 <= i <= d else F(0)
    H = [[a(2 * (i + 1) - (j + 1)) for j in range(d)] for i in range(d)]
    for k in range(d):
        piv = H[k][k]
        if piv == 0: return None
        if piv < 0: return False     # leading minors' ratio negative -> some minor negative (product of pivots)
        for i in range(k + 1, d):
            f = H[i][k] / piv
            if f:
                for j in range(k, d): H[i][j] -= f * H[k][j]
    return True

def polydiv(a, b):
    a = [F(v) for v in a]; q = []
    while len(a) >= len(b):
        c = a[0] / b[0]; q.append(c)
        a = [x - c * y for x, y in zip(a, b + [0] * (len(a) - len(b)))][1:]
    while a and a[0] == 0: a = a[1:]
    return q, a

def limit_ok(A, beta):
    gamma = 1 - 1 / F(beta)
    return hurwitz_all_positive(shift(psi_A(A), -gamma))   # roots x = lambda + gamma

def boundary_ok(A):
    p = psi_A(A); phi3 = [1, 1, 1]
    q, r = polydiv(p, phi3)
    if r: return False, 'phi3 does not divide'
    q2, r2 = polydiv(q, phi3)
    if not r2: return False, 'phi3^2 divides'
    return hurwitz_all_positive(shift(q, F(-1, 2))) is True, 'ok'

def own_point_check(A, eps, delta):
    n = len(A)
    W = ctln.weight_matrix(A, eps, delta)
    M = [[(F(-1) if i == j else F(0)) + W[i][j] for j in range(n)] for i in range(n)]
    cp = berkowitz(M)
    return hurwitz_all_positive(cp)

def check_point(d6, eps, delta, beta=None, q=None, t=None, label=''):
    A = my_d6(d6); assert A == repo_from_d6(d6), 'decoder mismatch ' + d6
    n = len(A)
    eps = F(eps); delta = F(delta)
    errs = []
    if not (eps > 0 and delta > 0 and delta * (1 - eps) / eps > 1): errs.append('illegal')
    if q is not None and F(q) != delta * (1 - eps) / eps: errs.append('q mismatch')
    if beta is not None and F(beta) != 1 + delta / eps: errs.append('beta mismatch %s vs %s' % (beta, 1 + delta / eps))
    if t is not None and F(t) != (1 - eps) / eps: errs.append('t mismatch')
    if all(A[i][j] and A[j][i] for i in range(n) for j in range(n) if i != j): errs.append('is clique')
    r = ctln.analyze_support(A, list(range(n)), eps, delta)
    if not (r['positive'] and r['stable'] is True and r['stable_fixed_point']): errs.append('ctln: pos=%s stable=%s' % (r['positive'], r['stable']))
    if n <= 20:
        own = own_point_check(A, eps, delta)
        if own is not True: errs.append('own Hurwitz check: %s' % own)
    return A, errs

results = []
def rec(label, A_errs, extra=None):
    A, errs = A_errs
    if extra: errs = errs + extra
    results.append((label, len(A), errs))
    print(('OK  ' if not errs else 'FAIL') + ' n=%d %s %s' % (len(A), label, errs), flush=True)

J = lambda f: json.load(open(TH + f))
# F certificate witnesses (limit claimed at beta)
for w in J('family_F_certificate.json')['F6_witnesses']:
    if w['n'] > MAXN or MAXN >= 30: continue
    Ae = check_point(w['d6'], w['eps'], w['delta'], w['beta'], w['q'], w['t'])
    rec('F cert k=%d' % w['k'], Ae, [] if limit_ok(Ae[0], w['beta']) is True else ['limit fails at beta'])
for w in J('F_family_witnesses.json'):
    if w['n'] > MAXN or MAXN >= 30: continue
    Ae = check_point(w['d6'], w['eps'], w['delta'], w['beta'], w['q'], w['t'])
    ex = []
    if w.get('limit_ok_at_beta_exact'):
        if limit_ok(Ae[0], w['beta']) is not True: ex = ['limit fails at beta']
    rec('F early %s k=%d' % (w['family'][:20], w['k']), Ae, ex)
for r in J('family_E_witnesses.json'):
    w = r['witness']
    if r['n'] > MAXN or MAXN >= 30: continue
    Ae = check_point(r['d6'], w['eps'], w['delta'], w['beta'], w['q'], w['t'])
    ok, why = boundary_ok(Ae[0])
    ex = [] if ok else ['boundary: ' + why]
    if limit_ok(Ae[0], w['beta']) is True: ex.append('unexpectedly limit-stable at witness beta')
    rec('E a=%d' % r['a'], Ae, ex)
for r in (J("cube_n8_witnesses.json") if MAXN < 30 else []):
    w = r['witness']
    Ae = check_point(r['d6'], w['eps'], w['delta'], w['beta'], w['q'], w['t'])
    ok, why = boundary_ok(Ae[0])
    rec('cube ' + r['d6'], Ae, [] if ok else ['boundary: ' + why])
for f in (["top_n8.json", "cert_q_910.json", "cert_q_11.json"] if MAXN < 30 else []):
    for r in J(f):
        if 'eps' not in r: continue
        b = r.get('beta_exact', r.get('beta'))
        Ae = check_point(r['d6'], r['eps'], r['delta'], b, r['q'], r['t'])
        rec(f + ' ' + r['d6'], Ae, [] if limit_ok(Ae[0], b) is True else ['limit fails at beta'])
for f, key in ([('family_Y_certificate.json', 'Y5_witness')] if len(sys.argv) > 2 else [('family_Q2_certificate.json', 'W5_witnesses'), ('family_Q_certificate.json', 'Q5_witness'), ('family_Y_certificate.json', 'Y5_witness')]):
    ws = J(f)[key]
    if isinstance(ws, dict): ws = [ws]
    for w in ws:
        if w['n'] > MAXN or (MAXN >= 30 and w['n'] <= 20): continue
        Ae = check_point(w['d6'], w['eps'], w['delta'], w['beta'], w['q'], w['t'])
        rec(f + ' n=%d' % w['n'], Ae, [] if limit_ok(Ae[0], w['beta']) is True else ['limit fails at beta'])
nf = sum(1 for r in results if r[2])
print('TOTAL %d records, %d failures' % (len(results), nf))
