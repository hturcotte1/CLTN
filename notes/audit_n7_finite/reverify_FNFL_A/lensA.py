"""Independent Lens-A re-verification of certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json.
Builds W from CLAUDE.md sec. 2 directly in sympy; exact rationals only."""
import json, sys
import sympy as sp
from sympy import Rational as R, Matrix, symbols, Poly, eye

CERT = '/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'
d = json.load(open(CERT))['graph']
ok_all = True
def check(name, cond):
    global ok_all
    print(('PASS ' if cond else 'FAIL ') + name)
    if not cond: ok_all = False

# --- adjacency from rows
rows = d['adjacency_rows']
A = [[int(c) for c in r] for r in rows]
n = len(A)
check('n == 7', n == 7)
check('zero diagonal', all(A[i][i] == 0 for i in range(n)))

# --- adjacency from digraph6 (own decoder): '&', N(n)=chr(63+n), then n*n bits row-major, 6 per char, MSB first
s = d['digraph6_as_enumerated']
assert s[0] == '&'
n6 = ord(s[1]) - 63
bits = []
for ch in s[2:]:
    v = ord(ch) - 63
    bits += [(v >> (5 - b)) & 1 for b in range(6)]
A6 = [[bits[i * n6 + j] for j in range(n6)] for i in range(n6)]
check('digraph6 decodes to n=7', n6 == 7)
check('digraph6 adjacency == adjacency_rows', A6 == A)
check('padding bits zero', all(b == 0 for b in bits[n6 * n6:]))
E = sum(map(sum, A)); check(f'edge count {E} == cert {d["edges"]}', E == d['edges'])

# --- graph notions
clique = all(A[i][j] == 1 and A[j][i] == 1 for i in range(n) for j in range(n) if i != j)
check('not a clique', not clique)
missing = [(i, j) for i in range(n) for j in range(n) if i != j and A[i][j] == 0]
print('  non-edges:', len(missing), missing[:6], '...')
# own Def 1.1 (Curto): k dominates j iff (i) all i notin{j,k}: i->j => i->k ; (ii) j->k and not k->j
def dom(k, j):
    if not (A[j][k] == 1 and A[k][j] == 0): return False
    return all(not (A[i][j] == 1 and A[i][k] == 0) for i in range(n) if i not in (j, k))
mydom = [(k, j) for j in range(n) for k in range(n) if k != j and dom(k, j)]
check(f'domination-free (own Lemma D reading): dominating pairs={mydom}', mydom == [])
sys.path.insert(0, '/home/user/CLTN')
import graphs, ctln
check('graphs.is_domination_free agrees', graphs.is_domination_free(A) is True)
# Also: Lemma D' inside sigma = full set is the same condition (sigma = V).
check('full support: no node outside sigma (condition 2 vacuous, no targets possible)', sorted(range(n)) == list(range(7)))

# --- symbolic W
e, dl = symbols('eps delta', positive=True)
def Wsym():
    return Matrix(n, n, lambda i, j: 0 if i == j else (-1 + e if A[j][i] == 1 else -1 - dl))
Ws = Wsym()
# closed form cross-check
Am = Matrix(A); J = sp.ones(n, n)
closed = -dl * eye(n) + (1 + dl) * J - (e + dl) * Am.T   # entry (i,j) of I-W reads A[j][i]
check('I - W == -delta I + (1+delta) J - (eps+delta) A^T  (symbolic identity)', sp.simplify(eye(n) - Ws - closed) == sp.zeros(n, n))
# NOTE: CLAUDE.md writes A_sigma; with W[i][j] reading A[j][i] the matrix is A^T entrywise. Record whether A vs A^T matters:
closedA = -dl * eye(n) + (1 + dl) * J - (e + dl) * Am
print('  closed form with A (not A^T) equals I-W?', sp.simplify(eye(n) - Ws - closedA) == sp.zeros(n, n))

def hurwitz_minors(c):
    """c = [a0,...,am] coefficients high->low of monic-ish poly; return [Delta_1..Delta_m]"""
    m = len(c) - 1
    a = lambda k: c[k] if 0 <= k <= m else 0
    H = Matrix(m, m, lambda i, j: a(2 * (j + 1) - (i + 1)))
    return [H[:k, :k].det() for k in range(1, m + 1)]

def rhp_check(cp):
    """independent stability check: exact count of roots in closed right half-plane box
    [0,B]x[-B,B] (B = Cauchy bound) via sympy's exact complex root counting; returns (left,right)."""
    lam = symbols('lam')
    P = Poly(cp, lam)
    c = P.all_coeffs()
    B = 1 + max(abs(R(ci) / R(c[0])) for ci in c[1:])
    B = sp.ceiling(B) + 1
    right = P.count_roots(0 - B*sp.I, B + B*sp.I)
    total = P.degree()
    return (total - right, right)

def analyze(eps, delta, label):
    print(f'\n=== {label}: eps={eps} delta={delta}')
    t = (1 - eps) / eps; beta = 1 + delta / eps; q = delta * (1 - eps) / eps
    legal = eps > 0 and delta > 0 and q > 1
    print(f'  t={t} beta={beta} q={q} legal={legal}  t(beta-2)={t*(beta-2)}')
    W = Ws.subs({e: eps, dl: delta})
    IW = eye(n) - W
    x = IW.LUsolve(sp.ones(n, 1))
    x = [sp.cancel(sp.together(v)) for v in x]
    assert all(v.is_Rational for v in x)
    resid = IW * Matrix(x) - sp.ones(n, 1)
    fp = all(r == 0 for r in resid)
    pos = all(v > 0 for v in x)
    M = -eye(n) + W
    lam = symbols('lam')
    cp = (lam * eye(n) - M).det(method='berkowitz')
    cp = Poly(sp.expand(cp), lam)
    coeffs = [R(c) for c in cp.all_coeffs()]
    Hs = hurwitz_minors(coeffs)
    signs = ['+' if h > 0 else ('0' if h == 0 else '-') for h in Hs]
    allcoefpos = all(c > 0 for c in coeffs)
    stable = allcoefpos and all(h > 0 for h in Hs)
    print('  x =', [str(v) for v in x]); print('  x positive:', pos, ' residual zero:', fp)
    print('  charpoly:', [str(c) for c in coeffs])
    print('  Hurwitz signs:', signs, ' stable(all>0):', stable)
    for k, h in enumerate(Hs, 1):
        if h <= 0: print(f'   Delta_{k} = {h}')
    # independent: root isolation
    rr = rhp_check(cp.as_expr())
    print('  root isolation (left,right):', rr)
    # ctln cross-check (Fractions)
    from fractions import Fraction as F
    rec = ctln.analyze_support(A, list(range(n)), F(eps.p, eps.q), F(delta.p, delta.q))
    cx = [R(str(v)) for v in rec['x']] if rec.get('x') is not None else None
    print('  ctln: stable=', rec.get('stable'), '|', rec.get('stability_detail'), '| x match:', cx == x)
    return dict(legal=legal, fp=fp, pos=pos, coeffs=coeffs, Hs=Hs, signs=signs, stable=stable, rr=rr, rec=rec, x=x, cx=cx, beta=beta, t=t, q=q)

wit = d['witnesses']
out = {}
for (eps, delta, lab, want_stable) in [(R(1, 49), R(33, 1568), 'W1 stable', True),
                                       (R(1, 36), R(17, 588), 'W2 stable', True),
                                       (R(1, 101), R(33, 3232), 'W3 unstable same beta', False)]:
    r = analyze(eps, delta, lab); out[lab] = r
    check(f'{lab}: legal', r['legal'])
    check(f'{lab}: fixed-point residual 0', r['fp'])
    check(f'{lab}: x > 0', r['pos'])
    if want_stable:
        check(f'{lab}: all charpoly coeffs >0 and all Hurwitz minors >0', r['stable'])
        check(f'{lab}: root isolation: 7 roots in open LHP', r['rr'] == (7, 0))
    else:
        check(f'{lab}: NOT stable', not r['stable'])
        check(f'{lab}: no vanishing minor (strictly decided)', all(h != 0 for h in r['Hs']) and all(c != 0 for c in r['coeffs']))
        check(f'{lab}: some minor strictly negative', any(h < 0 for h in r['Hs']))
        check(f'{lab}: root isolation finds a root with Re>0 or not all left', r['rr'] is not None and r['rr'] != (7, 0))
    check(f'{lab}: ctln.analyze_support stable == {want_stable}', r['rec'].get('stable') is want_stable)
    check(f'{lab}: ctln x == sympy x', r['cx'] == r['x'])
    # compare with certificate witness
    pool = wit['stable_positive_legal'] + wit['unstable_same_beta_larger_t']
    cw = [w for w in pool if R(w['eps']) == eps and R(w['delta']) == delta]
    check(f'{lab}: certificate witness present', len(cw) == 1)
    if cw:
        cw = cw[0]
        check(f'{lab}: cert x matches', [R(v) for v in cw['x']] == r['x'])
        check(f'{lab}: cert charpoly matches', [R(v) for v in cw['charpoly']] == r['coeffs'])
        check(f'{lab}: cert t,beta,q match', (R(cw['t']), R(cw['beta']), R(cw['q'])) == (r['t'], r['beta'], r['q']))
        check(f'{lab}: cert stable_detail == ctln detail', cw['stable_detail'] == r['rec'].get('stability_detail'))
check('W1 and W3 share beta', out['W1 stable']['beta'] == out['W3 unstable same beta']['beta'])
print('\nALL LENS-A CHECKS PASS' if ok_all else '\nSOME CHECK FAILED')
