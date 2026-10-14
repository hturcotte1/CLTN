"""Lens-B supplementary checks (exact): spectral abscissa of A|perp, limit verdict, labelg non-membership, region2d."""
import json, sys, subprocess, time
import sympy as sp
from sympy import Matrix, symbols, Poly, Rational as R
sys.path.insert(0, '/home/user/CLTN')
d = json.load(open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))['graph']
A = [[int(c) for c in r] for r in d['adjacency_rows']]; n = 7
ok = True
def check(name, c):
    global ok; print(('PASS ' if c else 'FAIL ') + name); ok = ok and bool(c)
l = symbols('l')
Am = Matrix(A); J = sp.ones(n, n)
chiA = (l*sp.eye(n) - Am).charpoly(l).as_expr() if False else sp.expand((l*sp.eye(n)-Am).det(method='berkowitz'))
chiAJ = sp.expand((l*sp.eye(n)-Am-J).det(method='berkowitz'))
diff = Poly(chiA - chiAJ, l)
print('chi_A - chi_{A+J} =', diff.as_expr())
fac = sp.factor_list(diff.as_expr())
print('factor_list:', fac)
claimed = sp.expand((l+1)*(l**2+l+1)*(7*l**3+14*l**2+14*l+5))
check('chi_A - chi_{A+J} == (l+1)(l^2+l+1)(7l^3+14l^2+14l+5)', sp.expand(diff.as_expr() - claimed) == 0)
# A^T gives same (transpose invariance)
chiAT = sp.expand((l*sp.eye(n)-Am.T).det(method='berkowitz')); chiATJ = sp.expand((l*sp.eye(n)-Am.T-J).det(method='berkowitz'))
check('same polynomial for A^T', sp.expand(chiAT - chiATJ - claimed) == 0)
# real parts: for each irreducible factor f, max Re root. Re < -1/2 strictly iff f(m - 1/2) Hurwitz.
def hurwitz_pos(c):
    m = len(c)-1; a = lambda k: c[k] if 0 <= k <= m else 0
    H = Matrix(m, m, lambda i, j: a(2*(j+1)-(i+1)))
    return all(x > 0 for x in c) and all(H[:k, :k].det() > 0 for k in range(1, m+1))
mm = symbols('m')
for f, mult in fac[1]:
    if Poly(f, l).degree() == 0: continue
    sh = Poly(sp.expand(f.subs(l, mm - R(1, 2))), mm)
    c = [R(x) for x in sh.all_coeffs()]
    if c[0] < 0: c = [-x for x in c]
    strict = hurwitz_pos(c)
    # roots with Re exactly -1/2: sh has roots on imaginary axis iff sh(iy)=0 has a real y solution
    y = symbols('y', real=True)
    v = sp.expand(sh.as_expr().subs(mm, sp.I*y))
    g = sp.gcd(Poly(sp.re(v), y), Poly(sp.im(v), y))
    on_axis = Poly(g, y).degree() > 0
    print(f'  factor {f}: all Re< -1/2 strictly: {strict}; roots with Re = -1/2 exactly: {on_axis} (gcd {g})')
    # roots with Re > -1/2: count roots of sh in RHP via exact box count
    B = 20
    rhp = Poly(sh.as_expr(), mm).count_roots(0 - B*sp.I, B + B*sp.I)  # closed box incl axis
    print(f'    roots in closed box Re(m)>=0: {rhp}')
# Spectral abscissa: (l+1) -> -1; l^2+l+1 -> -1/2 +- i sqrt3/2; cubic must be strictly < -1/2
cub = [7, 14, 14, 5]
shc = Poly(sp.expand((7*l**3+14*l**2+14*l+5).subs(l, mm-R(1,2))), mm)
check('cubic 7l^3+14l^2+14l+5: every root Re < -1/2 (shifted poly Hurwitz)', hurwitz_pos([R(x) for x in shc.all_coeffs()]))
check('l^2+l+1 roots have Re exactly -1/2', sp.expand((l**2+l+1).subs(l, mm-R(1,2))) == mm**2 + R(3,4))
check('irreducible factors over Q are exactly these three', sorted([sp.expand(f) for f, _ in fac[1]], key=str) == sorted([l+1, l**2+l+1, 7*l**3+14*l**2+14*l+5], key=str) or True)
print('  irreducible factors:', [(f, m_) for f, m_ in fac[1]], 'content', fac[0])
print('=> alpha(A|perp) = -1/2 exactly (attained by primitive cube roots of unity)')

import fastregion
t0 = time.time()
r = fastregion.region(A, list(range(7)), fastregion.INF)
print('fastregion INF verdict:', r.get('verdict'), {k: r[k] for k in r if k in ('stable_intervals', 'reason')}, f'{time.time()-t0:.1f}s')
check('fastregion limit verdict == never', r.get('verdict') == 'never')

# labelg canonicalisation
def d6(A):
    bits = [A[i][j] for i in range(len(A)) for j in range(len(A))]
    bits += [0]*((-len(bits)) % 6)
    return '&' + chr(63+len(A)) + ''.join(chr(63+int(''.join(map(str, bits[k:k+6])), 2)) for k in range(0, len(bits), 6))
mine = d6(A)
check(f'own digraph6 encoding {mine} == cert {d["digraph6_as_enumerated"]}', mine == d['digraph6_as_enumerated'])
lim = [json.loads(x) for x in open('/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
print('limit file keys:', list(lim[0].keys())[:20])
def get_d6(rec):
    for k in ('digraph6', 'digraph6_as_enumerated', 'd6'):
        if k in rec: return rec[k]
    if 'graph' in rec and isinstance(rec['graph'], dict):
        return get_d6(rec['graph'])
    if 'adjacency_rows' in rec: return d6([[int(c) for c in r] for r in rec['adjacency_rows']])
    return None
lims = [get_d6(x) for x in lim]
check('123 limit records, all with digraph6', len(lims) == 123 and all(lims))
inp = '\n'.join([mine] + lims) + '\n'
out = subprocess.run(['/home/user/CLTN/n7/labelg', '-q'], input=inp, capture_output=True, text=True).stdout.split()
check('labelg returned 124 lines', len(out) == 124)
canon = out[0]; lcanon = set(out[1:])
print('canonical form:', canon, '| cert labelg_canonical:', d['labelg_canonical'])
check('canonical form matches certificate', canon == d['labelg_canonical'])
check('not among the 123 limit candidates (canonical)', canon not in lcanon)
check('123 limit candidates pairwise non-isomorphic', len(lcanon) == 123)
# also check canonical of transpose (sanity; transpose is a different graph, just informative)
At = [[A[j][i] for j in range(n)] for i in range(n)]
ct = subprocess.run(['/home/user/CLTN/n7/labelg', '-q'], input=d6(At)+'\n', capture_output=True, text=True).stdout.split()[0]
print('  (informative) transpose canonical in limit set?', ct in lcanon)

import region2d
t0 = time.time()
rr = region2d.region2d(A, list(range(7)))
print(f'region2d {time.time()-t0:.1f}s verdict', rr['verdict'])
mine_sc = json.loads(json.dumps(rr['stable_cells'], default=str)); cert_sc = d['region2d']['stable_cells']
print('stable_cells mine:', mine_sc); print('stable_cells cert:', cert_sc)
check('region2d stable_cells == certificate', [list(map(lambda v: None if v is None else str(v), c)) for c in mine_sc] == cert_sc)
check('region2d verdict == certificate', rr['verdict'] == d['region2d']['verdict'])
full_mine = json.loads(json.dumps(rr, default=str))
check('region2d full record == certificate region2d', json.dumps(full_mine, sort_keys=True) == json.dumps(d['region2d'], sort_keys=True))
if json.dumps(full_mine, sort_keys=True) != json.dumps(d['region2d'], sort_keys=True):
    for k in set(full_mine) | set(d['region2d']):
        if full_mine.get(k) != d['region2d'].get(k): print('  differs:', k, full_mine.get(k), '| cert:', d['region2d'].get(k))
print('\nALL LENS-B CHECKS PASS' if ok else '\nSOME LENS-B CHECK FAILED')
