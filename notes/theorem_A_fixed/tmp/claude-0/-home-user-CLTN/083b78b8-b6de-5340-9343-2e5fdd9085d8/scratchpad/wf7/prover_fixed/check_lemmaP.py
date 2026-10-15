"""
check_lemmaP.py -- EXACT checks of the proved identities against explicit graphs (tests of the proofs, not the proofs):
 (C1) Lemma P: det(zI_n - aJ_n - D') = z^(p-1) Psi(z), checked at n + 2 exact rational z (which proves the polynomial identity),
      for R, Q, F (k = 1..3), Q2 (k = 1, 2), Y (k = 1) and 12 random
      (core, in-sets) pairs with 1-2 classes, at a in {3/7, 5, 200/21}.
 (C2) Lemma S: for random graphs and random rational (beta, t), the Hurwitz verdict of K = I + beta D^T + tJ computed through
      det(zI - D - aJ) shifted by -1/beta equals ctln.analyze_support's verdict.
 (C3) Lemma X2': for random (core, in-sets, sizes, beta), positivity of the explicit graph (ctln) equals
      [w_G(beta) > 0 and sigma_j(beta) > 1 for all j].
Writes check_lemmaP.json.
"""
import json, random
import sympy as sp
from fractions import Fraction as F
import fixlib as fl
from theorem_B import FAMILIES

random.seed(20261015)
out = {'C1': [], 'C2': [], 'C3': []}
avals = (sp.Rational(3, 7), sp.Integer(5), sp.Rational(200, 21))
G1 = [[0, 1, 1, 0, 1, 1], [1, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 1], [1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0]]
cases = []
for nm, ks in (('R', [1, 2, 3]), ('Q', [1, 2, 3]), ('Q2', [1, 2]), ('Y', [1])):
    base, classes, rhos = FAMILIES[nm]
    for kk in ks: cases.append((nm, fl.from_d6(base), classes, [r * kk for r in rhos]))
for kk in (1, 2, 3): cases.append(('F', G1, [(1, 1, 1, 1, 1, 0)], [kk]))
for i in range(12):
    m = random.randint(3, 6)
    A = [[0 if a == b else random.randint(0, 1) for b in range(m)] for a in range(m)]
    s = random.randint(1, 2)
    classes = [tuple(random.randint(0, 1) for _ in range(m)) for _ in range(s)]
    sizes = [random.randint(1, 3) for _ in range(s)]
    cases.append(('random%d' % i, A, classes, sizes))
for nm, A, classes, sizes in cases:
    ok = fl.check_identity_fast(A, classes, sizes)
    out['C1'].append({'case': nm, 'core_d6': fl.to_d6(A), 'classes': [list(S) for S in classes], 'sizes': sizes, 'ok': ok})
    print('C1', nm, sizes, ok, flush=True)
    assert ok
# (C2)
z = fl.z
for i in range(40):
    n = random.randint(3, 7)
    A = [[0 if a == b else random.choice([0, 1, 1]) for b in range(n)] for a in range(n)]
    beta = F(random.randint(201, 400), 100); t = F(random.randint(1, 400), random.randint(1, 20))
    eps, delta = fl.eps_delta(beta, t)
    r = fl.ctln.analyze_support(A, list(range(n)), eps, delta)
    D = sp.Matrix(fl.nonedge(A)); av = sp.Rational(str(t / beta)); rv = sp.Rational(str(1 / beta))
    chi = sp.expand((z * sp.eye(n) - D - av * sp.ones(n, n)).det(method='berkowitz'))
    h = sp.Poly(sp.expand((-1) ** n * chi.subs(z, -rv - fl.x)), fl.x)
    v = fl.ctln.hurwitz_stable([F(str(sp.Rational(c))) for c in h.all_coeffs()])
    out['C2'].append({'d6': fl.to_d6(A), 'beta': str(beta), 't': str(t), 'ctln': r['stable'], 'lemmaS': v[0]})
    assert v[0] == r['stable'], (fl.to_d6(A), beta, t, v, r['stable'])
print('C2 ok', sum(1 for c in out['C2'] if c['ctln'] is True), 'stable of', len(out['C2']))
# (C3)
for i in range(60):
    m = random.randint(3, 6)
    A = [[0 if a == b else random.choice([0, 1, 1, 1]) for b in range(m)] for a in range(m)]
    s = random.randint(1, 2)
    classes = [tuple(random.choice([0, 1, 1]) for _ in range(m)) for _ in range(s)]
    sizes = [random.randint(1, 3) for _ in range(s)]
    beta = F(random.randint(201, 500), 100)
    pos, w, sig = fl.positivity_core(A, classes, beta)
    G = fl.build(A, classes, sizes); n = len(G)
    eps, delta = fl.eps_delta(beta, F(7))
    r = fl.ctln.analyze_support(G, list(range(n)), eps, delta)
    out['C3'].append({'core': fl.to_d6(A), 'classes': [list(S) for S in classes], 'sizes': sizes, 'beta': str(beta),
                      'pred': pos, 'ctln_positive': r['positive']})
    assert pos == r['positive'], (fl.to_d6(A), classes, sizes, beta)
print('C3 ok', sum(1 for c in out['C3'] if c['pred']), 'positive of', len(out['C3']))
json.dump(out, open('check_lemmaP.json', 'w'), indent=1)
print('ALL CHECKS PASSED')
