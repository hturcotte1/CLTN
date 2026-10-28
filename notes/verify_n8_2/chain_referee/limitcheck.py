"""Referee's own exact limit-stage checker (sympy, written from the definitions).
Limit ray eps = 1/(t+1), delta = (beta-1) eps.  Then I - W = eps K with K = C + t J, C = I + beta N,
N[i][j] = 1 iff i != j and A[j][i] == 0 (j does not send to i).  Rank-one facts (re-derived, not copied):
  adj(K) 1 = adj(C) 1 =: a(beta)   (t-free);  det K = det C + t * 1^T adj(C) 1 = d + t s;
  det(mu I + K) = det(mu I + C) + t * (det(mu I + C + J) - det(mu I + C)) = phi + t * omega.
Limit (t -> inf): fixed point x = a / (d + t s) -> positivity iff all a_i have the sign of s (Lemma A.2: stable => det K > 0 => s > 0 for large t => a_i > 0);
stability for large t iff omega(mu; beta) (degree m-1, leading coefficient m) is Hurwitz in mu (Lemma I.3: C|perp has Re spec > 0).
Critical set in beta > 2: real roots of a_i, s (constant term of omega) and the Hurwitz minor Delta_{m-2}(omega) (Orlando).
On every open cell between consecutive critical values the verdict is constant (root continuity; leading coefficient m fixed)."""
import json, sys, time
from fractions import Fraction as Fr
import sympy as sp
from sympy import Matrix, symbols, Poly, eye, ones, Rational, QQ
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/chain_referee')
from myfilter2 import decode_d6, hurwitz_minors
beta, mu = symbols('beta mu')
def hurwitz3(coeffs):
    """3-valued Hurwitz test on exact coefficients (highest first, leading > 0): True / False / None (a vanishing minor or coefficient)."""
    c = [Fr(int(x.p), int(x.q)) if hasattr(x, 'p') else Fr(x) for x in coeffs]
    assert c[0] > 0
    if any(x < 0 for x in c): return False
    if any(x == 0 for x in c): return None
    m = hurwitz_minors(c)
    if any(x < 0 for x in m): return False
    if any(x == 0 for x in m): return None
    return True
def analyse(d6):
    n, A = decode_d6(d6)
    N = Matrix(n, n, lambda i, j: 0 if i == j else 1 - A[j][i])
    C = eye(n) + beta * N
    a = (C.adjugate() * ones(n, 1))
    a = [sp.expand(x) for x in a]
    s = sp.expand(sum(a))
    phi = (mu * eye(n) + C).det(method='berkowitz')
    phiJ = (mu * eye(n) + C + ones(n, n)).det(method='berkowitz')
    omega = Poly(sp.expand(phiJ - phi), mu)
    oc = omega.all_coeffs()                         # highest first, polynomials in beta
    assert omega.degree() == n - 1 and sp.expand(oc[0]) == n and sp.expand(oc[-1] - s) == 0
    d = n - 1
    H = Matrix(d - 1, d - 1, lambda i, j: oc[2 * i - j + 1] if 0 <= 2 * i - j + 1 <= d else 0)
    Delta = sp.expand(H.det(method='berkowitz'))
    polys = {f'a{i}': a[i] for i in range(n)}; polys['s'] = s; polys['Delta'] = Delta
    return n, A, a, s, oc, Delta, polys
def roots_gt2(polys):
    """distinct real roots > 2 of the product (squarefree part), as exact isolating intervals (sympy)."""
    prod = Poly(1, beta, domain=QQ)
    for p in polys.values():
        P = Poly(p, beta, domain=QQ)
        assert not P.is_zero
        prod = prod * P
    sqf = Poly(sp.sqf_part(prod.as_expr()), beta, domain=QQ)
    ivs = sqf.intervals(inf=Rational(2), sqf=True)
    ivs = [tuple(iv) for iv in ivs]
    at2 = sqf.count_roots(Rational(2), Rational(2))
    assert sum(1 for iv in ivs if iv[0] == iv[1] == 2) == at2, 'root at beta=2 not isolated exactly'
    ivs = [iv for iv in ivs if iv[1] > 2]       # strict: roots exactly at 2 excluded
    return sqf, ivs
def verdict_at(n, a, s, oc, b0):
    b0 = Rational(b0)
    sv = s.subs(beta, b0)
    if sv == 0: return None, None
    pos = all(x.subs(beta, b0) / sv > 0 for x in a)
    st = hurwitz3([c.subs(beta, b0) for c in oc])
    return pos, st
def check_record(rec):
    t0 = time.time(); out = {'d6': rec['d6'], 'verdict_recorded': rec['verdict']}
    n, A, a, s, oc, Delta, polys = analyse(rec['d6'])
    cert = rec['cert']; crit = cert['critical']; cells = cert['cells']
    sqf, ivs = roots_gt2(polys)
    out['my_critical_count'] = len(ivs); out['cert_critical_count'] = len(crit)
    # each cert critical interval contains exactly one of my roots (closed interval root count)
    per = []
    for c in crit:
        lo, hi = Rational(c['interval'][0]), Rational(c['interval'][1])
        per.append(sqf.count_roots(lo, hi))
    out['roots_per_cert_interval'] = per
    # cells chain from beta_min=2 to infinity through consecutive critical indices
    order = sorted(range(len(crit)), key=lambda i: Rational(crit[i]['interval'][0]))
    chain = [(-1 if cells[0]['lo'] == -1 else None)]
    chain_ok = cells[0]['lo'] == -1 and cells[-1]['hi'] is None and all(cells[i]['hi'] == cells[i + 1]['lo'] for i in range(len(cells) - 1))
    chain_ok = chain_ok and [c['hi'] for c in cells[:-1]] == order
    out['cells_chain_ok'] = chain_ok
    mine_cells = []; any_stable = False; mism = 0
    for c in cells:
        b0 = Rational(c['sample_beta'])
        # sample strictly inside its cell
        lo_ok = b0 > 2 if c['lo'] == -1 else b0 > Rational(crit[c['lo']]['interval'][1])
        hi_ok = True if c['hi'] is None else b0 < Rational(crit[c['hi']]['interval'][0])
        pos, st = verdict_at(n, a, s, oc, b0)
        mine_cells.append({'sample_beta': str(b0), 'inside': bool(lo_ok and hi_ok), 'positive': pos, 'stable': st, 'cert_positive': c['positive'], 'cert_stable': c['stable']})
        if pos != c['positive'] or st != c['stable'] or not (lo_ok and hi_ok): mism += 1
        if pos and st: any_stable = True
    out['cells'] = mine_cells; out['cell_mismatches'] = mism
    out['my_verdict'] = 'stable_somewhere' if any_stable else ('never' if all(c['positive'] is not None and c['stable'] is not None for c in mine_cells) else 'escalate')
    out['agree'] = (out['my_verdict'] == rec['verdict'] and mism == 0 and chain_ok and len(ivs) == len(crit) and all(p == 1 for p in per))
    out['seconds'] = round(time.time() - t0, 1)
    return out
if __name__ == '__main__':
    which = sys.argv[1]; k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    recs = []
    for l in open('n8/n8_limit_results.jsonl'):
        r = json.loads(l)
        if which == 'never_full' and r['verdict'] == 'never' and 'cert' in r: recs.append(r)
        elif which == 'stable' and r['verdict'] == 'stable_somewhere' and not r['complete']: recs.append(r)
        elif which == 'complete' and r['complete']: recs.append(r)
        if len(recs) >= k: break
    res = []
    for r in recs:
        o = check_record(r); res.append(o)
        print(json.dumps({k2: v for k2, v in o.items() if k2 != 'cells'}), flush=True)
        print('   cells:', [(c['sample_beta'], c['positive'], c['stable'], c['inside']) for c in o['cells']], flush=True)
    json.dump(res, open(f'/home/user/CLTN/notes/verify_n8_2/chain_referee/limitcheck_{which}.json', 'w'), indent=1, default=str)
    print('ALL AGREE' if all(o['agree'] for o in res) else 'DISAGREEMENT', len(res))
