"""Independent exact verifier (sympy Rationals only, no floats).
Usage: python3 -I verify.py <repo_root> <out_dir>
"""
import sys, json
from sympy import Rational, Matrix, eye, zeros, symbols, Poly, S

ROOT, OUT = sys.argv[1], sys.argv[2]
lam = symbols('lam')

def decode_d6(s):
    assert s[0] == '&', s
    n = ord(s[1]) - 63
    assert 0 <= n <= 62
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        assert 0 <= v < 64, (s, ch)
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    need = n * n
    assert len(bits) >= need and len(bits) - need < 6, (s, len(bits), need)
    assert all(b == 0 for b in bits[need:]), ('nonzero padding', s)
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]

def build_W(A, eps, delta):
    n = len(A)
    W = zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i == j:
                W[i, j] = S(0)
            elif A[j][i] == 1:          # j sends to i: weak
                W[i, j] = -1 + eps
            else:
                W[i, j] = -1 - delta
    return W

def hurwitz_minors(coeffs):
    """coeffs = [a0, a1, ..., an] of monic p(lam)=a0 lam^n + ... + an."""
    n = len(coeffs) - 1
    a = lambda k: coeffs[k] if 0 <= k <= n else S(0)
    H = Matrix(n, n, lambda i, j: a(2 * (j + 1) - (i + 1)))
    return [H[:k, :k].det(method='bareiss') for k in range(1, n + 1)]

def routh_first_column(coeffs):
    """Independent cross-check: Routh array first column (exact).
    Returns list or None if a zero pivot is hit (singular case)."""
    n = len(coeffs) - 1
    r0 = [coeffs[k] for k in range(0, n + 1, 2)]
    r1 = [coeffs[k] for k in range(1, n + 1, 2)]
    width = len(r0)
    r1 = r1 + [S(0)] * (width - len(r1))
    rows = [r0, r1]
    for _ in range(n - 1):
        p, q = rows[-2], rows[-1]
        if q[0] == 0:
            return None
        new = [(q[0] * p[k + 1] - p[0] * q[k + 1]) / q[0] if k + 1 < width else S(0) for k in range(width)]
        rows.append(new)
    return [r[0] for r in rows]

def det_eval(M, x0):
    n = M.shape[0]
    return (x0 * eye(n) - M).det(method='bareiss')

def check(A, eps, delta):
    n = len(A)
    eps, delta = Rational(eps), Rational(delta)
    res = {'eps': str(eps), 'delta': str(delta)}
    q = delta * (1 - eps) / eps
    res['q'] = str(q)
    res['legal'] = bool(eps > 0 and eps < 1 and delta > 0 and q > 1)
    W = build_W(A, eps, delta)
    M = eye(n) - W
    # closed-form cross-check: I - W = -delta I + (1+delta) J - (eps+delta) A^T ... per CLAUDE.md
    # CLAUDE.md writes A_sigma; entry (i,j) of I-W reads A[j][i], so the matrix form uses A^T.
    At = Matrix(A).T
    J = Matrix(n, n, lambda i, j: 1)
    res['closed_form_ok'] = bool((-delta * eye(n) + (1 + delta) * J - (eps + delta) * At) - M == zeros(n, n))
    D = M.det(method='bareiss')
    res['det_I_minus_W'] = str(D)
    if D == 0:
        res['positive'] = None
        res['x'] = None
    else:
        x = M.LUsolve(Matrix([1] * n))
        assert M * x == Matrix([1] * n)
        res['x'] = [str(v) for v in x]
        res['positive'] = bool(all(v > 0 for v in x))
    Jac = -eye(n) + W
    cp = Jac.charpoly(lam)               # Berkowitz (sympy default)
    coeffs = [S(c) for c in cp.all_coeffs()]
    assert coeffs[0] == 1
    # cross-check char poly by evaluation at n+2 integer points
    res['charpoly_crosscheck'] = all(Poly(cp.as_expr(), lam).eval(k) == det_eval(Jac, k) for k in range(-2, n + 1))
    res['charpoly_hi_to_lo'] = [str(c) for c in coeffs]
    mins = hurwitz_minors(coeffs)
    res['hurwitz_minors_sign'] = [(1 if bool(m > 0) else (-1 if bool(m < 0) else 0)) for m in mins]
    if any(m == 0 for m in mins):
        res['stable'] = 'UNDECIDED'
    elif all(m > 0 for m in mins):
        res['stable'] = True
    else:
        res['stable'] = False
    rf = routh_first_column(coeffs)
    if rf is None:
        res['routh'] = 'zero_pivot'
    else:
        res['routh'] = bool(all(v > 0 for v in rf))
    res['routh_agrees'] = (res['routh'] == res['stable']) if res['stable'] != 'UNDECIDED' else None
    return res

def is_complete_bidirected(A):
    n = len(A)
    return all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)

def main():
    # decoder self-test
    G = decode_d6('&EZjtwWE')
    arcs = sum(map(sum, G))
    gen_file = open(f'{ROOT}/n7/geneson.d6').read().split()
    selftest = {'geneson_rows': [''.join(map(str, r)) for r in G], 'arcs': arcs,
                'node0_out': [j for j in range(6) if G[0][j]],
                'geneson_d6_file': gen_file}
    F = json.load(open(f'{ROOT}/certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))
    fn_rows = F['graph']['adjacency_rows']
    selftest['FNFL_rows_match'] = [''.join(map(str, r)) for r in decode_d6('&FNFL~[udo?')] == fn_rows
    # positive control: Geneson at eps=1/100, delta=29/2500 (test_ctln.py claims stable FP, not clique)
    selftest['geneson_control'] = {k: v for k, v in check(G, '1/100', '29/2500').items() if k in ('legal','positive','stable','routh')}
    selftest['geneson_control_eps_1_10'] = {k: v for k, v in check(G, '1/10', '29/250').items() if k in ('legal','positive','stable','routh')}
    K7 = [[0 if i == j else 1 for j in range(7)] for i in range(7)]
    selftest['K7_control'] = {k: v for k, v in check(K7, '1/49', '33/1568').items() if k in ('legal','positive','stable','routh')}
    print('SELFTEST', selftest, flush=True)

    fin = [json.loads(l) for l in open(f'{ROOT}/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
    lim = {r['d6']: r for r in (json.loads(l) for l in open(f'{ROOT}/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl'))}
    rows, failures = [], []
    for k, r in enumerate(fin):
        d6 = r['d6']
        A = decode_d6(d6)
        n = len(A)
        meta = {'d6': d6, 'n': n, 'arcs': sum(map(sum, A)), 'loops': sum(A[i][i] for i in range(n))}
        meta['arcs_match_record'] = meta['arcs'] == r['edges']
        cb = is_complete_bidirected(A)
        meta['clique'] = cb
        meta['clique_match_record'] = (cb == r['clique'])
        if d6 in lim:
            L = lim[d6]
            meta['A_match_limit_record'] = (L['A'] == A)
            w = L['pointwise_witnesses'][0]
            pa = (w['eps'], w['delta'], 'limit file pointwise_witnesses[0]')
        else:
            meta['A_match_limit_record'] = None
            pa = ('1/49', '33/1568', 'FNFL json witness 1')
        c = r['stable_cells'][0]
        beta = Rational(c[2]); tlo = Rational(c[3])
        t = tlo + 1 if c[4] is None else (tlo + Rational(c[4])) / 2
        epsb = 1 / (t + 1); deltab = (beta - 1) * epsb
        pts = [pa, (str(epsb), str(deltab), f'first stable cell: beta={beta}, t={t}')]
        if d6 not in lim:
            pts.append(('1/36', '17/588', 'FNFL json witness 2'))
        for (e, dl, src) in pts:
            res = check(A, e, dl)
            row = dict(meta); row['point_source'] = src; row.update(res)
            ok = (res['legal'] is True and res['positive'] is True and res['stable'] is True
                  and not cb and res['closed_form_ok'] and res['charpoly_crosscheck']
                  and res['routh_agrees'] is True and meta['arcs_match_record']
                  and meta['clique_match_record'] and meta['loops'] == 0
                  and meta['A_match_limit_record'] in (True, None))
            row['verdict'] = 'COUNTEREXAMPLE_CONFIRMED' if ok else 'CHECK_FAILED'
            if not ok:
                failures.append({k2: row[k2] for k2 in ('d6', 'point_source', 'eps', 'delta', 'legal', 'positive', 'stable', 'clique')})
            rows.append(row)
        print(k, d6, [x['verdict'] for x in rows if x['d6'] == d6], flush=True)
    graphs_ok = sorted({x['d6'] for x in rows} - {f['d6'] for f in failures})
    summary = {'n_graphs': len(fin), 'n_points': len(rows), 'n_graphs_all_points_confirmed': len(graphs_ok),
               'n_failed_points': len(failures), 'failures': failures, 'selftest': selftest,
               'n_undecided': sum(1 for x in rows if x['stable'] == 'UNDECIDED')}
    json.dump({'summary': summary, 'rows': rows}, open(f'{OUT}/v_sympy_results.json', 'w'), indent=1)
    print('SUMMARY', json.dumps({k: v for k, v in summary.items() if k != 'selftest'}))

main()
