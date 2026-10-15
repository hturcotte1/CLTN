import json, sys, time
from fractions import Fraction as F
sys.path.insert(0, '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf6/v_structure')
from exact import *

R = '/home/user/CLTN/'
OUT = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf6/v_structure/'
problems = []
def flag(msg):
    problems.append(msg); print('DISCREPANCY:', msg)

# ---- self tests of primitives ----
assert hurwitz([1, 3, 3, 1])[0] == 'stable'
assert hurwitz([1, 0, 1])[0] == 'undecided'
assert hurwitz([1, -1, 1])[0] == 'unstable'
assert hurwitz([1, 1, 1, 6])[0] == 'unstable'      # (x+2)(x^2 - x + 3)
assert charpoly([[0, 1], [1, 0]]) == [1, 0, -1]
assert poly_shift([1, 0, -1], F(1)) == [1, 2, 0]  # (x+1)^2-1
for s in ['&EZjtwWE', '&FNFL~[udo?']:
    assert d6_encode(d6_decode(s)) == s

# ---- digraph6 decoder check on Geneson ----
G = d6_decode('&EZjtwWE')
gen_file = open(R + 'n7/geneson.d6').read().split()[0]
geneson = {'d6': '&EZjtwWE', 'n': len(G), 'arcs': n_arcs(G), 'rows': [''.join(map(str, r)) for r in G],
           'node0_out': [j for j in range(6) if G[0][j]], 'matches_n7_geneson_d6': gen_file == '&EZjtwWE',
           'loops': [i for i in range(6) if G[i][i]]}
print('Geneson', geneson)
# independent conversion from the paper's nonedge matrix N (N[a][b]=1 iff b -/-> a), as quoted verbatim in notes
Npaper = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]]
Agen = [[0 if a == b else 1 - Npaper[b][a] for b in range(6)] for a in range(6)]   # A[a][b]=1 iff a->b iff N[b][a]=0
geneson['equals_paper_N_converted'] = (G == Agen)
geneson['isomorphic_to_paper'] = canon(G) == canon(Agen)
if not (len(G) == 6 and n_arcs(G) == 19 and G == Agen and not geneson['loops'] and geneson['matches_n7_geneson_d6']):
    flag('Geneson decode mismatch')

fin = [json.loads(l) for l in open(R + 'certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl')]
lim = [json.loads(l) for l in open(R + 'certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
fnfl = json.load(open(R + 'certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json'))
print(len(fin), len(lim))
if len(fin) != 124: flag('finite count %d' % len(fin))
if len(lim) != 123: flag('limit count %d' % len(lim))

# FNFL adjacency rows vs decode
Af = d6_decode('&FNFL~[udo?')
if [''.join(map(str, r)) for r in Af] != fnfl['graph']['adjacency_rows']: flag('FNFL rows mismatch')

# ---- (i) structure ----
t0 = time.time()
rows = []
canon_fin = {}
for rec in fin:
    A = d6_decode(rec['d6'])
    n = len(A)
    row = {'d6': rec['d6'], 'n': n}
    if n != 7: flag('%s n=%d' % (rec['d6'], n))
    loops = [i for i in range(n) if A[i][i]]
    if loops: flag('%s has loops %s' % (rec['d6'], loops))
    row['arcs'] = n_arcs(A); row['missing_arcs'] = missing_arcs(A)
    row['complete_bidirected'] = is_complete_bidirected(A)
    row['clique'] = row['complete_bidirected']
    if row['arcs'] != rec['edges']: flag('%s edges field %s != %d' % (rec['d6'], rec['edges'], row['arcs']))
    if rec['clique'] != row['clique']: flag('%s clique field' % rec['d6'])
    if row['complete_bidirected']: flag('%s IS complete bidirected' % rec['d6'])
    dp = dominated_pairs(A)
    row['dominated_pairs_k_dominates_j'] = dp
    row['domination_free'] = not dp
    if rec['domination_free'] != row['domination_free']: flag('%s domination_free field %s but pairs %s' % (rec['d6'], rec['domination_free'], dp))
    c = canon(A)
    row['canon_d6'] = d6_encode([list(c[i * n:(i + 1) * n]) for i in range(n)])
    if c in canon_fin: flag('%s isomorphic to %s' % (rec['d6'], canon_fin[c]))
    canon_fin[c] = rec['d6']
    # the record's own 'canonical' field should be isomorphic to d6
    if 'canonical' in rec and rec['canonical']:
        if canon(d6_decode(rec['canonical'])) != c: flag('%s canonical field %s not isomorphic' % (rec['d6'], rec['canonical']))
    row['_A'] = A
    rows.append(row)
print('structure done %.1fs' % (time.time() - t0))

canon_lim = {}
for L in lim:
    A = d6_decode(L['d6'])
    if A != L['A']: flag('limit %s: A field != decode' % L['d6'])
    if n_arcs(A) != L['edges']: flag('limit %s edges' % L['d6'])
    if L.get('complete') or L.get('clique'): flag('limit %s flagged complete/clique' % L['d6'])
    c = canon(A)
    if c in canon_lim: flag('limit file duplicate %s ~ %s' % (L['d6'], canon_lim[c]))
    canon_lim[c] = L['d6']
fin_in_lim = {c: d for c, d in canon_fin.items()}
for row, rec in zip(rows, fin):
    c = canon(row['_A']) if False else None
inlim_set = set()
for rec in fin:
    c = canon(d6_decode(rec['d6']))
    member = c in canon_lim
    if member != rec['in_limit_candidates']:
        flag('%s in_limit_candidates=%s but membership=%s' % (rec['d6'], rec['in_limit_candidates'], member))
    if member: inlim_set.add(c)
if set(canon_lim) != inlim_set:
    flag('limit graphs not covered by finite file: %s' % [canon_lim[c] for c in set(canon_lim) - inlim_set])
fnfl_in_lim = canon(Af) in canon_lim
if fnfl_in_lim: flag('FNFL is among limit candidates')
lim_by_canon = {canon(d6_decode(L['d6'])): L for L in lim}
print('isomorphism checks done %.1fs' % (time.time() - t0))

# ---- exact CTLN at a finite point ----
def W_of(A, eps, delta):
    n = len(A)
    return [[F(0) if i == j else ((-1 + eps) if A[j][i] == 1 else (-1 - delta)) for j in range(n)] for i in range(n)]

def finite_check(A, eps, delta):
    eps = F(eps); delta = F(delta); n = len(A)
    q = delta * (1 - eps) / eps
    legal = (0 < eps < 1) and delta > 0 and q > 1
    W = W_of(A, eps, delta)
    # cross-check closed form  I - W = -delta I + (1+delta) J - (eps+delta) A^T
    IW = [[F(int(i == j)) - W[i][j] for j in range(n)] for i in range(n)]
    CF = [[-delta * int(i == j) + (1 + delta) - (eps + delta) * A[j][i] for j in range(n)] for i in range(n)]
    assert IW == CF
    x = solve(IW, [1] * n)
    positive = x is not None and all(v > 0 for v in x)
    Jac = [[W[i][j] - int(i == j) for j in range(n)] for i in range(n)]
    st, minors = hurwitz(charpoly(Jac))
    return {'eps': str(eps), 'delta': str(delta), 'q': str(q), 'legal': legal,
            'singular': x is None, 'positive': positive, 'stable': st,
            'min_x': str(min(x)) if x else None}

# ---- limit check ----
def limit_check(A, beta):
    beta = F(beta); n = len(A)
    gamma = 1 - 1 / beta
    J = [[1] * n for _ in range(n)]
    AJ = [[A[i][j] + 1 for j in range(n)] for i in range(n)]
    cA = charpoly(A); cAJ = charpoly(AJ)
    mchi = poly_sub(cA, cAJ)
    assert len(mchi) == n and mchi[0] == n, mchi   # degree n-1, leading coefficient n
    chi = [c / n for c in mchi]
    assert all(c.denominator == 1 for c in mchi), mchi
    g = poly_shift(chi, -gamma)          # g(x) = chi(x - gamma)
    st, minors = hurwitz(g)
    AT = [[A[j][i] + (gamma if i == j else 0) for j in range(n)] for i in range(n)]
    w = solve(AT, [1] * n)
    wpos = w is not None and all(v > 0 for v in w)
    return {'beta': str(beta), 'gamma': str(gamma), 'legal_limit': beta > 2,
            'm_chi_perp': [int(c) for c in mchi], 'shifted_hurwitz': st,
            'w_positive': wpos, 'w_singular': w is None}, chi

# sanity: limit chi_perp via an independent route (basis e_i - e_n of 1^perp, G_V = (V^T V)^-1 V^T A V)
def chi_perp_basis(A):
    n = len(A)
    V = [[(1 if i == k else 0) - (1 if i == n - 1 else 0) for k in range(n - 1)] for i in range(n)]
    VT = [list(r) for r in zip(*V)]
    VTV = matmul(VT, V); VTAV = matmul(matmul(VT, A), V)
    # solve VTV X = VTAV column by column
    cols = [solve(VTV, [VTAV[i][k] for i in range(n - 1)]) for k in range(n - 1)]
    X = [[cols[k][i] for k in range(n - 1)] for i in range(n - 1)]
    return charpoly(X)

table = []
for row, rec in zip(rows, fin):
    A = row.pop('_A')
    d6 = rec['d6']
    out = dict(row)
    out['record_limit_verdict'] = rec['limit_verdict']
    out['in_limit_candidates_record'] = rec['in_limit_candidates']
    cells = rec['stable_cells']
    # every stable cell sample point at finite (beta, t), with t_sample from region2d cells
    tsamp = {}
    for c in rec['region2d']['cells']:
        for tc in c['t_cells']:
            tsamp[(c['beta_sample'], tc['t_lo'], tc['t_hi'])] = tc['t_sample']
    cell_checks = []
    for (blo, bhi, bs, tlo, thi) in cells:
        ts = tsamp.get((bs, tlo, thi))
        if ts is None: flag('%s: no t_sample for stable cell %s' % (d6, (blo, bhi, bs, tlo, thi))); continue
        beta = F(bs); t = F(ts)
        if not (F(blo) <= beta <= F(bhi)): flag('%s beta_sample outside cell' % d6)
        if not (F(tlo) <= t and (thi is None or t <= F(thi))): flag('%s t_sample outside cell' % d6)
        eps = 1 / (t + 1); delta = (beta - 1) * eps
        fc = finite_check(A, eps, delta)
        fc['beta'] = bs; fc['t'] = ts
        fc['legal_t_beta_form'] = t * (beta - 2) > 1
        assert fc['legal_t_beta_form'] == fc['legal']
        ok = fc['legal'] and fc['positive'] and fc['stable'] == 'stable'
        fc['counterexample_confirmed'] = ok
        if not ok: flag('%s stable cell sample (beta=%s,t=%s) fails: %s' % (d6, bs, ts, fc))
        cell_checks.append(fc)
    out['finite_cell_samples'] = cell_checks
    out['n_cell_samples_confirmed'] = sum(c['counterexample_confirmed'] for c in cell_checks)
    if rec['in_limit_candidates']:
        # first stable cell that is unbounded in t
        first = cells[0]
        if first[4] is not None:
            flag('%s first stable cell bounded in t' % d6)
        lc, chi = limit_check(A, first[2])
        chi2 = chi_perp_basis(A)
        if [F(c) for c in chi] != chi2: flag('%s chi_perp two routes disagree' % d6)
        lc['ok'] = lc['legal_limit'] and lc['shifted_hurwitz'] == 'stable' and lc['w_positive']
        if not lc['ok']: flag('%s limit check failed %s' % (d6, lc))
        out['limit_check'] = lc
        # pointwise witnesses from the limit file
        L = lim_by_canon[canon(A)]
        A_L = d6_decode(L['d6'])
        pw = []
        for wtn in L['pointwise_witnesses']:
            fc = finite_check(A_L, F(wtn['eps']), F(wtn['delta']))
            fc['record_says'] = {'legal': wtn['legal'], 'is_fixed_point': wtn['is_fixed_point'], 'stable': wtn['stable']}
            # compare x
            fc['counterexample_confirmed'] = fc['legal'] and fc['positive'] and fc['stable'] == 'stable'
            claimed = wtn['legal'] and wtn['is_fixed_point'] and wtn['stable'] and not wtn['agrees_with_conjecture']
            if claimed != fc['counterexample_confirmed']:
                flag('%s (limit d6 %s) witness eps=%s delta=%s: record claims %s, mine %s' % (d6, L['d6'], wtn['eps'], wtn['delta'], claimed, fc))
            pw.append(fc)
        out['limit_file_d6'] = L['d6']
        out['limit_file_witnesses'] = pw
        out['n_limit_witnesses_confirmed'] = sum(p['counterexample_confirmed'] for p in pw)
    table.append(out)
print('per-graph checks done %.1fs' % (time.time() - t0))

json.dump({'geneson': geneson, 'table': table, 'problems': problems,
           'fnfl_in_limit_candidates': fnfl_in_lim}, open(OUT + 'verify_out.json', 'w'), indent=1)
print('problems:', len(problems))
