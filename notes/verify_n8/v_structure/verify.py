"""v_structure verifier for the n = 8 limit candidates (independent implementation, exact arithmetic only).
Runs the full method on ALL records; 3 worker processes."""
import json, sys, time, math, subprocess, collections
from fractions import Fraction as F
from multiprocessing import Pool
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8/v_structure')
from exact8 import *

R = '/home/user/CLTN/'
OUT = '/home/user/CLTN/notes/verify_n8/v_structure/'
IN = R + 'n8/n8_limit_HITS.jsonl'

def simplest_between(L, U):
    """the rational with smallest denominator strictly inside (L, U), L < U (Stern-Brocot)."""
    L = F(L); U = F(U)
    assert L < U
    fl = math.floor(L)
    if fl + 1 < U and fl + 1 > L: return F(fl + 1)
    if L == fl:  # L integer: look in (L, U)
        # continue on reciprocal of fractional parts
        pass
    # recursive continued-fraction construction
    if math.floor(L) + 1 <= math.floor(U) and math.floor(L) + 1 < U:
        return F(math.floor(L) + 1) if math.floor(L) + 1 > L else F(math.floor(L) + 2)
    # same integer part
    i = math.floor(L)
    l2 = L - i; u2 = U - i
    if l2 == 0:
        # (0, u2): simplest is 1/ceil(1/u2 + tiny)
        k = math.floor(1 / u2) + 1
        return i + F(1, k)
    return i + 1 / simplest_between(1 / u2, 1 / l2)

def root_side(beta, a, b, p):
    """position of rational beta relative to the unique root r of p in [a, b]: -1 (beta < r), 0 (==), +1 (beta > r).
    a == b means the exact rational root a (p(a) == 0 is checked by check_isolating)."""
    beta = F(beta); a = F(a); b = F(b)
    if a == b: return -1 if beta < a else (0 if beta == a else 1)
    if beta <= a: return -1
    if beta >= b: return 1
    pa = poly_eval(p, a); pb = poly_eval(p, beta)
    if pb == 0: return 0
    return -1 if (pb > 0) == (pa > 0) else 1

def endpoint_lo_bound(lo):
    """(certified lower bound L for 'strictly above lo', root-side test data)."""
    if lo['kind'] == 'beta_min':
        return F(lo['beta']), None
    a, b = map(F, lo['beta_interval'])
    return b, (a, b, lo['beta_minpoly_hi_to_lo'])

def strictly_inside(beta, iv):
    lo, hi = iv['lo'], iv['hi']
    beta = F(beta)
    if lo['kind'] == 'beta_min':
        above = beta > F(lo['beta'])
    else:
        a, b = map(F, lo['beta_interval'])
        above = root_side(beta, a, b, lo['beta_minpoly_hi_to_lo']) > 0
    a, b = map(F, hi['beta_interval'])
    below = root_side(beta, a, b, hi['beta_minpoly_hi_to_lo']) < 0
    return above and below

def check_isolating(a, b, p):
    """exact: p(a) != 0, p(b) != 0, sign change, and Sturm count of distinct roots in (a, b] equals 1."""
    a = F(a); b = F(b)
    if a == b:   # exact rational root
        return poly_eval(p, a) == 0
    pa = poly_eval(p, a); pb = poly_eval(p, b)
    if a > b or pa == 0 or pb == 0 or (pa > 0) == (pb > 0): return False
    return sturm_count(sturm_sequence(p), a, b) == 1

def refine_root_interval(a, b, p, steps=1):
    """exact bisection of an isolating interval (a, b) of p (unique simple root, sign change): returns (a', b')."""
    a = F(a); b = F(b)
    for _ in range(steps):
        if a == b: break
        m = (a + b) / 2; pm = poly_eval(p, m)
        if pm == 0: a = b = m
        elif (pm > 0) == (poly_eval(p, a) > 0): a = m
        else: b = m
    return a, b

def work(line):
    rec = json.loads(line)
    d6 = rec['d6']; A = d6_decode(d6); n = len(A)
    probs = []
    row = {'d6': d6, 'n': n}
    if n != 8: probs.append('n=%d' % n)
    loops = [i for i in range(n) if A[i][i]]
    if loops: probs.append('loops %s' % loops)
    row['arcs'] = n_arcs(A); row['missing_arcs'] = missing_arcs(A)
    row['complete_bidirected'] = is_complete_bidirected(A)
    if row['arcs'] != rec['edges']: probs.append('edges field %s != arcs %d' % (rec['edges'], row['arcs']))
    if row['complete_bidirected']: probs.append('IS complete bidirected')
    if rec.get('complete') is not False: probs.append('complete field is %r' % rec.get('complete'))
    if rec.get('verdict') != 'stable_somewhere': probs.append('verdict %r' % rec.get('verdict'))
    cert = rec['cert']
    if cert.get('sigma') != list(range(8)): probs.append('sigma %r' % cert.get('sigma'))
    if cert.get('beta_min') != '2': probs.append('beta_min %r' % cert.get('beta_min'))
    if cert.get('verdict') != 'stable_somewhere': probs.append('cert verdict %r' % cert.get('verdict'))
    dp = dominated_pairs(A)
    row['dominated_pairs_k_dominates_j'] = dp
    row['domination_free'] = not dp
    if dp: probs.append('NOT domination free: %s' % dp)
    row['canon'] = ''.join(map(str, canon(A)))
    # ---- stable intervals
    ivs = cert['stable_intervals']
    row['n_stable_intervals'] = len(ivs)
    if len(ivs) != 1: probs.append('%d stable intervals' % len(ivs))
    ivinfo = []
    for iv in ivs:
        lo, hi = iv['lo'], iv['hi']
        info = {'lo_kind': lo['kind'], 'hi_kind': hi['kind']}
        if hi['kind'] != 'root': probs.append('hi kind %r' % hi['kind'])
        a, b = map(F, hi['beta_interval']); p = hi['beta_minpoly_hi_to_lo']
        info['hi_beta_interval'] = [str(a), str(b)]; info['hi_deg'] = len(p) - 1
        info['hi_isolating_ok'] = check_isolating(a, b, p)
        if not info['hi_isolating_ok']: probs.append('hi isolating interval not certified')
        qa, qb = map(F, hi['q_interval'])
        info['hi_q_interval'] = [str(qa), str(qb)]
        if (qa, qb) != (a - 1, b - 1): probs.append('hi q_interval != beta_interval - 1')
        if lo['kind'] == 'beta_min':
            if F(lo['beta']) != 2 or F(lo['q']) != 1: probs.append('lo beta_min fields %r' % lo)
            L = F(2); info['lo_beta'] = '2'
        elif lo['kind'] == 'root':
            c, d = map(F, lo['beta_interval']); pl = lo['beta_minpoly_hi_to_lo']
            info['lo_beta_interval'] = [str(c), str(d)]; info['lo_deg'] = len(pl) - 1
            info['lo_isolating_ok'] = check_isolating(c, d, pl)
            if not info['lo_isolating_ok']: probs.append('lo isolating interval not certified')
            if not (c >= 2): probs.append('lo root interval below beta_min')
            L = d
            if 'q_interval' in lo:
                if tuple(map(F, lo['q_interval'])) != (c - 1, d - 1): probs.append('lo q_interval != beta_interval - 1')
        else:
            probs.append('lo kind %r' % lo['kind']); L = None
        info['hi_q_upper_bound'] = str(b - 1)
        info['certified_gap'] = (L is not None and L < a)
        info['note'] = None if (L is not None and L < a) else 'isolating intervals touch/overlap (L=%s, a=%s): interior certified by sign test' % (L, a)
        info['_L'] = L; info['_a'] = a
        ivinfo.append(info)
    row['intervals'] = ivinfo
    # ---- cells vs intervals (light consistency check)
    cells = cert.get('cells', [])
    bad_cells = 0
    for c in cells:
        sb = F(c['sample_beta'])
        claimed = bool(c.get('stable')) and bool(c.get('positive')) and c.get('stable') is True
        inside = any(strictly_inside(sb, iv) for iv in ivs)
        if claimed != inside: bad_cells += 1
    row['n_cells'] = len(cells); row['cells_inconsistent_with_intervals'] = bad_cells
    if bad_cells: probs.append('%d cells inconsistent with stable_intervals' % bad_cells)
    # ---- witnesses
    wl = []
    first_illegal = None
    betas = set()
    for k, w in enumerate(rec['witnesses']):
        t = w['t']; beta = F(w['beta']); betas.add(beta)
        eps = F(1, t + 1); delta = (beta - 1) * eps; q = delta * (1 - eps) / eps
        wi = {'t': t, 'beta': str(beta)}
        cons = (F(w['eps']) == eps and F(w['delta']) == delta and F(w['q']) == q and isinstance(t, int))
        wi['fields_consistent'] = cons
        if not cons: probs.append('witness %d fields inconsistent' % k)
        wi['beta_strictly_inside'] = strictly_inside(beta, ivs[0])
        if not wi['beta_strictly_inside'] and w['stable']: probs.append('witness %d beta not strictly inside stable interval but stored stable=true' % k)
        legal = q > 1
        wi['legal'] = legal
        if k == 0: first_illegal = not legal
        r = ctln_check(A, eps, delta)
        wi.update({'eps': r['eps'], 'delta': r['delta'], 'q': r['q'], 'stable': r['stable'], 'positive': r['positive']})
        wi['stored_stable'] = w['stable']; wi['stored_positive'] = w['positive']
        if (r['stable'] == 'stable') != bool(w['stable']): probs.append('witness %d stable: stored %s mine %s' % (k, w['stable'], r['stable']))
        if r['positive'] != bool(w['positive']): probs.append('witness %d positive: stored %s mine %s' % (k, w['positive'], r['positive']))
        wi['counterexample_confirmed'] = legal and r['positive'] and r['stable'] == 'stable'
        wl.append(wi)
    row['witnesses'] = wl
    row['witness_betas_identical'] = len(betas) == 1
    row['n_witnesses_outside_interval'] = sum(not x['beta_strictly_inside'] for x in wl)
    row['first_witness_illegal'] = first_illegal
    row['n_witnesses_confirmed'] = sum(x['counterexample_confirmed'] for x in wl)
    leg = next((x for x in wl if x['counterexample_confirmed']), None)
    row['stored_legal_witness'] = {k: leg[k] for k in ('eps', 'delta', 'q', 'beta', 't')} if leg else None
    # ---- recommended rational interior point (certified by the isolating intervals alone when possible)
    info = ivinfo[0]; L = info.pop('_L'); a = info.pop('_a')
    for inf in ivinfo[1:]: inf.pop('_L'); inf.pop('_a')
    rec_pt = None
    iv0 = ivs[0]; hia, hib = map(F, iv0['hi']['beta_interval']); hp = iv0['hi']['beta_minpoly_hi_to_lo']
    if iv0['lo']['kind'] == 'root':
        loa, lob = map(F, iv0['lo']['beta_interval']); lp = iv0['lo']['beta_minpoly_hi_to_lo']
    else:
        loa = lob = F(2); lp = None
    nref = 0
    while not (lob < hia):
        hia, hib = refine_root_interval(hia, hib, hp)
        if lp is not None: loa, lob = refine_root_interval(loa, lob, lp)
        nref += 1
        assert nref < 200
    beta_rec = simplest_between(lob, hia)
    assert lob < beta_rec < hia
    cert_how = 'isolating intervals only' if nref == 0 else 'isolating intervals after %d exact bisection steps' % nref
    row['refined_intervals'] = {'lo': [str(loa), str(lob)], 'hi': [str(hia), str(hib)], 'bisection_steps': nref}
    row['beta_rec_strictly_inside'] = strictly_inside(beta_rec, ivs[0])
    if not row['beta_rec_strictly_inside']: probs.append('recommended beta not strictly inside')
    tried = []
    for t in [10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]:
        if t * (beta_rec - 2) <= 1: continue           # q <= 1 : illegal
        eps = F(1, t + 1); delta = (beta_rec - 1) * eps
        r = ctln_check(A, eps, delta)
        tried.append((t, r['stable'], r['positive']))
        if r['legal'] and r['positive'] and r['stable'] == 'stable':
            rec_pt = {'eps': r['eps'], 'delta': r['delta'], 'beta': str(beta_rec), 't': str(t), 'q': r['q'],
                      'beta_certified_by': cert_how, 'stable': 'stable', 'positive': True}
            break
    row['recommended_point'] = rec_pt; row['recommended_tried'] = [list(map(str, x)) for x in tried]
    if rec_pt is None: probs.append('no stable legal point found at the recommended beta')
    row['problems'] = probs
    return row

if __name__ == '__main__':
    t0 = time.time()
    lines = open(IN).read().splitlines()
    lines = [l for l in lines if l.strip()]
    print('records', len(lines), flush=True)
    with Pool(3) as pool:
        rows = []
        for i, row in enumerate(pool.imap(work, lines, chunksize=20)):
            rows.append(row)
            if (i + 1) % 1000 == 0: print(i + 1, '%.0fs' % (time.time() - t0), flush=True)
    print('workers done %.0fs' % (time.time() - t0), flush=True)
    # ---- pairwise non-isomorphism, own canonical form
    cc = collections.Counter(r['canon'] for r in rows)
    dup = {c: [r['d6'] for r in rows if r['canon'] == c] for c, k in cc.items() if k > 1}
    print('own canonical forms: distinct', len(cc), 'duplicates', len(dup))
    # ---- labelg second opinion
    open(OUT + 'input.d6', 'w').write('\n'.join(r['d6'] for r in rows) + '\n')
    labelg = {}
    try:
        subprocess.run([R + 'n7/labelg', '-q', OUT + 'input.d6', OUT + 'labelg_out.d6'], check=True, timeout=600)
        outl = open(OUT + 'labelg_out.d6').read().split()
        labelg = {'ran': True, 'n_lines': len(outl), 'n_distinct': len(set(outl))}
        # and labelg's classes must coincide with mine
        same = collections.defaultdict(set)
        for r, o in zip(rows, outl): same[o].add(r['canon'])
        labelg['classes_agree_with_own'] = all(len(s) == 1 for s in same.values()) and len(same) == len(cc)
    except Exception as e:
        labelg = {'ran': False, 'error': repr(e)}
    print('labelg', labelg)
    problems = [(r['d6'], p) for r in rows for p in r['problems']]
    for c, ds in dup.items(): problems.append((ds[0], 'isomorphic duplicates %s' % ds))
    json.dump({'rows': rows, 'own_canon_distinct': len(cc), 'duplicates': dup, 'labelg': labelg,
               'problems': problems, 'elapsed_s': time.time() - t0}, open(OUT + 'verify_out.json', 'w'))
    print('problems', len(problems), 'elapsed %.0fs' % (time.time() - t0))
