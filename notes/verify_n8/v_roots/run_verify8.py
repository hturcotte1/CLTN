"""v_roots verifier for the n = 8 limit candidates.  Per graph:
  P1: first stored legal witness (q > 1), recomputed from scratch;
  P2: own interior point of the certified stable beta-interval at t = 10000 (eps = 1/10001, delta = (beta-1) eps).
Cheap exact methods (Routh, Cauchy index) on every point; sympy isolation additionally on every K-th graph (index % K == 0).
"""
import sys, json, time, os
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8/v_roots')
from fractions import Fraction as F
from multiprocessing import Pool
from vroots8 import d6_decode, d6_encode, is_complete_bidirected, check_point, refine_root

K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
T2 = 10000
IN = '/home/user/CLTN/n8/n8_limit_HITS.jsonl'
OUT = '/home/user/CLTN/notes/verify_n8/v_roots/'

def tdpt(beta, t):
    eps = F(1, t + 1); delta = (beta - 1) * eps
    return eps, delta

def summarize(tag, p):
    keep = ('eps', 'delta', 'q', 'legal', 'det_IminusW', 'nondegenerate', 'positive', 'stable',
            'rhp_routh', 'rhp_cauchy', 'rhp_sympy', 'sympy_run', 'methods_agree', 'counterexample_point')
    d = {'point': tag}; d.update({k: p.get(k) for k in keep})
    for k in ('routh_undecided_first_col', 'cauchy_note', 'sympy_note'):
        if k in p: d[k] = p[k]
    return d

def do_graph(args):
    idx, line = args
    r = json.loads(line)
    wits0 = r['witnesses']
    no_claimed = not any(F(w['delta']) * (1 - F(w['eps'])) / F(w['eps']) > 1 and w['stable'] and w['positive'] for w in wits0)
    use_sympy = (idx % K == 0) or no_claimed      # sympy on every K-th graph, and on every graph whose stored witnesses never claim a legal stable point
    A = d6_decode(r['d6'])
    row = {'idx': idx, 'd6': r['d6'], 'edges_decoded': sum(map(sum, A)), 'edges_stored': r['edges'],
           'complete_bidirected': is_complete_bidirected(A), 'sympy_sample': use_sympy, 'points': [], 'issues': []}
    try:
        if d6_encode(A) != r['d6']: row['issues'].append('d6 does not round-trip')
        if row['edges_decoded'] != r['edges']: row['issues'].append('edge count mismatch')
        if row['complete_bidirected'] or r['complete']: row['issues'].append('complete bidirected graph (not a counterexample)')
        if r['cert']['sigma'] != list(range(8)): row['issues'].append('sigma is not the full support')
        # ---- P1: first legal stored witness ----
        wits = r['witnesses']
        legal_flags = []
        first = None
        for w in wits:
            eps, delta = F(w['eps']), F(w['delta'])
            e2, d2 = tdpt(F(w['beta']), int(w['t']))
            if (e2, d2) != (eps, delta): row['issues'].append('witness t=%s: beta/t inconsistent with eps/delta' % w['t'])
            q = delta * (1 - eps) / eps
            if q != F(w['q']): row['issues'].append('witness t=%s: stored q differs from delta(1-eps)/eps' % w['t'])
            legal_flags.append(q > 1)
            if q > 1 and first is None: first = w
        row['witness_legal_flags'] = legal_flags
        row['witness_stored_stable_flags'] = [bool(w['stable']) for w in wits]
        row['n_witnesses_illegal'] = sum(1 for f in legal_flags if not f)
        row['first_stored_witness_illegal'] = not legal_flags[0]
        first_claimed = next((w for w in wits if F(w['delta']) * (1 - F(w['eps'])) / F(w['eps']) > 1 and w['stable'] and w['positive']), None)
        row['has_stored_legal_stable_witness'] = first_claimed is not None
        if first is None:
            row['issues'].append('no legal stored witness')
        else:
            row['first_legal_witness_t'] = int(first['t'])
            todo = [('witness', first)]
            if first_claimed is not None and first_claimed is not first:
                todo.append(('claimed_witness', first_claimed))
            for tag, w in todo:
                p = check_point(A, F(w['eps']), F(w['delta']), use_sympy=use_sympy)
                s = summarize('%s t=%s beta=%s' % (tag, w['t'], w['beta']), p)
                s['stored_stable_positive'] = [w.get('stable'), w.get('positive')]
                if (p['stable'], p['positive']) != (w.get('stable'), w.get('positive')):
                    row['issues'].append('stored %s t=%s flags (stable,positive)=%s != recomputed %s' % (tag, w['t'], (w.get('stable'), w.get('positive')), (p['stable'], p['positive'])))
                if not p['nondegenerate']: row['issues'].append('det(I-W)=0 at %s' % tag)
                if not p['methods_agree']: row['issues'].append('methods disagree at %s' % tag)
                row['points'].append(s)
        # ---- P2: own interior point of the certified stable interval at t = T2 ----
        si = r['cert']['stable_intervals']
        if len(si) != 1: row['issues'].append('expected exactly one stable interval, got %d' % len(si))
        iv = si[0]
        lo, hi = iv['lo'], iv['hi']
        assert hi['kind'] == 'root', hi['kind']
        beta_legal = 1 + F(1) / (1 - F(1, T2 + 1))        # q > 1 at t = T2  <=>  beta > beta_legal (= 2 + 1/T2)
        if lo['kind'] == 'beta_min':
            lo_cert = F(lo['beta']); lo_desc = 'beta_min=%s' % lo['beta']
        else:
            assert lo['kind'] == 'root', lo['kind']
            la, lb = F(lo['beta_interval'][0]), F(lo['beta_interval'][1])
            la, lb = refine_root(lo['beta_minpoly_hi_to_lo'], la, lb, lambda a, b: b - a < F(1, 10**6))
            lo_cert = lb; lo_desc = 'root in [%s,%s]' % (la, lb)
        lo_use = max(lo_cert, beta_legal)
        ha, hb = F(hi['beta_interval'][0]), F(hi['beta_interval'][1])
        ha, hb = refine_root(hi['beta_minpoly_hi_to_lo'], ha, hb,
                             lambda a, b: (a - lo_use > (b - a)) or (b - a < F(1, 10**30)))
        hi_cert = ha
        row['interval_certified'] = [str(lo_use), str(hi_cert)]
        row['interval_lo'] = lo_desc; row['interval_hi'] = 'root in [%s,%s]' % (ha, hb)
        wb = F(wits[0]['beta'])
        row['stored_witness_beta_inside_certified_limit_interval'] = bool(lo_cert < wb < hi_cert)
        if hi_cert <= lo_use:
            row['issues'].append('no legal interior point at t=%d: certified interval empty (%s, %s)' % (T2, lo_use, hi_cert))
        else:
            beta2 = lo_use + (hi_cert - lo_use) * F(3, 7)
            if beta2 == wb: beta2 = lo_use + (hi_cert - lo_use) * F(4, 9)
            assert lo_use < beta2 < hi_cert and beta2 > beta_legal
            eps, delta = tdpt(beta2, T2)
            p = check_point(A, eps, delta, use_sympy=use_sympy)
            s = summarize('interior t=%d beta=%s' % (T2, beta2), p)
            if not p['legal']: row['issues'].append('interior point illegal (should not happen)')
            if not p['counterexample_point']: row['issues'].append('interior point does not check: legal=%s positive=%s stable=%s' % (p['legal'], p['positive'], p['stable']))
            if not p['nondegenerate']: row['issues'].append('det(I-W)=0 at interior point')
            if not p['methods_agree']: row['issues'].append('methods disagree at interior point')
            row['points'].append(s)
    except Exception as e:
        row['issues'].append('EXCEPTION %r' % (e,))
    pts = row['points']
    # stored witnesses that CLAIM stability (stored stable & positive) must be verified counterexample points;
    # a stored witness that the producing stage itself recorded as unstable is only checked for consistency.
    claimed = [s for s in pts if not s['point'].startswith('interior') and s.get('stored_stable_positive') == [True, True]]
    interior = [s for s in pts if s['point'].startswith('interior')]
    claimed_ok = all(s['counterexample_point'] for s in claimed)
    interior_ok = len(interior) == 1 and interior[0]['counterexample_point']
    decided_false = lambda s: (s['stable'] is False) or (s['positive'] is False)
    row['n_verified_counterexample_points'] = sum(s['counterexample_point'] for s in pts)
    if claimed_ok and interior_ok and not row['issues']:
        row['verdict'] = 'CONFIRMED'
    elif row['complete_bidirected'] or any(decided_false(s) for s in claimed) or (interior and decided_false(interior[0]) and row['n_verified_counterexample_points'] == 0):
        row['verdict'] = 'REFUTED'
    else:
        row['verdict'] = 'UNDECIDED'
    return row

if __name__ == '__main__':
    lines = [l for l in open(IN) if l.strip()]
    print('graphs', len(lines), 'K', K, 'pid', os.getpid()); sys.stdout.flush()
    t0 = time.time()
    rows = []
    with Pool(3) as pool, open(OUT + 'v_roots_n8_table.jsonl', 'w') as fo:
        for k, row in enumerate(pool.imap(do_graph, list(enumerate(lines)), chunksize=4)):
            rows.append(row)
            fo.write(json.dumps(row) + '\n')
            if row['verdict'] != 'CONFIRMED' or k % 200 == 0:
                print(k, row['d6'], row['verdict'], row['issues'], '%.0fs' % (time.time() - t0)); sys.stdout.flush()
    rows.sort(key=lambda r: r['idx'])
    pts = [s for r in rows for s in r['points']]
    summary = {
        'verifier': 'v_roots', 'K_sympy_every_kth_graph': K, 'T2_interior_point': T2,
        'n_graphs': len(rows),
        'n_confirmed': sum(r['verdict'] == 'CONFIRMED' for r in rows),
        'n_refuted': sum(r['verdict'] == 'REFUTED' for r in rows),
        'n_undecided': sum(r['verdict'] == 'UNDECIDED' for r in rows),
        'n_points_checked': len(pts),
        'n_points_counterexample': sum(s['counterexample_point'] for s in pts),
        'n_points_undecided': sum(s['stable'] is None for s in pts),
        'n_points_method_disagreement': sum(not s['methods_agree'] for s in pts),
        'n_points_degenerate': sum(not s['nondegenerate'] for s in pts),
        'n_points_sympy_run': sum(bool(s['sympy_run']) for s in pts),
        'n_graphs_sympy_sample': sum(r['sympy_sample'] for r in rows),
        'n_points_repeated_roots': sum(1 for s in pts if 'repeated' in str(s.get('sympy_note'))),
        'n_graphs_complete': sum(r['complete_bidirected'] for r in rows),
        'n_distinct_d6': len(set(r['d6'] for r in rows)),
        'first_legal_witness_t_histogram': {},
        'n_stored_witnesses_illegal': sum(r.get('n_witnesses_illegal', 0) for r in rows),
        'n_graphs_with_illegal_stored_witness': sum(1 for r in rows if r.get('n_witnesses_illegal', 0) > 0),
        'n_graphs_first_stored_witness_illegal': sum(1 for r in rows if r.get('first_stored_witness_illegal')),
        'n_graphs_without_stored_legal_stable_witness': sum(1 for r in rows if r.get('has_stored_legal_stable_witness') is False),
        'n_graphs_first_legal_witness_stored_unstable': sum(1 for r in rows if r['points'] and r['points'][0].get('stored_stable_positive') == [False, True]),
        'n_points_by_tag': {},
        'n_witness_beta_outside_certified_limit_interval': sum(1 for r in rows if r.get('stored_witness_beta_inside_certified_limit_interval') is False),
        'wall_seconds': round(time.time() - t0, 1),
        'not_confirmed': [(r['d6'], r['verdict'], r['issues']) for r in rows if r['verdict'] != 'CONFIRMED'],
    }
    for s in pts:
        tag = s['point'].split()[0]
        summary['n_points_by_tag'][tag] = summary['n_points_by_tag'].get(tag, 0) + 1
    for r in rows:
        t = str(r.get('first_legal_witness_t'))
        summary['first_legal_witness_t_histogram'][t] = summary['first_legal_witness_t_histogram'].get(t, 0) + 1
    json.dump(summary, open(OUT + 'v_roots_n8_summary.json', 'w'), indent=1)
    print(json.dumps(summary, indent=1))
