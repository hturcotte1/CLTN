"""n = 8 candidate verification, verifier v_fraction. Pure Fraction arithmetic, no project imports."""
import json, sys, os, time, random
from fractions import Fraction as F
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vf8 import *

REPO = '/home/user/CLTN'
OUT = os.path.dirname(os.path.abspath(__file__))
HITS = f'{REPO}/n8/n8_limit_HITS.jsonl'
RESULTS = f'{REPO}/n8/n8_limit_results.jsonl'
S = str
WIDTH = F(1, 10**9)
T_CHOICES = [10000, 100000, 1000000, 10**7]

def ed_from_bt(beta, t):
    eps = F(1, t + 1)
    return eps, (beta - 1) * eps

def point(A, eps, delta, label):
    r = analyze(A, eps, delta)
    q = q_of(eps, delta)
    row = dict(label=label, eps=S(eps), delta=S(delta), q=S(q), legal=r['legal'], positive=r['positive'],
               stable=r['stability'], minor_signs=[sgn(d) for d in r['minors']],
               min_x=S(min(r['x'])) if r['x'] else None)
    return row, r

def endpoint_bounds(ep, which):
    """Return (lower_bound, upper_bound, note) rationals enclosing the endpoint beta of a stable interval."""
    if ep['kind'] == 'beta_min':
        b = F(ep['beta']); return b, b, 'beta_min'
    assert ep['kind'] == 'root', ep
    if 'beta' in ep:
        b = F(ep['beta'])
        assert peval(ep['beta_minpoly_hi_to_lo'], b) == 0 if 'beta_minpoly_hi_to_lo' in ep else True
        return b, b, 'rational_root'
    a, b = F(ep['beta_interval'][0]), F(ep['beta_interval'][1])
    assert a <= b
    lo, hi, note = refine_root(ep['beta_minpoly_hi_to_lo'], a, b, WIDTH)
    if lo is None:
        return a, b, 'NO_SIGN_CHANGE_fallback_to_stored_interval'
    return lo, hi, note

def verify_one(line):
    rec = json.loads(line)
    d6 = rec['d6']
    ent = dict(d6=d6, failures=[], notes=[])
    try:
        A = decode_d6(d6)
        ent['n'] = len(A); ent['edges'] = sum(map(sum, A)); ent['loops'] = has_loops(A)
        ent['complete'] = is_complete_bidirected(A)
        if ent['n'] != 8: ent['failures'].append('n != 8')
        if ent['edges'] != rec['edges']: ent['failures'].append(f"edge count {ent['edges']} != stored {rec['edges']}")
        if ent['loops']: ent['failures'].append('loops')
        if ent['complete'] != rec['complete'] or ent['complete']: ent['failures'].append('complete flag')
        # ---- stored witnesses ----
        ws = []
        mism = []
        first_legal = None
        for k, w in enumerate(rec['witnesses']):
            beta, t = F(w['beta']), int(w['t'])
            eps, delta = ed_from_bt(beta, t)
            if (eps, delta) != (F(w['eps']), F(w['delta'])): ent['failures'].append(f'witness {k}: eps/delta not reproduced from beta,t')
            row, r = point(A, eps, delta, f'stored_witness_{k}')
            row['beta'] = S(beta); row['t'] = t
            if F(w['q']) != q_of(eps, delta): ent['failures'].append(f'witness {k}: stored q mismatch')
            row['stored_stable'] = w['stable']; row['stored_positive'] = w['positive']; row['stored_detail'] = w['detail']
            agree_pos = (row['positive'] == w['positive'])
            agree_stab = None if row['stable'] == 'undecided' else ((row['stable'] == 'stable') == bool(w['stable']))
            row['agree_positive'] = agree_pos; row['agree_stable'] = agree_stab
            if not agree_pos or agree_stab is False:
                mism.append(k); ent['failures'].append(f'witness {k}: disagreement mine=({row["positive"]},{row["stable"]}) stored=({w["positive"]},{w["stable"]})')
            if row['legal'] and first_legal is None: first_legal = k
            ws.append(row)
        ent['witnesses'] = ws
        ent['witness_mismatches'] = mism
        ent['first_witness_illegal'] = (not ws[0]['legal']) if ws else None
        ent['first_legal_witness_index'] = first_legal
        # ---- (b) interior point of the stable beta-interval at large t ----
        si = rec['cert']['stable_intervals']
        assert len(si) == 1, 'expected exactly one stable interval'
        lo_lb, lo_ub, lo_note = endpoint_bounds(si[0]['lo'], 'lo')
        hi_lb, hi_ub, hi_note = endpoint_bounds(si[0]['hi'], 'hi')
        ent['interval'] = dict(lo_bounds=[S(lo_lb), S(lo_ub)], lo_note=lo_note, hi_bounds=[S(hi_lb), S(hi_ub)], hi_note=hi_note)
        if 'NO_SIGN_CHANGE' in lo_note or 'NO_SIGN_CHANGE' in hi_note: ent['notes'].append('endpoint minpoly had no sign change on stored interval')
        if not (lo_ub < hi_lb):
            ent['failures'].append(f'cannot separate interval endpoints: lo<= {lo_ub}, hi>= {hi_lb}')
            ent['point_b'] = None
        else:
            chosen = None
            for t in T_CHOICES:
                L = max(lo_ub, 2 + F(1, t))   # q = (beta-1)*t/(t+1) > 1  iff  beta > 2 + 1/t
                H = hi_lb
                if L < H:
                    g = (H - L) / 4
                    beta = simplest_between(L + g, H - g)
                    assert L < beta < H and lo_ub < beta < hi_lb
                    chosen = (beta, t); break
            if chosen is None:
                ent['failures'].append('no legal window inside the stable interval for any t in T_CHOICES'); ent['point_b'] = None
            else:
                beta, t = chosen
                eps, delta = ed_from_bt(beta, t)
                row, r = point(A, eps, delta, 'interior_point')
                row['beta'] = S(beta); row['t'] = t
                row['strictly_inside'] = bool(lo_ub < beta < hi_lb)
                if not row['legal']: ent['failures'].append('interior point illegal')
                ent['point_b'] = row
        # ---- (a) first stored witness with q > 1 ----
        ent['point_a'] = ws[first_legal] if first_legal is not None else None
        if first_legal is None: ent['notes'].append('no stored witness with q > 1')
        # ---- verdict ----
        # CONFIRMED: my interior point (b) is legal, x > 0 and Hurwitz-stable, every stored witness is re-decided exactly
        #   with the same stable/positive flags, the graph is not complete and no consistency check failed.
        # REFUTED: a decided contradiction (point (b) decided unstable or non-positive, a stored flag contradicted, graph complete).
        # UNDECIDED: otherwise (a vanishing Hurwitz minor at (b), or an exception / unreproducible point).
        # Point (a) (first stored witness with q > 1) is reported separately as point_a_verdict: a stored witness that the
        # file itself marks unstable is a finite-t point, not a claim, so an exact 'unstable' there that agrees with the
        # stored flag is not a refutation.
        def ok(row): return row is not None and row['legal'] and row['positive'] and row['stable'] == 'stable'
        def bad(row): return row is not None and row['legal'] and (not row['positive'] or row['stable'] == 'unstable')
        pa, pb = ent['point_a'], ent['point_b']
        ent['point_a_verdict'] = None if pa is None else ('stable_positive' if ok(pa) else ('unstable_or_nonpositive' if bad(pa) else 'undecided'))
        contradiction = ent['complete'] or any('disagreement' in f or 'edge count' in f for f in ent['failures'])
        if contradiction or bad(pb):
            v = 'REFUTED'
        elif ok(pb) and not ent['failures']:
            v = 'CONFIRMED'
        else:
            v = 'UNDECIDED'
        if v == 'CONFIRMED' and pa is None:
            ent['notes'].append('no legal stored witness; confirmed on my interior point')
        if v == 'CONFIRMED' and pa is not None and not ok(pa):
            ent['notes'].append('first legal stored witness is decided unstable, in agreement with its stored flag (finite-t point); confirmed on my interior point')
        ent['verdict'] = v
    except Exception as e:
        ent['failures'].append(f'EXCEPTION {type(e).__name__}: {e}'); ent['verdict'] = 'UNDECIDED'
    return ent

def controls():
    out = {}
    # decoder checks
    G = decode_d6('&EZjtwWE')
    out['geneson_decode'] = dict(n=len(G), arcs=sum(map(sum, G)), roundtrip=encode_d6(G) == '&EZjtwWE')
    g = analyze(G, F(1, 100), F(29, 2500))
    out['geneson_n6_control'] = dict(legal=g['legal'], positive=g['positive'], stability=g['stability'])
    K8 = [[int(i != j) for j in range(8)] for i in range(8)]
    out['K8_d6'] = encode_d6(K8)
    k = []
    for eps, delta, lab in [(F(1, 4), F(1, 2), 'eps1/4_delta1/2'), (F(1, 10001), F(3, 2) * F(1, 10001), 'ray_beta5/2_t10000'),
                            (F(1, 101), F(1001, 101000), 'ray_beta2001/1000_t100_ILLEGAL'), (F(1, 3), F(2), 'eps1/3_delta2')]:
        row, _ = point(K8, eps, delta, lab); k.append(row)
    out['K8_control'] = k
    out['K8_control_ok'] = all(r['positive'] and r['stable'] == 'stable' for r in k)
    C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    c = analyze(C3, F(1, 4), F(1, 2))
    out['C3_control_expect_unstable'] = dict(positive=c['positive'], stability=c['stability'])
    # 20 random 'never' strict survivors
    nev = []
    for l in open(RESULTS):
        r = json.loads(l)
        if r['verdict'] == 'never': nev.append((r['d6'], r.get('reason')))
    rng = random.Random(20261027)
    sample = rng.sample(nev, 20)
    rows = []
    n_stable_pos_legal = 0; n_pts = 0; n_unst_or_nonpos = 0; n_und = 0
    betas = [F(2) + F(1, 1000), F(21, 10), F(5, 2), F(3), F(4), F(6), F(10)]
    for d6, reason in sample:
        A = decode_d6(d6)
        pts = []
        for t in (10000, 100000):
            for beta in betas:
                eps, delta = ed_from_bt(beta, t)
                row, _ = point(A, eps, delta, f'ray_beta{beta}_t{t}')
                n_pts += 1
                if row['legal'] and row['positive'] and row['stable'] == 'stable': n_stable_pos_legal += 1
                elif not row['positive'] or row['stable'] == 'unstable': n_unst_or_nonpos += 1
                else: n_und += 1
                pts.append(dict(beta=S(beta), t=t, q=row['q'], legal=row['legal'], positive=row['positive'], stable=row['stable']))
        rows.append(dict(d6=d6, reason=reason, edges=sum(map(sum, A)), complete=is_complete_bidirected(A), points=pts,
                         any_stable_positive_legal=any(p['legal'] and p['positive'] and p['stable'] == 'stable' for p in pts)))
    out['never_sample'] = dict(seed=20261027, n_never_total=len(nev), n_sampled=len(sample), points_per_graph=len(betas) * 2,
                               betas=[S(b) for b in betas], ts=[10000, 100000], points_total=n_pts,
                               points_unstable_or_nonpositive=n_unst_or_nonpos, points_undecided=n_und,
                               points_stable_positive_legal=n_stable_pos_legal,
                               graphs_with_a_stable_positive_legal_point=sum(r['any_stable_positive_legal'] for r in rows), rows=rows)
    return out

if __name__ == '__main__':
    t0 = time.time()
    lines = open(HITS).read().splitlines()
    print('candidates:', len(lines), flush=True)
    ctl = controls()
    print('controls done', time.time() - t0, json.dumps({k: v for k, v in ctl.items() if k != 'never_sample'}), flush=True)
    print('never sample:', json.dumps({k: v for k, v in ctl['never_sample'].items() if k != 'rows'}), flush=True)
    results = []
    with Pool(3) as pool, open(f'{OUT}/per_graph.jsonl', 'w') as fo:
        for i, ent in enumerate(pool.imap(verify_one, lines, chunksize=20)):
            results.append(ent); fo.write(json.dumps(ent) + '\n')
            if (i + 1) % 500 == 0: print(i + 1, 'done', round(time.time() - t0, 1), 's', flush=True)
    assert len(results) == len(lines) and [r['d6'] for r in results] == [json.loads(l)['d6'] for l in lines]
    from collections import Counter
    verd = Counter(r['verdict'] for r in results)
    table = {}
    for r in results:
        pb = r.get('point_b'); pa = r.get('point_a')
        table[r['d6']] = dict(point=dict(beta=pb['beta'], t=pb['t'], eps=pb['eps'], delta=pb['delta'], q=pb['q']) if pb else None,
                              x_positive=pb['positive'] if pb else None, hurwitz=pb['minor_signs'] if pb else None,
                              hurwitz_verdict=pb['stable'] if pb else None,
                              first_legal_witness=dict(index=r['first_legal_witness_index'], q=pa['q'], x_positive=pa['positive'], hurwitz=pa['stable']) if pa else None,
                              first_witness_illegal=r['first_witness_illegal'], witness_mismatches=r['witness_mismatches'],
                              verdict=r['verdict'], failures=r['failures'], notes=r['notes'])
    json.dump(table, open(f'{OUT}/v_fraction_n8_table.json', 'w'), indent=0)
    npts = sum(len(r.get('witnesses', [])) + (1 if r.get('point_b') else 0) for r in results)
    summary = dict(verifier='v_fraction', n_candidates=len(results), distinct_d6=len({r['d6'] for r in results}),
                   verdicts=dict(verd), points_checked=npts,
                   witness_points=sum(len(r.get('witnesses', [])) for r in results), interior_points=sum(1 for r in results if r.get('point_b')),
                   first_witness_illegal=sum(1 for r in results if r.get('first_witness_illegal')),
                   no_legal_stored_witness=sum(1 for r in results if r.get('point_a') is None),
                   witness_disagreements=sum(len(r.get('witness_mismatches', [])) for r in results),
                   witness_agree_stable_both_decided=sum(1 for r in results for w in r.get('witnesses', []) if w['agree_stable'] is True),
                   witness_undecided_mine=sum(1 for r in results for w in r.get('witnesses', []) if w['stable'] == 'undecided'),
                   stored_unstable_witnesses=sum(1 for r in results for w in r.get('witnesses', []) if not w['stored_stable']),
                   my_unstable_witnesses=sum(1 for r in results for w in r.get('witnesses', []) if w['stable'] == 'unstable'),
                   interior_t_used=dict(Counter(r['point_b']['t'] for r in results if r.get('point_b'))),
                   interior_strictly_inside=sum(1 for r in results if r.get('point_b') and r['point_b']['strictly_inside']),
                   endpoint_notes=dict(Counter(r['interval']['lo_note'] + '|' + r['interval']['hi_note'] for r in results if 'interval' in r)),
                   graphs_with_failures=sum(1 for r in results if r['failures']),
                   not_confirmed=[r['d6'] for r in results if r['verdict'] != 'CONFIRMED'],
                   seconds=round(time.time() - t0, 1))
    json.dump(dict(summary=summary, controls=ctl), open(f'{OUT}/v_fraction_n8_results.json', 'w'), indent=1)
    print(json.dumps(summary, indent=1))
    for r in results:
        if r['failures'] or r['verdict'] != 'CONFIRMED': print('NOTCONF', r['d6'], r['verdict'], r['failures'], r['notes'])
