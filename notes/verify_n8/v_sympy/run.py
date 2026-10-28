"""v_sympy n = 8 driver.  Usage: python3 -I run.py <input.jsonl> <out_dir> <chunk_idx> <n_chunks>
Sample rules (0-based line index i in the input file):
  every graph        : decode, arcs == edges, not complete, no loops; exact check at the FIRST LEGAL stored witness
                       (legal iff q = delta(1-eps)/eps > 1 recomputed exactly); Sturm root count on the stored
                       hi-endpoint isolating interval (must be exactly 1).
  i % 10 == 0        : re-decide ALL FOUR stored witnesses and compare stable/positive with stored fields.
  i % 50 == 0        : sharpness of the hi endpoint: refine the root by exact Sturm bisection, then at t = 100000
                       beta_out = r_hi + 1/1000 must NOT be stable-and-positive and beta_in = r_lo - d_in
                       (d_in = min(1/1000, (r_lo - lo)/2)) must be stable-and-positive.
  first legal witness not stable-and-positive: check all four witnesses and, if none confirms, an OWN witness
                       beta = (lo + r_lo)/2 at t = 100000 (then 10^6, 10^7) inside the certified interval.
"""
import sys, json, time, os
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vsympy8 import *
from sympy import Rational

INP, OUT, CH, NCH = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

def lo_value(lo):
    if lo['kind'] == 'beta_min':
        return Rational(lo['beta'])
    assert lo['kind'] == 'root', lo['kind']
    return Rational(lo['beta_interval'][1])      # upper bound of the lo root's isolating interval: inside is > this

def brief(res):
    return {'legal': res['legal'], 'positive': res['positive'], 'stable': res['stable'],
            'routh_agrees': res['routh_agrees'], 'q': res['q']}

def sp(res):
    return res['legal'] is True and res['positive'] is True and res['stable'] is True

def work(item):
    i, line = item
    t0 = time.time()
    r = json.loads(line)
    d6 = r['d6']
    row = {'idx': i, 'd6': d6}
    try:
        A = decode_d6(d6)
        n = len(A)
        row['n'] = n
        row['arcs'] = sum(map(sum, A))
        row['arcs_match_record'] = (row['arcs'] == r['edges'])
        row['loops'] = sum(A[k][k] for k in range(n))
        row['complete'] = is_complete_bidirected(A)
        row['complete_match_record'] = (row['complete'] == r['complete'])
        row['reencode_ok'] = (encode_d6(A) == d6)
        row['sigma_full'] = (r['cert']['sigma'] == list(range(n)))
        W = r['witnesses']
        # ---- first legal witness ----
        fl = None
        legal_flags = []
        for k, w in enumerate(W):
            eps, delta = Rational(w['eps']), Rational(w['delta'])
            q = delta * (1 - eps) / eps
            assert q == Rational(w['q']), ('stored q mismatch', d6, k)
            # witness consistency with the limit ray
            assert eps == Rational(1, w['t'] + 1) and delta == (Rational(w['beta']) - 1) * eps, ('witness ray mismatch', d6, k)
            leg = bool(eps > 0 and eps < 1 and delta > 0 and q > 1)
            legal_flags.append(leg)
            if leg and fl is None:
                fl = k
        row['witness_legal'] = legal_flags
        row['first_legal_idx'] = fl
        if fl is None:
            row['first_legal'] = None
            row['verdict_first_legal'] = 'NO_LEGAL_WITNESS'
        else:
            w = W[fl]
            res = check(A, w['eps'], w['delta'], crosscheck=True)
            b = brief(res); b['t'] = w['t']; b['beta'] = w['beta']
            b['closed_form_ok'] = res['closed_form_ok']; b['charpoly_crosscheck'] = res['charpoly_crosscheck']
            b['hurwitz_signs'] = res['hurwitz_minors_sign']; b['x'] = res['x']
            b['stored'] = {'stable': w['stable'], 'positive': w['positive']}
            b['agrees_with_stored'] = (res['stable'] == w['stable'] and res['positive'] == w['positive'])
            row['first_legal'] = b
            if res['stable'] == 'UNDECIDED' or res['positive'] is None:
                row['verdict_first_legal'] = 'UNDECIDED'
            elif sp(res) and not row['complete'] and res['closed_form_ok'] and res['charpoly_crosscheck'] and res['routh_agrees'] is True:
                row['verdict_first_legal'] = 'COUNTEREXAMPLE_CONFIRMED'
            else:
                row['verdict_first_legal'] = 'NOT_STABLE_POSITIVE_HERE'
        # ---- Sturm count on the stored isolating interval (every graph) ----
        si = r['cert']['stable_intervals']
        assert len(si) == 1
        lo, hi = si[0]['lo'], si[0]['hi']
        assert hi['kind'] == 'root'
        lov = lo_value(lo)
        a, b_ = Rational(hi['beta_interval'][0]), Rational(hi['beta_interval'][1])
        P, seq = sturm_seq(hi['beta_minpoly_hi_to_lo'])
        row['hi_minpoly_degree'] = P.degree()
        if a == b_:
            row['hi_root_rational'] = True
            row['sturm_count_hi'] = 1 if P.eval(a) == 0 else 0
            r_lo = r_hi = a
        else:
            row['hi_root_rational'] = False
            if P.eval(a) == 0 or P.eval(b_) == 0:
                row['sturm_count_hi'] = 'ENDPOINT_IS_ROOT'
                r_lo = r_hi = None
            else:
                row['sturm_count_hi'] = count_roots_open(P, seq, a, b_)
                r_lo, r_hi = None, None
        row['interval_lo_value'] = str(lov)
        row['interval_nonempty'] = bool(a >= lov)
        # ---- all four witnesses (every 10th, and whenever the first legal witness did not confirm) ----
        need_all = (i % 10 == 0) or (row['verdict_first_legal'] != 'COUNTEREXAMPLE_CONFIRMED')
        if need_all:
            allw = []
            for k, w in enumerate(W):
                res = check(A, w['eps'], w['delta'], crosscheck=False)
                b = brief(res); b['t'] = w['t']
                b['stored'] = {'stable': w['stable'], 'positive': w['positive']}
                b['agrees_with_stored'] = (res['stable'] == w['stable'] and res['positive'] == w['positive'])
                b['confirms'] = bool(sp(res) and not row['complete'])
                allw.append(b)
            row['all_witnesses'] = allw
            row['all_witnesses_agree_with_stored'] = all(x['agrees_with_stored'] for x in allw)
            row['some_stored_witness_confirms'] = any(x['confirms'] for x in allw)
        # ---- refinement (needed for sharpness sample and for own witness) ----
        need_refine = (i % 50 == 0) or (row['verdict_first_legal'] != 'COUNTEREXAMPLE_CONFIRMED' and not row.get('some_stored_witness_confirms', False))
        if need_refine and row['sturm_count_hi'] == 1:
            if r_lo is None:
                r_lo, r_hi = refine_root(P, seq, a, b_, iters=40)
            row['hi_root_refined'] = [str(r_lo), str(r_hi)]
            row['hi_root_refined_float_display_only'] = [float(r_lo), float(r_hi)]
        if i % 50 == 0 and row['sturm_count_hi'] == 1:
            d_out = Rational(1, 1000)
            d_in = min(Rational(1, 1000), (r_lo - lov) / 2)
            sh = {'d_out': str(d_out), 'd_in': str(d_in)}
            T = 100000
            bo = r_hi + d_out
            eps, delta = point_from_beta(bo, T)
            ro = check(A, eps, delta, crosscheck=False)
            sh['beyond'] = {'beta': str(bo), 't': T, **brief(ro)}
            sh['beyond_not_stable_positive'] = (not sp(ro)) and ro['stable'] != 'UNDECIDED'
            if sp(ro):   # possibly a later stable cell begins within d_out: shrink
                for d2, T2 in ((Rational(1, 10000), 1000000), (Rational(1, 100000), 10000000)):
                    bo2 = r_hi + d2
                    ro2 = check(A, *point_from_beta(bo2, T2), crosscheck=False)
                    sh.setdefault('beyond_retries', []).append({'beta': str(bo2), 't': T2, **brief(ro2)})
                    if not sp(ro2) and ro2['stable'] != 'UNDECIDED':
                        sh['beyond_not_stable_positive'] = True
                        break
            bi = r_lo - d_in
            ri = check(A, *point_from_beta(bi, T), crosscheck=False)
            sh['inside'] = {'beta': str(bi), 't': T, **brief(ri)}
            sh['inside_stable_positive'] = sp(ri)
            if not sp(ri):
                for T2 in (1000000, 10000000):
                    ri2 = check(A, *point_from_beta(bi, T2), crosscheck=False)
                    sh.setdefault('inside_retries', []).append({'beta': str(bi), 't': T2, **brief(ri2)})
                    if sp(ri2):
                        sh['inside_stable_positive'] = True
                        break
            sh['sharp'] = bool(sh['beyond_not_stable_positive'] and sh['inside_stable_positive'])
            row['sharpness'] = sh
        # ---- own witness for graphs with no confirming stored legal witness ----
        if row['verdict_first_legal'] != 'COUNTEREXAMPLE_CONFIRMED' and not row.get('some_stored_witness_confirms', False):
            own = []
            ok = False
            if row['sturm_count_hi'] == 1 and r_lo > lov:
                bw = (lov + r_lo) / 2
                for T in (100000, 1000000, 10000000):
                    eps, delta = point_from_beta(bw, T)
                    res = check(A, eps, delta, crosscheck=True)
                    e = {'beta': str(bw), 't': T, **brief(res), 'closed_form_ok': res['closed_form_ok'],
                         'charpoly_crosscheck': res['charpoly_crosscheck'], 'x': res['x']}
                    own.append(e)
                    if sp(res) and res['closed_form_ok'] and res['charpoly_crosscheck'] and res['routh_agrees'] is True and not row['complete']:
                        ok = True
                        break
            row['own_witness'] = own
            row['own_witness_confirms'] = ok
        # ---- overall verdict ----
        if row['verdict_first_legal'] == 'COUNTEREXAMPLE_CONFIRMED':
            row['verdict'] = 'CONFIRMED_AT_STORED_FIRST_LEGAL_WITNESS'
        elif row.get('some_stored_witness_confirms'):
            row['verdict'] = 'CONFIRMED_AT_LATER_STORED_WITNESS'
        elif row.get('own_witness_confirms'):
            row['verdict'] = 'CONFIRMED_AT_OWN_WITNESS'
        elif row['verdict_first_legal'] == 'UNDECIDED' or any(x['stable'] == 'UNDECIDED' for x in row.get('all_witnesses', [])):
            row['verdict'] = 'UNDECIDED'
        else:
            row['verdict'] = 'NOT_CONFIRMED'
        row['integrity_ok'] = bool(row['arcs_match_record'] and row['complete_match_record'] and not row['complete']
                                   and row['loops'] == 0 and row['reencode_ok'] and row['sigma_full']
                                   and row['sturm_count_hi'] == 1 and row['interval_nonempty'])
    except Exception as ex:
        row['verdict'] = 'ERROR'
        row['error'] = repr(ex)
    row['seconds'] = round(time.time() - t0, 3)
    return row

def main():
    lines = open(INP).read().splitlines()
    items = [(i, l) for i, l in enumerate(lines) if i % NCH == CH]
    out = open(f'{OUT}/table_part{CH}.jsonl', 'w')
    t0 = time.time()
    with Pool(3) as pool:
        for k, row in enumerate(pool.imap(work, items, chunksize=4)):
            out.write(json.dumps(row) + '\n')
            if k % 200 == 0 or row['verdict'] not in ('CONFIRMED_AT_STORED_FIRST_LEGAL_WITNESS',):
                print(k, row['idx'], row['d6'], row['verdict'], row.get('verdict_first_legal'), row['seconds'], f'{time.time()-t0:.0f}s', flush=True)
    out.close()
    print('DONE chunk', CH, len(items), f'{time.time()-t0:.0f}s', flush=True)

if __name__ == '__main__':
    main()
