"""driver: per-candidate exact all-t certification. Usage: run.py OUT WALLSECONDS [WORKERS]"""
import sys, json, time, os, traceback
from fractions import Fraction as F
from multiprocessing import Pool
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/interval_all_t')
import engine as E

HITS = '/home/user/CLTN/n8/n8_limit_HITS.jsonl'
CW = '/home/user/CLTN/certificates/n8_candidates_legal_witnesses.json'
TS = [100, 1000, 10000, 100000]

def fs(x):
    return str(x)

def root_rec(root, facts):
    lo, hi, fi = root
    return {'lo': fs(lo), 'hi': fs(hi), 'factor_index': fi}

def analyse_beta(A, beta, tag):
    """symbolic-in-t analysis at a fixed rational beta; returns dict (all strings exact)."""
    conds, extra = E.symbolic_conditions(A, beta)
    # set including legality
    res = E.sign_set(dict(conds))
    # set without legality (stability + positivity only), for comparing with stored stable/positive flags
    conds_sp = {k: v for k, v in conds.items() if k != 'L'}
    res_sp = E.sign_set(dict(conds_sp))
    def pack(res):
        roots = [dict(root_rec(r, res['facts']), vanishing=v) for r, v in zip(res['roots'], res['root_conds'])]
        S = []
        for left, right in res['S']:
            S.append({'left': None if left is None else root_rec(left, res['facts']),
                      'right': None if right is None else root_rec(right, res['facts'])})
        unbounded = bool(res['S']) and res['S'][-1][1] is None
        # t* = largest positive root of any condition polynomial (None if no positive roots)
        tstar = root_rec(res['roots'][-1], res['facts']) if res['roots'] else None
        # the set beyond the largest root is all-positive iff unbounded
        return {'factors_hi_to_lo': [[int(c) for c in f.all_coeffs()] for f in res['facts']],
                'roots': roots, 'S': S, 'S_contains_all_large_t': unbounded,
                't_star': tstar, 'n_roots': len(roots)}
    out = {'beta': fs(beta), 'degrees': {k: len(v) - 1 for k, v in conds.items()},
           'legal_threshold_t': fs(F(beta.denominator, beta.numerator - 2 * beta.denominator)) if beta > 2 else ('never' if beta <= 2 else None),
           'full': pack(res), 'stab_pos_only': pack(res_sp)}
    return out, conds, extra

def in_S(t, packed):
    """membership of rational t in S (list of open intervals with root-interval endpoints). Returns True/False/None(undecided: t inside a root's isolating interval)."""
    t = F(t)
    for r in packed['roots']:
        if F(r['lo']) <= t <= F(r['hi']):
            return None
    for iv in packed['S']:
        lo_ok = iv['left'] is None or t > F(iv['left']['hi'])
        hi_ok = iv['right'] is None or t < F(iv['right']['lo'])
        if lo_ok and hi_ok:
            return True
    return False

def job(args):
    idx, rec, cwg = args
    t0 = time.time()
    out = {'idx': idx, 'd6': rec['d6'], 'edges': rec['edges'], 'problems': []}
    try:
        n, A = E.parse_digraph6(rec['d6'])
        assert n == 8
        if sum(map(sum, A)) != rec['edges']:
            out['problems'].append('edge count mismatch')
        out['complete'] = all(A[i][j] == 1 for i in range(n) for j in range(n) if i != j)
        if out['complete']:
            out['problems'].append('graph is complete (not a counterexample candidate)')
        cell = cwg['witness_in_first_stable_cell']
        beta = F(cell['beta'])
        assert F(cell['eps']) == F(1, int(cell['t']) + 1) and F(cell['delta']) == (beta - 1) * F(cell['eps'])
        out['beta_cell'] = fs(beta)
        out['W_definition_check'] = E.build_K_vs_W_check(A, beta.numerator, beta.denominator)
        if not out['W_definition_check']:
            out['problems'].append('K formula disagrees with W definition')
        # beta strictly inside the stored stable interval
        iv = rec['cert']['stable_intervals'][0]
        ok, info = E.beta_in_interval(beta, iv)
        out['beta_in_first_stable_interval'] = ok
        out['beta_interval_info'] = json.loads(json.dumps(info, default=str))
        if not ok:
            out['problems'].append('beta_cell not strictly inside first stable interval')
        # symbolic analysis at beta_cell
        an, conds, extra = analyse_beta(A, beta, 'cell')
        out['cell'] = an
        if not an['full']['S_contains_all_large_t']:
            out['problems'].append('FATAL: at beta_cell the all-conditions set does not contain all large t')
        # cell witness t must be in S
        tc = int(cell['t'])
        m = in_S(tc, an['full'])
        out['cell_t_in_S'] = m
        if m is not True:
            out['problems'].append('cell witness t=%d not certified inside S (membership=%s)' % (tc, m))
        pw = E.pointwise(A, F(cell['eps']), F(cell['delta']))
        out['cell_pointwise'] = {'stable': pw['stable'], 'positive': pw['positive'], 'detail': pw['detail'], 'q': fs(pw['q']), 'legal': pw['legal']}
        if not (pw['stable'] is True and pw['positive'] is True and pw['legal']):
            out['problems'].append('cell witness fails pointwise re-decision')
        if F(cell['q']) != pw['q']:
            out['problems'].append('cell witness q mismatch')
        # pointwise/symbolic sign agreement at beta_cell for the 4 t values (+ membership consistency)
        agree = []
        for t in TS:
            eps = F(1, t + 1); delta = (beta - 1) * eps
            p = E.pointwise(A, eps, delta)
            sy_H = [E.peval(extra['H'][k], F(t)) for k in range(n)]
            sy_D = E.peval(extra['D'], F(t)); sy_N = [E.peval(extra['N'][i], F(t)) for i in range(n)]
            hs = all(((a > 0) == (b > 0)) and ((a < 0) == (b < 0)) for a, b in zip(p['minors'], sy_H))
            xs = (p['x'] is not None and sy_D != 0 and all(((xi > 0) == (sy_N[i] * sy_D > 0)) for i, xi in enumerate(p['x']))) or (p['x'] is None and sy_D == 0)
            mem = in_S(t, an['full'])
            flag = (p['stable'] is True and p['positive'] is True and p['legal'])
            cons = (mem is None) or (mem == flag)
            agree.append({'t': t, 'hurwitz_signs_agree': hs, 'x_signs_agree': xs, 'in_S': mem, 'pointwise_ok': flag, 'consistent': cons,
                          'stable': p['stable'], 'positive': p['positive']})
            if not (hs and xs and cons):
                out['problems'].append('symbolic/pointwise disagreement at beta_cell t=%d' % t)
        out['cell_crosscheck'] = agree
        # stored witnesses (their own beta)
        wits = rec['witnesses']
        betas = set(w['beta'] for w in wits)
        out['beta_stored'] = sorted(betas)
        sw = []
        for w in wits:
            eps = F(w['eps']); delta = F(w['delta']); bw = F(w['beta'])
            assert eps == F(1, w['t'] + 1) and delta == (bw - 1) * eps
            p = E.pointwise(A, eps, delta)
            ok = (p['stable'] == w['stable']) and (p['positive'] == w['positive']) and (p['q'] == F(w['q']))
            sw.append({'t': w['t'], 'beta': w['beta'], 'stored_stable': w['stable'], 'stored_positive': w['positive'],
                       'my_stable': p['stable'], 'my_positive': p['positive'], 'my_detail': p['detail'], 'agree': ok})
            if not ok:
                out['problems'].append('stored witness flags disagree with exact re-decision at t=%d' % w['t'])
        out['stored_pointwise'] = sw
        # symbolic analysis at stored beta(s), compare flags with membership in the stability+positivity set
        out['stored'] = {}
        for bs in sorted(betas, key=F):
            bsF = F(bs)
            an2, conds2, extra2 = analyse_beta(A, bsF, 'stored')
            ok_iv, info2 = E.beta_in_interval(bsF, iv)
            chk = []
            for w in wits:
                if w['beta'] != bs: continue
                mem = in_S(w['t'], an2['stab_pos_only'])
                flag = bool(w['stable'] and w['positive'])
                cons = (mem is None) or (mem == flag)
                chk.append({'t': w['t'], 'in_S_stab_pos': mem, 'stored_stable_positive': flag, 'consistent': cons})
                if not cons:
                    out['problems'].append('stored witness flag inconsistent with exact t-set at beta=%s t=%d' % (bs, w['t']))
            # expectation: beta strictly inside the limit-stable interval <=> stab/pos set contains all large t
            unb = an2['stab_pos_only']['S_contains_all_large_t']
            # beta outside: strictly below lo or strictly above hi?
            lo_side = info2['lo'][-1]; hi_side = info2['hi'][-1]
            expect = None
            if ok_iv: expect = True
            elif (not lo_side) or (not hi_side):
                # strictly outside (compare functions give strict results unless beta == root exactly, which is impossible for irrational roots)
                expect = False
            out['stored'][bs] = {'analysis': an2, 'beta_in_first_stable_interval': ok_iv, 'interval_info': json.loads(json.dumps(info2, default=str)),
                                 'witness_checks': chk, 'stab_pos_set_contains_all_large_t': unb, 'expected_from_limit_interval': expect}
            if expect is not None and expect != unb:
                out['problems'].append('limit-interval membership of stored beta=%s (%s) contradicts exact t-set unboundedness (%s)' % (bs, expect, unb))
    except Exception as e:
        out['problems'].append('ERROR: ' + repr(e))
        out['traceback'] = traceback.format_exc()
    out['seconds'] = round(time.time() - t0, 3)
    return out

def main():
    outp = sys.argv[1]; wall = float(sys.argv[2]); workers = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    recs = [json.loads(l) for l in open(HITS)]
    cw = json.load(open(CW))['graphs']
    assert len(recs) == len(cw) and all(r['d6'] == g['d6'] for r, g in zip(recs, cw))
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.loads(l)['idx'])
            except Exception: pass
    nulls = [i for i, g in enumerate(cw) if g['witness_stored_legal'] is None]
    order = [i for i in range(len(recs)) if i % 4 == 0] + [i for i in nulls if i % 4 != 0]
    order += [i for i in range(len(recs)) if i % 4 == 2 and i not in nulls]
    order += [i for i in range(len(recs)) if i % 2 == 1 and i not in nulls]
    assert sorted(order) == list(range(len(recs)))
    todo = [i for i in order if i not in done]
    print('done', len(done), 'todo', len(todo), flush=True)
    t0 = time.time(); k = 0
    with Pool(workers) as pool, open(outp, 'a') as f:
        for out in pool.imap(job, ((i, recs[i], cw[i]) for i in todo), chunksize=2):
            f.write(json.dumps(out) + '\n'); f.flush(); k += 1
            if out['problems']:
                print('PROBLEM idx', out['idx'], out['problems'], flush=True)
            if k % 100 == 0:
                print('progress', k, 'in', round(time.time() - t0), 's', flush=True)
            if time.time() - t0 > wall:
                print('wall limit reached after', k, flush=True)
                pool.terminate(); break
    print('finished batch', k, 'elapsed', round(time.time() - t0), flush=True)

if __name__ == '__main__':
    main()
