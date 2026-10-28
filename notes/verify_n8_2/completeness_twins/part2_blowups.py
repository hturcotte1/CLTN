"""Part 2: closed-twin blow-ups of the 123 n=7 limit counterexamples vs the 9,133 n=8 candidates."""
import sys, json, time, subprocess
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/completeness_twins')
from lib import *
from fractions import Fraction as F
R = '/home/user/CLTN/'
D = R + 'notes/verify_n8_2/completeness_twins/'
LABELG = R + 'n7/labelg'; FILTER8 = R + 'n7/filter8'

def labelg(lines):
    out = subprocess.run([LABELG, '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    res = out.stdout.split('\n')
    res = [x for x in res if x]
    assert len(res) == len(lines), (len(res), len(lines), out.stderr[:500])
    return res

t0 = time.time()
hits = [json.loads(l) for l in open(R + 'n8/n8_limit_HITS.jsonl')]
n7 = [json.loads(l) for l in open(R + 'certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
# Geneson positive control
G = json.load(open(R + 'certificates/stage0_geneson_2607.21396.json'))['A_ours']
dg = decide(G, F(1, 101), F(29, 2525)); print('Geneson control eps=1/101 delta=29/2525:', dg['stability'], dg['positive'])
dg = decide(G, F(1, 10), F(29, 250)); print('Geneson control eps=1/10 (paper says unstable, minor 3 negative):', dg['stability'], dg['positive'], [m < 0 for m in dg['minors']])

# canonical forms of the hits, both ways
hit_own = {}; hit_lab = {}
for r in hits:
    hit_own[canon(d6_decode(r['d6']))] = r
labs = labelg([r['d6'] for r in hits])
for r, l in zip(hits, labs): hit_lab[l] = r
print('hits: own canon distinct', len(hit_own), 'labelg distinct', len(hit_lab), round(time.time() - t0, 1), 's')
# results file labelg index
res_lab = {}
for line, lab in zip(open(R + 'n8/n8_limit_results.jsonl'), open(D + 'results_labelg.txt')):
    r = json.loads(line); res_lab[lab.strip()] = r
assert len(res_lab) == 347639
# the HITS labelg canons must all be in results with verdict stable_somewhere
assert all(l in res_lab and res_lab[l]['verdict'] == 'stable_somewhere' for l in hit_lab)

# blow-ups
blow = []  # (parent idx, v, B, d6)
for pi, r in enumerate(n7):
    A = d6_decode(r['d6']); assert not is_complete(A)
    for v in range(7):
        B = blow_up(A, v); blow.append((pi, v, B, d6_encode(B)))
own = [canon(B) for (_, _, B, _) in blow]
lab = labelg([s for (_, _, _, s) in blow])
print('blow-ups built:', len(blow), 'distinct own', len(set(own)), 'distinct labelg', len(set(lab)))
# own and labelg must induce the same partition
from collections import defaultdict
pa = defaultdict(set); pb = defaultdict(set)
for i, (a, b) in enumerate(zip(own, lab)): pa[a].add(i); pb[b].add(i)
same_partition = sorted(map(sorted, pa.values())) == sorted(map(sorted, pb.values()))
print('partition agreement own/labelg:', same_partition)
assert same_partition

rows = []; seen = set(); n_in = 0; n_out = 0
cmp_rows = []
for i, (pi, v, B, s) in enumerate(blow):
    if own[i] in seen: continue
    seen.add(own[i])
    in_own = own[i] in hit_own; in_lab = lab[i] in hit_lab
    assert in_own == in_lab, (s, in_own, in_lab)
    rec = {'parent_d6': n7[pi]['d6'], 'parent_idx': pi, 'v': v, 'blowup_d6': s, 'own_canon': own[i], 'labelg': lab[i],
           'arcs': arcs(B), 'underlying_edges': underlying_edges(B), 'in_hits_own': in_own, 'in_hits_labelg': in_lab}
    if in_own:
        n_in += 1
        h = hit_lab[lab[i]]
        assert canon(d6_decode(h['d6'])) == own[i]
        p = find_isomorphism(B, d6_decode(h['d6'])); assert p is not None
        rec['hit_d6'] = h['d6']
        # interval comparison parent vs stored n=8 interval
        piv = n7[pi]['cert']['stable_intervals']; hiv = h['cert']['stable_intervals']
        rec['parent_n_intervals'] = len(piv); rec['hit_n_intervals'] = len(hiv)
        def ep_sig(ep):
            if 'beta' in ep and ep['kind'] == 'beta_min': return ('beta_min', str(F(ep['beta'])))
            return ('root', tuple(ep['beta_minpoly_hi_to_lo']))
        def same_root(e1, e2):
            if ep_sig(e1) != ep_sig(e2): return False
            if 'beta' in e1 and e1['kind'] == 'beta_min': return True
            if 'beta' in e1 or 'beta' in e2:
                b1 = endpoint_bounds(e1); b2 = endpoint_bounds(e2)
            else:
                b1 = endpoint_bounds(e1); b2 = endpoint_bounds(e2)
            for _ in range(60):
                if b1[1] < b2[0] or b2[1] < b1[0]: return False
                if b1[0] == b1[1] and b2[0] == b2[1]: return b1[0] == b2[0]
                if b1[1] - b1[0] < F(1, 2**50) and b2[1] - b2[0] < F(1, 2**50): return True
                if 'beta' not in e1: b1 = refine_root(e1['beta_minpoly_hi_to_lo'], b1[0], b1[1], 1)
                if 'beta' not in e2: b2 = refine_root(e2['beta_minpoly_hi_to_lo'], b2[0], b2[1], 1)
            return True
        eq = len(piv) == len(hiv) and all(same_root(a['lo'], b['lo']) and same_root(a['hi'], b['hi']) for a, b in zip(piv, hiv))
        rec['interval_equal_to_parent'] = eq
        rec['parent_interval'] = [[ep_sig(iv['lo']), iv['lo'].get('beta_interval', iv['lo'].get('beta'))], [ep_sig(iv['hi']), iv['hi'].get('beta_interval', iv['hi'].get('beta'))]] if False else None
        rec['parent_hi_minpoly_deg'] = len(piv[0]['hi']['beta_minpoly_hi_to_lo']) - 1 if 'beta_minpoly_hi_to_lo' in piv[0]['hi'] else 0
        rec['hit_hi_minpoly_deg'] = len(hiv[0]['hi']['beta_minpoly_hi_to_lo']) - 1 if 'beta_minpoly_hi_to_lo' in hiv[0]['hi'] else 0
        rec['parent_hi_beta_interval'] = piv[0]['hi'].get('beta_interval', [piv[0]['hi'].get('beta')] * 2)
        rec['hit_hi_beta_interval'] = hiv[0]['hi'].get('beta_interval', [hiv[0]['hi'].get('beta')] * 2)
        rec['parent_lo'] = piv[0]['lo'].get('beta', piv[0]['lo'].get('beta_interval')); rec['hit_lo'] = hiv[0]['lo'].get('beta', hiv[0]['lo'].get('beta_interval'))
    else:
        n_out += 1
        rr = res_lab.get(lab[i])
        rec['in_results_file'] = rr is not None
        if rr is not None:
            p = find_isomorphism(B, d6_decode(rr['d6'])); assert p is not None, 'labelg match not confirmed by own isomorphism search'
            rec['results_d6'] = rr['d6']; rec['results_verdict'] = rr['verdict']; rec['results_reason'] = rr.get('reason', 'never_full' if rr['verdict'] == 'never' else None)
            rec['results_flag'] = rr['flag']
        rec['in_enumerated_class'] = rec['underlying_edges'] >= 15 and rec['arcs'] >= 29
        # filter8 verdict
        f8 = subprocess.run([FILTER8], input=s + '\n', capture_output=True, text=True)
        rec['filter8_stdout'] = f8.stdout.strip(); rec['filter8_stderr'] = f8.stderr.strip()
        # exact decisions at the parent's interior beta
        piv = n7[pi]['cert']['stable_intervals']
        dec = []
        for k, iv in enumerate(piv):
            beta, lob, hib, steps = interior_point(iv, extra_steps=2)
            # also quartile points
            betas = [beta, lob[1] + (hib[0] - lob[1]) / 4, lob[1] + 3 * (hib[0] - lob[1]) / 4]
            for bj, b in enumerate(betas):
                for t in (10000, 100000):
                    eps = F(1, t + 1); delta = (b - 1) * eps; q = delta * (1 - eps) / eps
                    assert legal(eps, delta), (str(eps), str(delta))
                    dB = decide(B, eps, delta); dA = decide(d6_decode(n7[pi]['d6']), eps, delta)
                    dec.append({'interval': k, 'which': ['mid', 'q1', 'q3'][bj], 'beta': str(b), 't': t, 'eps': str(eps), 'delta': str(delta), 'q': str(q),
                                'lo_bounds': [str(x) for x in lob], 'hi_bounds': [str(x) for x in hib], 'bisection_steps': steps + 2,
                                'blowup_stability': dB['stability'], 'blowup_positive': dB['positive'],
                                'blowup_first_nonpositive_minor': next((j + 1 for j, m in enumerate(dB['minors']) if m <= 0), None),
                                'blowup_min_x': str(min(dB['x'])) if dB['x'] else None,
                                'parent_stability': dA['stability'], 'parent_positive': dA['positive']})
        rec['decisions'] = dec
        rec['blowup_counterexample_at_some_tested_point'] = any(d['blowup_stability'] == 'stable' and d['blowup_positive'] for d in dec)
        rec['parent_counterexample_at_all_tested_points'] = all(d['parent_stability'] == 'stable' and d['parent_positive'] for d in dec)
        rec['any_undecided'] = any(d['blowup_stability'] == 'undecided' or d['parent_stability'] == 'undecided' for d in dec)
    rows.append(rec)
    if len(rows) % 100 == 0: print(len(rows), 'distinct blow-ups processed', round(time.time() - t0, 1), 's', flush=True)

with open(D + 'part2_blowups.jsonl', 'w') as f:
    for r in rows: f.write(json.dumps(r) + '\n')
summ = {'blowups_built': len(blow), 'distinct_blowups': len(rows), 'in_hits': n_in, 'not_in_hits': n_out,
        'not_in_hits_in_results': sum(1 for r in rows if not r['in_hits_own'] and r['in_results_file']),
        'not_in_hits_verdicts': {}, 'not_in_hits_outside_class': sum(1 for r in rows if not r['in_hits_own'] and not r['in_enumerated_class']),
        'not_in_hits_filter8_flags': {}, 'missed_counterexamples': [r['blowup_d6'] for r in rows if not r['in_hits_own'] and r['blowup_counterexample_at_some_tested_point']],
        'parent_control_failures': [r['blowup_d6'] for r in rows if not r['in_hits_own'] and not r['parent_counterexample_at_all_tested_points']],
        'undecided_points': sum(1 for r in rows if not r['in_hits_own'] and r['any_undecided']),
        'in_hits_interval_equal_to_parent': sum(1 for r in rows if r['in_hits_own'] and r['interval_equal_to_parent']),
        'in_hits_interval_differs': sum(1 for r in rows if r['in_hits_own'] and not r['interval_equal_to_parent']),
        'min_arcs_blowup': min(r['arcs'] for r in rows), 'min_underlying_edges_blowup': min(r['underlying_edges'] for r in rows),
        'seconds': round(time.time() - t0, 1)}
from collections import Counter
summ['not_in_hits_verdicts'] = dict(Counter((r.get('results_verdict'), r.get('results_reason')) and str((r.get('results_verdict'), r.get('results_reason'))) for r in rows if not r['in_hits_own']))
summ['not_in_hits_filter8_flags'] = dict(Counter((r['filter8_stdout'].split()[-1] if r['filter8_stdout'] else 'DROPPED') for r in rows if not r['in_hits_own']))
json.dump(summ, open(D + 'part2_summary.json', 'w'), indent=1)
print(json.dumps(summ, indent=1))
