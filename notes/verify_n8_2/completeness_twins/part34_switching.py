"""Part 3: Lemma J closure (bidirected in-twin switching) inside the 9,133 n=8 candidates.
   Part 4: the same closure on the deterministic sample of 'never' records (line index % 1692 == 0, verdict never)."""
import sys, json, time, subprocess, itertools
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/completeness_twins')
from lib import *
from fractions import Fraction as F
from collections import Counter
R = '/home/user/CLTN/'; D = R + 'notes/verify_n8_2/completeness_twins/'
LABELG = R + 'n7/labelg'; FILTER8 = R + 'n7/filter8'

def labelg(lines):
    if not lines: return []
    out = subprocess.run([LABELG, '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    res = [x for x in out.stdout.split('\n') if x]
    assert len(res) == len(lines), (len(res), len(lines), out.stderr[:500]); return res

t0 = time.time()
hits = [json.loads(l) for l in open(R + 'n8/n8_limit_HITS.jsonl')]
res_lab = {}; res_lines = []
for line, lab in zip(open(R + 'n8/n8_limit_results.jsonl'), open(D + 'results_labelg.txt')):
    r = json.loads(line); res_lab[lab.strip()] = r; res_lines.append(r)
hit_own = {canon(d6_decode(r['d6'])): r for r in hits}
hit_lab = dict(zip(labelg([r['d6'] for r in hits]), hits))
assert len(hit_own) == len(hit_lab) == 9133

def interval_sig(r):
    out = []
    for iv in r['cert']['stable_intervals']:
        s = []
        for ep in (iv['lo'], iv['hi']):
            if ep['kind'] == 'beta_min': s.append(('beta_min', str(F(ep['beta']))))
            else: s.append(('root', tuple(ep['beta_minpoly_hi_to_lo']), tuple(str(F(x)) for x in ep['beta_interval'])))
        out.append(tuple(s))
    return tuple(out)

def same_intervals(r1, r2):
    """same number of intervals, same endpoint minpolys, and isolating intervals isolating the same root (refined)."""
    i1 = r1['cert']['stable_intervals']; i2 = r2['cert']['stable_intervals']
    if len(i1) != len(i2): return False
    for a, b in zip(i1, i2):
        for e1, e2 in ((a['lo'], b['lo']), (a['hi'], b['hi'])):
            if e1['kind'] != e2['kind']: return False
            if e1['kind'] == 'beta_min':
                if F(e1['beta']) != F(e2['beta']): return False
                continue
            if e1['beta_minpoly_hi_to_lo'] != e2['beta_minpoly_hi_to_lo']: return False
            b1 = endpoint_bounds(e1); b2 = endpoint_bounds(e2)
            ok = None
            for _ in range(80):
                if b1[1] < b2[0] or b2[1] < b1[0]: ok = False; break
                if (b1[1] - b1[0] < F(1, 2**60)) and (b2[1] - b2[0] < F(1, 2**60)): ok = True; break
                b1 = refine_root(e1['beta_minpoly_hi_to_lo'], b1[0], b1[1], 1); b2 = refine_root(e2['beta_minpoly_hi_to_lo'], b2[0], b2[1], 1)
            if ok is None: ok = True
            if not ok: return False
    return True

def all_switches(A):
    """every graph obtained by switching a nonempty subset of the switchable third vertices of one bidirected in-twin pair."""
    out = []
    for (u, v) in bidirected_in_twins(A):
        ws = [w for w in range(len(A)) if w not in (u, v) and A[u][w] + A[v][w] == 1]
        for k in range(1, len(ws) + 1):
            for sub in itertools.combinations(ws, k):
                B = [row[:] for row in A]
                for w in sub: B[u][w], B[v][w] = A[v][w], A[u][w]
                out.append(((u, v), sub, B))
    return out

def closure_check(records, name, expect_in_hits):
    rows = []; stats = Counter(); exceptions = []
    allq = []  # (rec index, pair, sub, B, d6)
    for ri, r in enumerate(records):
        A = d6_decode(r['d6'])
        tw = bidirected_in_twins(A); stats['records'] += 1; stats['records_with_twins'] += bool(tw); stats['twin_pairs'] += len(tw)
        sw = all_switches(A); stats['switches'] += len(sw); stats['single_switches'] += sum(1 for s in sw if len(s[1]) == 1)
        for (pair, sub, B) in sw: allq.append((ri, pair, sub, B, d6_encode(B)))
    labs = labelg([q[4] for q in allq])
    for (ri, pair, sub, B, s), lab in zip(allq, labs):
        r = records[ri]; c = canon(B)
        row = {'source_d6': r['d6'], 'twins': pair, 'switched_at': sub, 'switch_d6': s, 'own_canon': c, 'labelg': lab,
               'arcs': arcs(B), 'underlying_edges': underlying_edges(B), 'in_class': underlying_edges(B) >= 15 and arcs(B) >= 29}
        assert row['arcs'] == arcs(d6_decode(r['d6'])) and row['underlying_edges'] == underlying_edges(d6_decode(r['d6']))
        in_own = c in hit_own; in_lab = lab in hit_lab; assert in_own == in_lab
        row['in_hits'] = in_own
        rr = res_lab.get(lab); row['in_results'] = rr is not None
        if rr is not None:
            row['results_d6'] = rr['d6']; row['results_verdict'] = rr['verdict']; row['results_reason'] = rr.get('reason', 'never_full' if rr['verdict'] == 'never' else None)
            if expect_in_hits:
                assert find_isomorphism(B, d6_decode(rr['d6'])) is not None
                row['same_interval_as_source'] = same_intervals(r, rr) if rr['verdict'] == 'stable_somewhere' else None
                row['same_D_data_consistent'] = row['same_interval_as_source']
            else:
                row['same_verdict'] = rr['verdict'] == r['verdict']
                if rr['verdict'] == r['verdict'] == 'never':
                    row['same_reason'] = rr.get('reason') == r.get('reason')
        ok = (in_own and row.get('same_interval_as_source', False)) if expect_in_hits else (rr is not None and rr['verdict'] == 'never')
        row['ok'] = ok
        if not ok:
            # exception: decide exactly
            f8 = subprocess.run([FILTER8], input=s + '\n', capture_output=True, text=True)
            row['filter8'] = f8.stdout.strip() or 'DROPPED'; row['filter8_stderr'] = f8.stderr.strip()
            dec = []
            if expect_in_hits:
                for k, iv in enumerate(r['cert']['stable_intervals']):
                    beta, lob, hib, steps = interior_point(iv, extra_steps=2)
                    for t in (10000, 100000):
                        eps = F(1, t + 1); delta = (beta - 1) * eps
                        dB = decide(B, eps, delta); dA = decide(d6_decode(r['d6']), eps, delta)
                        dec.append({'beta': str(beta), 't': t, 'eps': str(eps), 'delta': str(delta), 'q': str(delta * (1 - eps) / eps), 'switch_stability': dB['stability'], 'switch_positive': dB['positive'], 'source_stability': dA['stability'], 'source_positive': dA['positive'], 'same_charpoly': dB['charpoly'] == dA['charpoly'], 'same_x': dB['x'] == dA['x']})
            else:
                for t in (10000, 100000):
                    for beta in (F(2001, 1000), F(21, 10), F(5, 2), F(4)):
                        eps = F(1, t + 1); delta = (beta - 1) * eps
                        dB = decide(B, eps, delta); dA = decide(d6_decode(r['d6']), eps, delta)
                        dec.append({'beta': str(beta), 't': t, 'switch_stability': dB['stability'], 'switch_positive': dB['positive'], 'source_stability': dA['stability'], 'source_positive': dA['positive'], 'same_charpoly': dB['charpoly'] == dA['charpoly'], 'same_x': dB['x'] == dA['x']})
            row['decisions'] = dec
            row['missed_counterexample'] = any(d['switch_stability'] == 'stable' and d['switch_positive'] for d in dec)
            exceptions.append(row)
        rows.append(row)
    stats['switch_results'] = len(rows); stats['ok'] = sum(1 for x in rows if x['ok']); stats['exceptions'] = len(exceptions)
    stats['distinct_switch_canons'] = len(set(x['own_canon'] for x in rows))
    stats['switch_in_hits'] = sum(1 for x in rows if x['in_hits']); stats['switch_in_results'] = sum(1 for x in rows if x['in_results'])
    stats['leaves_class'] = sum(1 for x in rows if not x['in_class'])
    if expect_in_hits:
        stats['same_interval'] = sum(1 for x in rows if x.get('same_interval_as_source'))
        stats['switch_equals_source_canon'] = sum(1 for x in rows if x['own_canon'] == canon(d6_decode(x['source_d6'])))
    else:
        stats['same_verdict'] = sum(1 for x in rows if x.get('same_verdict')); stats['same_reason'] = sum(1 for x in rows if x.get('same_reason'))
        stats['verdict_pairs'] = dict(Counter(str((records[0]['verdict'], x.get('results_verdict'), x.get('results_reason'))) for x in rows))
    stats['missed_counterexamples'] = sum(1 for x in exceptions if x.get('missed_counterexample'))
    with open(D + name + '_rows.jsonl', 'w') as f:
        for x in rows: f.write(json.dumps(x) + '\n')
    with open(D + name + '_exceptions.jsonl', 'w') as f:
        for x in exceptions: f.write(json.dumps(x) + '\n')
    stats['seconds'] = round(time.time() - t0, 1)
    print(name, json.dumps(dict(stats), indent=1), flush=True)
    return dict(stats)

s3 = closure_check(hits, 'part3_hits_closure', True)
# Part 4 sample: every 1,692nd line (0-based line index multiple of 1692) with verdict never
sample = [r for i, r in enumerate(res_lines) if i % 1692 == 0 and r['verdict'] == 'never']
print('part4 sample size', len(sample), 'of', len([r for i, r in enumerate(res_lines) if i % 1692 == 0]), 'sampled lines')
s4 = closure_check(sample, 'part4_never_closure', False)
s4['sample_rule'] = "0-based line index i of n8/n8_limit_results.jsonl with i % 1692 == 0 and verdict == 'never'"; s4['sample_size'] = len(sample)
json.dump({'part3': s3, 'part4': s4}, open(D + 'part34_summary.json', 'w'), indent=1)
