"""(a) resolve the 19 own/filter8 flag mismatches on blow-ups by exact root isolation of chi_perp(w - 1/2);
   (b) independent root-count cross-check of the Hurwitz decisions on the 60 S-flagged non-hit blow-ups (Jacobian charpoly);
   (c) for in-hits blow-ups whose stored interval differs from the parent's: exact decision at a beta that is inside one interval
       and outside the other, at t = 10^5 (consistency of the stored n = 8 interval with the blow-up, and Lemma W non-equality)."""
import sys, json, time
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/completeness_twins')
from lib import *
from rootcount import root_signs
from part2b_spectral import chi_perp_shifted, domination_free, classify_limit
from fractions import Fraction as F
from collections import Counter
R = '/home/user/CLTN/'; D = R + 'notes/verify_n8_2/completeness_twins/'
t0 = time.time()
rows = [json.loads(l) for l in open(D + 'part2_blowups.jsonl')]
flags = {json.loads(l)['blowup_d6']: json.loads(l) for l in open(D + 'part2_blowups_flags.jsonl')}
mism = json.load(open(D + 'part2b_summary.json'))['mismatches']
out = {'a_mismatch_resolution': [], 'b_rootcount_crosscheck': {}, 'c_interval_difference_checks': []}
for s in mism:
    B = d6_decode(s); c = chi_perp_shifted(B); neg, zero, pos = root_signs(c)
    fl = flags[s]; r = next(x for x in rows if x['blowup_d6'] == s)
    out['a_mismatch_resolution'].append({'blowup_d6': s, 'own_flag': fl['own_filter_flag'], 'filter8': (r['filter8_stdout'].split()[-1] if r['filter8_stdout'] else 'DROPPED'),
        'roots_of_chi_perp_by_Re_vs_minus_half': {'Re<-1/2': int(neg), 'Re=-1/2': int(zero), 'Re>-1/2': int(pos)}, 'alpha_vs_minus_half': ('<' if pos == 0 and zero == 0 else ('=' if pos == 0 else '>'))})
print('a done', round(time.time() - t0, 1), flush=True)
# (b)
cnt = Counter(); nb = 0
for r in rows:
    if r['in_hits_own'] or not r['filter8_stdout'] or not r['filter8_stdout'].endswith(' S'): continue
    B = d6_decode(r['blowup_d6'])
    for d in r['decisions']:
        if d['t'] != 100000 or d['which'] != 'mid': continue
        dec = decide(B, F(d['eps']), F(d['delta'])); neg, zero, pos = root_signs(dec['charpoly'])
        rc = 'stable' if (zero == 0 and pos == 0) else ('undecided' if zero > 0 else 'unstable')
        cnt[(dec['stability'], rc, dec['stability'] == rc)] += 1; nb += 1
out['b_rootcount_crosscheck'] = {'points': nb, 'hurwitz_vs_rootcount': {str(k): v for k, v in cnt.items()}}
print('b done', round(time.time() - t0, 1), flush=True)
# (c)
n7 = {json.loads(l)['d6']: json.loads(l) for l in open(R + 'certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')}
hits = {json.loads(l)['d6']: json.loads(l) for l in open(R + 'n8/n8_limit_HITS.jsonl')}
def bounds(iv, steps):
    lo = endpoint_bounds(iv['lo']); hi = endpoint_bounds(iv['hi'])
    for _ in range(steps):
        if 'beta' not in iv['lo']: lo = refine_root(iv['lo']['beta_minpoly_hi_to_lo'], lo[0], lo[1], 1)
        if 'beta' not in iv['hi']: hi = refine_root(iv['hi']['beta_minpoly_hi_to_lo'], hi[0], hi[1], 1)
    return lo, hi
ccount = Counter()
for r in rows:
    if not r['in_hits_own'] or r['interval_equal_to_parent']: continue
    P = n7[r['parent_d6']]; H = hits[r['hit_d6']]; A = d6_decode(P['d6']); B = d6_decode(r['blowup_d6'])
    piv = P['cert']['stable_intervals'][0]; hiv = H['cert']['stable_intervals'][0]
    plo, phi = bounds(piv, 30); hlo, hhi = bounds(hiv, 30)
    rec = {'blowup_d6': r['blowup_d6'], 'parent_d6': P['d6'], 'hit_d6': H['d6'], 'parent_interval_bounds': [str(plo[0]), str(phi[1])], 'hit_interval_bounds': [str(hlo[0]), str(hhi[1])], 'tests': []}
    # a beta strictly inside the parent interval and strictly outside the hit interval, or vice versa
    cands = []
    if phi[0] > hhi[1]:  # parent extends further right
        b = (hhi[1] + phi[0]) / 2
        if plo[1] < b: cands.append(('in_parent_out_hit', b))
    if hhi[0] > phi[1]:
        b = (phi[1] + hhi[0]) / 2
        if hlo[1] < b: cands.append(('in_hit_out_parent', b))
    if plo[0] > hlo[1] and hlo[1] < plo[0]:
        b = (hlo[1] + plo[0]) / 2
        if b < hhi[0]: cands.append(('in_hit_out_parent_left', b))
    if hlo[0] > plo[1]:
        b = (plo[1] + hlo[0]) / 2
        if b < phi[0]: cands.append(('in_parent_out_hit_left', b))
    for kind, b in cands:
        t = 100000; eps = F(1, t + 1); delta = (b - 1) * eps
        dB = decide(B, eps, delta); dA = decide(A, eps, delta)
        exp_b = kind.startswith('in_hit'); exp_a = kind.startswith('in_parent')
        okB = ((dB['stability'] == 'stable' and dB['positive']) == exp_b); okA = ((dA['stability'] == 'stable' and dA['positive']) == exp_a)
        rec['tests'].append({'kind': kind, 'beta': str(b), 't': t, 'blowup': [dB['stability'], dB['positive']], 'parent': [dA['stability'], dA['positive']], 'consistent_with_stored_intervals': okB and okA})
        ccount[(kind, okB and okA)] += 1
    if not cands: ccount[('no_separating_beta_found', None)] += 1
    out['c_interval_difference_checks'].append(rec)
out['c_summary'] = {str(k): v for k, v in ccount.items()}
json.dump(out, open(D + 'part5_crosschecks.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'c_interval_difference_checks'}, indent=1)); print(out['c_summary'], round(time.time() - t0, 1))
