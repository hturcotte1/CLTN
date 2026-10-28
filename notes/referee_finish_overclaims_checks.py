"""Exact checks backing notes/referee_finish_overclaims.json (overclaim referee, 2026-10-28).

Run from anywhere:  nice -n 10 python3 /home/user/CLTN/notes/referee_finish_overclaims_checks.py
Writes notes/referee_finish_overclaims_checks.json.  Exact arithmetic only (Fraction, integer polynomials,
three-valued Hurwitz test of ctln.hurwitz_stable; nauty labelg via notes/n8_structure_tools/lib.py).

A. One-vertex deletions of the 9,133 equal to the finite-eps-only n = 7 counterexample &FNFL~[udo?
   (the structure tools compare deletions with the 123 limit counterexamples only).
B. The 16 known eight-node finite-eps-only counterexamples (E_1 of Theorem E, the 15 Lemma B 'cube' graphs):
   witnesses re-decided with ctln.analyze_support, pairwise non-isomorphic, absent from the limit list.
C. Theorem Q's k -> infinity limit versus the n = 8 record (where the minimum real parts lie).
D. Plastic number versus the decimal 1.3247; champion's lower q end versus rho.
"""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN/notes/n8_structure_tools')
sys.path.insert(0, '/home/user/CLTN')
import lib, small_lists as sl
import ctln

R = '/home/user/CLTN/'
out = {}

# ---------------- A
hits = [json.loads(l) for l in open(R + 'n8/n8_limit_HITS.jsonl')]
d6s = [h['d6'] for h in hits]
FNFL = '&FNFL~[udo?'
fc = lib.labelg_canon([FNFL])[0]
n7set = set(lib.labelg_canon(sl.n7_list()))
n6set = set(lib.labelg_canon([d for _, d in sl.n6_list()]))
dels, idx = [], []
for i, d in enumerate(d6s):
    A = lib.from_d6(d)
    for v in range(8):
        dels.append(lib.to_d6(lib.induced(A, [x for x in range(8) if x != v])))
        idx.append((i, v))
can = lib.labelg_canon(dels)
fn = {}
for (i, v), c in zip(idx, can):
    if c == fc:
        fn.setdefault(i, []).append(v)
recs = []
for i in sorted(fn):
    A = lib.from_d6(d6s[i])
    sinks = [v for v in range(8) if all(A[v][j] == 1 for j in range(8) if j != v)]
    sink_dels = lib.labelg_canon([lib.to_d6(lib.induced(A, [x for x in range(8) if x != v])) for v in sinks])
    closed = lib.twin_pairs(A)['closed']
    six = []
    for a in range(8):
        for b in range(a + 1, 8):
            six.append(lib.to_d6(lib.induced(A, [x for x in range(8) if x not in (a, b)])))
    recs.append({'d6': d6s[i], 'FNFL_deletion_vertices': fn[i], 'sink_nodes': sinks,
                 'FNFL_vertex_is_sink': [v in sinks for v in fn[i]],
                 'FNFL_vertex_in_closed_twin_pair': [any(v in p for p in closed) for v in fn[i]],
                 'some_sink_deletion_in_123_limit_list': any(c in n7set for c in sink_dels),
                 'contains_an_n6_counterexample': any(c in n6set for c in lib.labelg_canon(six))})
out['A_FNFL_deletions'] = {
    'graphs_with_a_one_vertex_deletion_equal_to_FNFL': len(recs),
    'records': recs,
    'of_these_counted_in_the_1972_no_sink_deletion_gives_an_n7_cex': sum(1 for r in recs if not r['some_sink_deletion_in_123_limit_list']),
    'of_these_with_FNFL_at_a_closed_twin_vertex': sum(1 for r in recs if any(r['FNFL_vertex_in_closed_twin_pair'])),
}

# ---------------- B
E = json.load(open(R + 'notes/research_wf3/theory/family_E_witnesses.json'))
C = json.load(open(R + 'notes/research_wf3/theory/cube_n8_witnesses.json'))
cand = [('E_1', r['d6'], r['witness'], r.get('A_rows')) for r in E if r['n'] == 8]
cand += [('cube', r['d6'], r['witness'], None) for r in C]
hitcan = set(lib.labelg_canon(d6s))
cc = lib.labelg_canon([d for _, d, _, _ in cand])
brec = []
for (kind, d, w, rows), c in zip(cand, cc):
    A = lib.from_d6(d)
    if rows is not None:
        assert [''.join(str(x) for x in r) for r in A] == rows, d
    eps, delta = F(w['eps']), F(w['delta'])
    rec = ctln.analyze_support(A, list(range(8)), eps, delta)
    brec.append({'kind': kind, 'd6': d, 'eps': str(eps), 'delta': str(delta), 'q': str(ctln.q_of(eps, delta)),
                 'legal': ctln.legal(eps, delta), 'positive': rec['positive'], 'stable': rec['stable'],
                 'detail': rec['stability_detail'], 'clique': rec['clique'], 'in_limit_list': c in hitcan})
out['B_finite_eps_only_n8'] = {
    'count': len(brec), 'pairwise_non_isomorphic': len(set(cc)) == len(cc),
    'distinct_isomorphism_classes': len(set(cc)),
    'isomorphic_groups': [[d for (_, d, _, _), c2 in zip(cand, cc) if c2 == c] for c in sorted(set(cc)) if cc.count(c) > 1],
    'all_legal_positive_stable_nonclique': all(r['legal'] and r['positive'] and r['stable'] is True and not r['clique'] for r in brec),
    'none_in_limit_list': not any(r['in_limit_list'] for r in brec), 'records': brec}

# ---------------- C
def shifted_left_test(coeffs_hi_to_lo, c):
    """True/False/None: do all roots z of p satisfy Re z > -c ?  (Hurwitz test of r(w) = p(-w - c))."""
    m = len(coeffs_hi_to_lo) - 1
    # expand p(-w - c) = sum a_k (-w - c)^(m-k)
    res = [F(0)] * (m + 1)  # lo..hi in w
    for k, a in enumerate(coeffs_hi_to_lo):
        e = m - k
        # (-w - c)^e = (-1)^e (w + c)^e
        from math import comb
        for j in range(e + 1):
            res[j] += F(a) * (-1) ** e * comb(e, j) * F(c) ** (e - j)
    hi_to_lo = list(reversed(res))
    if hi_to_lo[0] < 0:
        hi_to_lo = [-x for x in hi_to_lo]
    v, det = ctln.hurwitz_stable(hi_to_lo)
    return v, det
qlim = [1, 1, 3, 1, 1]                  # z^4+z^3+3z^2+z+1, notes_theory.md line 357 (Theorem Q limit)
rec8 = [8, 16, 34, 33, 30, 13, 5]       # psi_D / z of the extremal pair, notes/n8_structure.md 1.8
out['C_theorem_Q_vs_record'] = {
    'Q_limit_poly_all_roots_Re_gt_minus_5_over_14': shifted_left_test(qlim, F(5, 14)),
    'record_sextic_all_roots_Re_gt_minus_4_over_11': shifted_left_test(rec8, F(4, 11)),
    'record_sextic_all_roots_Re_gt_minus_3_over_8': shifted_left_test(rec8, F(3, 8)),
    'reading': 'If the first is True and the second False: record min Re < -4/11 < -5/14 < Q-limit min Re, so the record '
               'lies to the LEFT of (and closer to -1/2 than) the Theorem Q limit, and its q is smaller.',
    'q8_limit_upper_isolating_hi': '2755346911/1610612736',
    'q8_limit_hi_lt_7_4': F(2755346911, 1610612736) < F(7, 4),
    'theorem_Q_limit_stable_q_values_for_large_k': ['3/2 (k>=25)', '17/10 (k>=79)', '9/5 (k>=286)'],
    'nine_fifths_gt_7_4': F(9, 5) > F(7, 4),
}

# ---------------- D
f = lambda q: q ** 3 - q - 1
champ = [h for h in hits if h['d6'] == '&GMx|\\Aw^qz]w'][0]
lo = champ['cert']['stable_intervals'][0]['lo']
out['D_plastic_and_champion'] = {
    'f(13247/10000)': str(f(F(13247, 10000))), 'rho_gt_1.3247': f(F(13247, 10000)) < 0,
    'f(49/37)': str(f(F(49, 37))), 'rho_gt_49/37': f(F(49, 37)) < 0,
    'champion_lower_end_q_interval': lo['q_interval'], 'champion_lower_end_q_minpoly': lo['q_minpoly_hi_to_lo'],
    'champion_lower_q_hi_lt_49/37': F(lo['q_interval'][1]) < F(49, 37),
}
json.dump(out, open(R + 'notes/referee_finish_overclaims_checks.json', 'w'), indent=1)
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'records'} for k, v in out.items()}, indent=1))
