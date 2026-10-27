#!/usr/bin/env python3
"""write_summary.py -- certificates/n8_limit_summary.json from n8/n8_accounting.json, n8/n8_prefilter_TOTAL.txt,
n8/n8_limit_tally.json and n8/n8_limit_HITS.jsonl (mirrors certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit_summary.json)."""
import json, collections, os, sys
from fractions import Fraction as F
acc = json.load(open('n8/n8_accounting.json')); tally = json.load(open('n8/n8_limit_tally.json'))
tot = dict(zip(*[iter(open('n8/n8_prefilter_TOTAL.txt').read().split())] * 2)); tot = {k: int(v) for k, v in tot.items()}
hits = [json.loads(l) for l in open('n8/n8_limit_HITS.jsonl')]
edge_counts = sorted(collections.Counter(h['edges'] for h in hits).items())
# stable beta-intervals of the candidates: lower end (beta_min, always 2 = the legality limit q -> 1) and upper end
def ends(h):
    out = []
    for iv in h['cert']['stable_intervals']:
        lo, hi = iv['lo'], iv['hi']
        out.append((lo.get('beta') or lo.get('beta_interval'), hi.get('beta') or hi.get('beta_interval'), hi.get('kind') or ('unbounded' if hi.get('beta') == 'inf' else 'algebraic' if 'beta_interval' in hi else 'rational')))
    return out
n_intervals = collections.Counter(len(h['cert']['stable_intervals']) for h in hits)
unbounded = sum(1 for h in hits if any(iv['hi'].get('beta') == 'inf' for iv in h['cert']['stable_intervals']))
start_at_2 = sum(1 for h in hits if str(h['cert']['stable_intervals'][0]['lo'].get('beta')) == '2')
wit_ok = sum(1 for h in hits if all(w['stable'] and w['positive'] for w in h['witnesses']))
wit_legal = sum(1 for h in hits if any(w['stable'] and w['positive'] and F(w['q']) > 1 for w in h['witnesses']))
summary = {
 'n': 8, 'date': '2026-10-28',
 'scope': 'eps -> 0 limit classification of the full support of every digraph on 8 vertices whose underlying graph has 15..28 edges (digraphs with sparser underlying graphs are excluded by the e_2 edge bound, Lemma A.2); finite-eps map not done',
 'enumeration': {'digraphs_read': tot['total'], 'exact_count': acc.get('grand_total'), 'domination_free': tot['domfree'], 'spectral_survivors': tot['survivors'], 'strict': tot['survivors'] - tot['boundary'], 'boundary_never_limit_stable': tot['boundary'], 'overflow_undecided': tot['overflow'],
                 'accounting': {k: acc[k] for k in acc if k in ('ok', 'kept_lines', 'jobs_checked', 'split_graphs_checked', 'shards_checked')},
                 'schemes': 'geng shards 0..8 (n8/parts); directg parts (n8/jobs); K8 by bidirected subgraph (n8/k8decomp, 40 batches, ledger check_k8decomp.py); K8-K3, K8-K_{1,5}, K8-(K2+K_{1,5}) by orbitenum (n8/orbitenum, 14 parts, ledger check_orbit.py)'},
 'limit_stage': {k: tally[k] for k in ('S', 'O', 'never_positivity', 'never_full', 'stable_somewhere', 'escalate', 'complete', 'stable_noncomplete', 'seconds', 'sharded', 'interrupted_single_run_records_cross_checked') if k in tally},
 'hits': len(hits), 'hits_file': 'n8/n8_limit_HITS.jsonl', 'results_file': 'n8/n8_limit_results.jsonl',
 'edge_counts': edge_counts,
 'stable_intervals_per_hit': sorted(n_intervals.items()), 'hits_with_unbounded_stable_interval': unbounded, 'hits_whose_first_interval_starts_at_beta_2': start_at_2,
 'hits_with_all_four_pointwise_witnesses_stable_and_positive': wit_ok, 'hits_with_a_legal_stable_positive_witness_q_gt_1': wit_legal,
 'status': 'CANDIDATES pending hand-verification by the user and independent machine confirmation (as for n = 7, task 28); a claimed counterexample stops the run and goes to the human before anything is built on it',
}
json.dump(summary, open('certificates/n8_limit_summary.json', 'w'), indent=1)
print(json.dumps({k: summary[k] for k in summary if k not in ('edge_counts',)}, indent=1)[:3000]); print('edge_counts', edge_counts)
