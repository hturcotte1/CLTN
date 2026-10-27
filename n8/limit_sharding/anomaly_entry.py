#!/usr/bin/env python3
"""anomaly_entry.py [--apply] -- compose (and with --apply append) the ANOMALIES.md entry for the n = 8 limit stage."""
import json, sys, collections, os, time
from fractions import Fraction as F
tally = json.load(open('n8/n8_limit_tally.json')); summ = json.load(open('certificates/n8_limit_summary.json'))
hits = [json.loads(l) for l in open('n8/n8_limit_HITS.jsonl')]
tot = summ['enumeration']; ec = summ['edge_counts']
qhi = []  # upper ends of the stable q-intervals (q = beta - 1 along the limit ray): rational bounds of the algebraic endpoint
for h in hits:
    iv = h['cert']['stable_intervals'][-1]['hi']; qi = iv.get('q_interval'); 
    if qi: qhi.append(F(qi[1]))
qmax = max(qhi) if qhi else None
lo_root = sum(1 for h in hits if h['cert']['stable_intervals'][0]['lo'].get('kind') == 'root')
entry = f"""
## 2026-10-28  n = 8 LIMIT STAGE COMPLETE: {len(hits):,} candidate counterexamples (full support, eps -> 0) -- RUN STOPPED for hand-verification

Every digraph on 8 vertices whose underlying graph has 15..28 edges ({tot['digraphs_read']:,} of them, the exact
Burnside count; sparser underlying graphs are excluded by the e_2 edge bound, Lemma A.2) was read by filter8
exactly once (n8/n8_accounting.json ok; per-B ledgers of the K8 batches and of the three orbitenum graphs ok):
{tot['domination_free']:,} domination-free, {tot['spectral_survivors']:,} spectral survivors of which
{tot['boundary_never_limit_stable']:,} sit on the boundary alpha = -1/2 (never limit-stable, exact) and
{tot['strict']:,} are strict; overflow 0.  The exact limit stage (fastregion.region with replay, no floating
point) classified all {tally['S']:,} strict survivors: never_positivity {tally['never_positivity']:,},
never_full {tally['never_full']:,}, stable_somewhere {tally['stable_somewhere']:,} (K8 itself plus
{tally['stable_noncomplete']:,} non-complete graphs = the candidates), escalate {tally['escalate']}.
Candidates by arc count: {', '.join(f'{e}: {c}' for e, c in ec)}.  {summ['hits_with_all_four_pointwise_witnesses_stable_and_positive']:,} of the
{len(hits):,} have all four stored pointwise witnesses (t = 100..100000 along beta = midpoint of the stable
interval) stable with a positive fixed point; the t = 100 witness is often illegal (q <= 1, as at n = 7), so a hand
check must use a witness with q > 1 ({summ['hits_with_a_legal_stable_positive_witness_q_gt_1']:,} have one stored).  Stable
intervals: one per candidate, {lo_root:,} start at an algebraic beta > 2, the rest at beta_min = 2 (q -> 1); the
largest upper end of a stable q-interval is below {qmax} (rational bound of the algebraic endpoint).
Scale against n = 7: 123 candidates there, {len(hits):,} here; all strict survivors come from 17 high-symmetry
underlying graphs (K8: 166,042 strict lines; K8-e: 143,266; K8-2K2: 19,701; ...), none from the 135 packed jobs.

Pipeline notes.  The single-threaded limit stage was killed by a container restart at 59,610 of 347,639 lines
(11:2x UTC Oct 27, the session itself was down until 22:58); it was finished as four shards of the UNMODIFIED
n8/run_n8_limit.py over four disjoint symlinked subsets of the job files (sharded, {tally['seconds']:,.0f} CPU-s), merged by
a checked script (disjoint d6 sets whose union is exactly the S/O line set; tally recomputed from the merged
records equals the sum of the shard tallies; every one of the {tally['interrupted_single_run_records_cross_checked']:,} complete records of the
interrupted single run equals the shard record for the same graph).  after_prefilter2.sh now refuses to run a
second time (a second pass would rewrite the tally as all-skipped).  Post-hoc content check of the orbit
enumerator on K8-e (the richest graph): 16 orbit parts re-enumerated the whole graph; parts 0..13
(4,666,650,090 digraphs, 125,256 strict lines) are all contained in the 143,266 canonical strict lines of the
321 directg parts; parts 14 and 15 were cut by the restart and are being redone for the full set equality.

STATUS: candidates, not results.  Per CLAUDE.md a claimed counterexample stops the run: nothing is built on
these until the user has hand-verified a sample and the independent machine confirmation (the four n = 7
verifiers of notes/verify_124, re-run on the n = 8 list) has passed.  Files: n8/n8_limit_HITS.jsonl (one record
per candidate with the exact stable interval, critical minimal polynomials and four pointwise witnesses),
n8/n8_limit_results.jsonl (every strict survivor), n8/n8_limit_tally.json, certificates/n8_limit_summary.json.
"""
print(entry)
if '--apply' in sys.argv:
    open('ANOMALIES.md', 'a').write(entry); print('appended')
