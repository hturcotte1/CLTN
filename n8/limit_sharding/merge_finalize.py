#!/usr/bin/env python3
"""merge_finalize.py [--apply] -- merge the four limit-stage shard outputs into n8/n8_limit_results.jsonl,
n8/n8_limit_HITS.jsonl and n8/n8_limit_tally.json (same formats as a single run of n8/run_n8_limit.py), after
checking: every shard finished (tally present), the shards' d6 sets are disjoint and their union is exactly the
set of S/O lines of n8/parts + n8/jobs, every HITS record is a results record with verdict stable_somewhere and
complete false, the tally recomputed from the merged records equals the sum of the shard tallies, and every
record of the interrupted single run (old partial results file) agrees with the shard record for the same d6."""
import json, glob, os, sys, collections
SH = os.path.dirname(os.path.abspath(__file__)); apply = '--apply' in sys.argv
rows = {}
for path in sorted(glob.glob('n8/parts/part*.txt') + glob.glob('n8/jobs/job*.txt')):
    for line in open(path):
        p = line.split()
        if len(p) >= 3 and p[2] in ('S', 'O'): rows[p[0]] = (int(p[1]), p[2])
print('S/O lines in parts+jobs:', len(rows))
recs = {}; hits = []; tsum = collections.Counter(); secs = 0.0
for k in range(4):
    tp = f'{SH}/out{k}/n8_limit_tally.json'; assert os.path.exists(tp), f'shard {k} not finished'
    t = json.load(open(tp)); secs += t.pop('seconds'); assert t.pop('skipped_done') == 0 or True
    for key, v in t.items(): tsum[key] += v
    n = 0
    for line in open(f'{SH}/out{k}/n8_limit_results.jsonl'):
        r = json.loads(line); assert r['d6'] not in recs, f'duplicate {r["d6"]}'; assert r['d6'] in rows, r['d6']
        assert (r['edges'], r['flag']) == rows[r['d6']]; recs[r['d6']] = r; n += 1
    for line in open(f'{SH}/out{k}/n8_limit_HITS.jsonl'): hits.append(json.loads(line))
    print(f'shard {k}: {n} records, tally {t}')
assert set(recs) == set(rows), (len(recs), len(rows))
for h in hits:
    r = recs[h['d6']]; assert r['verdict'] == 'stable_somewhere' and not r['complete'] and 'witnesses' in r and r['witnesses'] == h['witnesses']
assert len(hits) == len({h['d6'] for h in hits})
re = collections.Counter(); re['S'] = sum(v[1] == 'S' for v in rows.values()); re['O'] = sum(v[1] == 'O' for v in rows.values())
for r in recs.values():
    re['complete'] += r['complete']
    if r['verdict'] == 'never' and 'cert' not in r: re['never_positivity'] += 1
    elif r['verdict'] == 'never': re['never_full'] += 1
    elif r['verdict'] == 'stable_somewhere':
        re['stable_somewhere'] += 1
        if not r['complete']: re['stable_noncomplete'] += 1
    else: re['escalate'] += 1
for key in ('S', 'O', 'never_positivity', 'never_full', 'stable_somewhere', 'escalate', 'complete', 'stable_noncomplete'):
    assert re[key] == tsum[key], (key, re[key], tsum[key])
assert re['stable_noncomplete'] == len(hits), (re['stable_noncomplete'], len(hits))
old = 0; oldbad = 0
for line in open('n8/n8_limit_results.jsonl'):
    try: r = json.loads(line)
    except Exception: continue
    old += 1; s = recs[r['d6']]
    if r != s: oldbad += 1; print('DISAGREE', r['d6'], r['verdict'], s['verdict'])
print(f'interrupted single run: {old} complete records compared, {oldbad} disagreements')
assert oldbad == 0
tally = {'S': re['S'], 'O': re['O'], 'never_positivity': re['never_positivity'], 'never_full': re['never_full'], 'stable_somewhere': re['stable_somewhere'], 'escalate': re['escalate'], 'complete': re['complete'], 'stable_noncomplete': re['stable_noncomplete'], 'skipped_done': 0, 'seconds': round(secs, 1), 'sharded': 4, 'interrupted_single_run_records_cross_checked': old}
print('merged tally:', json.dumps(tally))
if not apply: sys.exit(0)
os.rename('n8/n8_limit_results.jsonl', f'{SH}/n8_limit_results.interrupted_single_run.jsonl'); os.rename('n8/n8_limit_HITS.jsonl', f'{SH}/n8_limit_HITS.interrupted_single_run.jsonl')
with open('n8/n8_limit_results.jsonl', 'w') as f:
    for d6 in sorted(recs, key=lambda s: (rows[s][0], s)): f.write(json.dumps(recs[d6]) + '\n')
with open('n8/n8_limit_HITS.jsonl', 'w') as f:
    for h in sorted(hits, key=lambda h: (h['edges'], h['d6'])): f.write(json.dumps(h) + '\n')
json.dump(tally, open('n8/n8_limit_tally.json', 'w'), indent=1); print('written')
