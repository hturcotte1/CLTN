#!/usr/bin/env python3
"""Close the n = 8 prefilter accounting chain.  Every digraph of the restricted space (8 nodes, >= 29 arcs,
underlying graph with 15..28 edges; 801,558,826,511 of them by the exact Burnside count in
n8/n8_restricted_counts.json) must have been read by filter8 exactly once:
  * geng shards 0..8 (n8/parts):  'total' in part<k>.stats == per-shard Burnside count;
  * jobs (n8/jobs): unsplit job 'total' == exact expected; a split graph's parts sum to its exact count;
  * grand total == 801,558,826,511.   Also reports survivors / boundary / overflow and the kept S/O line counts.
Exit 0 iff everything matches and every shard/job has its .done marker.  Writes n8/n8_accounting.json.
"""
import json, re, os, sys, glob
def stats(path):
    m = re.search(r'total (\d+) domfree (\d+) survivors (\d+) boundary (\d+) overflow (\d+)', open(path).read())
    return [int(x) for x in m.groups()] if m else None
def main():
    ok = True; problems = []
    shard = json.load(open('n8/n8_restricted_counts.json')); gc = json.load(open('n8/jobs/graph_counts.json'))
    agg = {'total': 0, 'domfree': 0, 'survivors': 0, 'boundary': 0, 'overflow': 0}
    def add(s):
        for k, v in zip(agg, s): agg[k] += v
    for k in gc['done_shards']:
        p = f'n8/parts/part{k}'
        if not os.path.exists(p + '.done') or not os.path.exists(p + '.stats'): ok = False; problems.append(f'shard {k} not done'); continue
        s = stats(p + '.stats')
        if s is None or s[0] != shard['per_shard'][str(k)]['digraphs']: ok = False; problems.append(f'shard {k} total mismatch {s}')
        if s: add(s)
    split_sum = {}
    for j in gc['jobs']:
        p = f"n8/jobs/job{j['job_id']}"
        if not os.path.exists(p + '.done') or not os.path.exists(p + '.stats'): ok = False; problems.append(f"job {j['job_id']} not done"); continue
        s = stats(p + '.stats')
        if s is None: ok = False; problems.append(f"job {j['job_id']} unreadable stats"); continue
        add(s)
        if j['exact']:
            if s[0] != j['expected']: ok = False; problems.append(f"job {j['job_id']} total {s[0]} != expected {j['expected']}")
        else:
            split_sum[j['graphs'][0]] = split_sum.get(j['graphs'][0], 0) + s[0]
    for g, tot in split_sum.items():
        if tot != gc['graph_counts'][g]: ok = False; problems.append(f'split graph {g}: parts sum {tot} != exact {gc["graph_counts"][g]}')
    if agg['total'] != gc['space_total']: ok = False; problems.append(f"grand total {agg['total']} != {gc['space_total']}")
    kept = {'S': 0, 'O': 0, 'other': 0}
    for f in glob.glob('n8/parts/part*.txt') + glob.glob('n8/jobs/job*.txt'):
        for line in open(f):
            p = line.split()
            if len(p) >= 3: kept[p[2] if p[2] in kept else 'other'] += 1
    out = {'ok': ok, 'problems': problems[:50], 'n_problems': len(problems), 'aggregate': agg, 'space_total': gc['space_total'],
           'kept_lines': kept, 'shards_checked': len(gc['done_shards']), 'jobs_checked': len(gc['jobs']), 'split_graphs_checked': len(split_sum)}
    json.dump(out, open('n8/n8_accounting.json', 'w'), indent=1); print(json.dumps(out, indent=1))
    sys.exit(0 if ok else 1)
if __name__ == '__main__':
    main()
