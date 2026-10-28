"""Referee check: do the jobs (jobs.tsv) + shards 0..8 cover every one of the 5,350 underlying graphs exactly once?
Split parts p/P: residue-class coverage of [0, lcm) exactly once per graph (directg -s nesting p/P == union (p+kP)/(mP))."""
import json, math, collections, os, sys
W = os.path.dirname(os.path.abspath(__file__))
full = [l.strip() for l in open(f'{W}/geng_8_15_28.g6') if l.strip()]
assert len(full) == 5350 and len(set(full)) == 5350
shard = {}
for k in range(9):
    for l in open(f'{W}/shard{k}.g6'):
        l = l.strip()
        if l: assert l not in shard; shard[l] = k
rc = json.load(open('n8/n8_restricted_counts.json'))
for k in range(9):
    assert sum(1 for g in shard if shard[g] == k) == rc['per_shard'][str(k)]['graphs'], k
rows = [l.rstrip('\n').split('\t') for l in open('n8/jobs/jobs.tsv')][1:]
gc = json.load(open('n8/jobs/graph_counts.json'))
unsplit = collections.Counter(); split = collections.defaultdict(list); decomp = collections.defaultdict(list)
ids = [int(r[0]) for r in rows]; assert len(ids) == len(set(ids)), 'duplicate job id'
for jid, ex, sp, gs in rows:
    gl = gs.split()
    if sp == '-':
        for g in gl: unsplit[g] += 1
    elif sp.startswith('decomp:'):
        assert len(gl) == 1; decomp[gl[0]].append(sp)
    else:
        assert len(gl) == 1
        kind = 'orbit' if sp.startswith('orbit:') else 'directg'
        p, P = map(int, sp.split(':')[-1].split('/'))
        split[gl[0]].append((kind, p, P))
    # every job has .done and .stats
    for ext in ('.done', '.stats'):
        assert os.path.exists(f'n8/jobs/job{jid}{ext}'), (jid, ext)
problems = []
for g, c in unsplit.items():
    if c != 1: problems.append(f'unsplit graph {g} appears {c} times')
cover = collections.Counter()
for g in unsplit: cover[g] += 1
for g, parts in split.items():
    kinds = {k for k, _, _ in parts}
    if len(kinds) != 1: problems.append(f'{g}: mixed enumerators {kinds}')
    L = 1
    for _, p, P in parts: L = L * P // math.gcd(L, P)
    cnt = [0] * L
    for _, p, P in parts:
        for x in range(p, L, P): cnt[x] += 1
    if any(c != 1 for c in cnt): problems.append(f'{g}: residue coverage not exactly once: min {min(cnt)} max {max(cnt)} parts {sorted(parts)}')
    cover[g] += 1
for g in decomp: cover[g] += 1
for g in shard: cover[g] += 1
missing = [g for g in full if cover[g] != 1]
if missing: problems.append(f'{len(missing)} graphs not covered exactly once: {missing[:10]} counts {[cover[g] for g in missing[:10]]}')
extra = [g for g in cover if g not in set(full)]
if extra: problems.append(f'graphs outside geng 8 15:28: {extra}')
# graph_counts.json consistency: every non-shard graph has an exact count, and sum == remaining_digraphs
rem = [g for g in full if g not in shard]
missing_cnt = [g for g in rem if g not in gc['graph_counts']]
if missing_cnt: problems.append(f'{len(missing_cnt)} remaining graphs without exact count')
s = sum(gc['graph_counts'][g] for g in rem if g in gc['graph_counts'])
shard_sum = sum(rc['per_shard'][str(k)]['digraphs'] for k in range(9))
print('unsplit graphs', len(unsplit), 'directg/orbit split graphs', len(split), 'decomp graphs', list(decomp), 'shard graphs', len(shard))
print('sum of exact counts of remaining graphs', s, 'shards 0..8', shard_sum, 'total', s + shard_sum, 'expected', rc['digraphs_total'])
if s + shard_sum != rc['digraphs_total']: problems.append('partition sums do not add to the space total')
# the jobs table in graph_counts.json agrees with jobs.tsv
tsv = {int(r[0]): (int(r[1]), r[2], r[3].split()) for r in rows}
gj = {j['job_id']: (j['expected'], j['split'], j['graphs']) for j in gc['jobs']}
if tsv != gj: problems.append(f'jobs.tsv and graph_counts.json jobs differ: {len(set(tsv) ^ set(gj))} ids differ, {sum(tsv[i] != gj[i] for i in set(tsv) & set(gj))} rows differ')
print('split-part multiplicities:', collections.Counter(len(v) for v in split.values()))
print('PROBLEMS:', problems if problems else 'none')
