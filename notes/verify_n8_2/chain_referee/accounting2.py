"""Referee's own accounting + content check of the n = 8 prefilter stage."""
import json, re, os, glob, subprocess, sys, collections
sys.path.insert(0, '/home/user/CLTN/notes/verify_n8_2/chain_referee')
from burnside import count, g6_edges
from myfilter2 import decode_d6
W = '/home/user/CLTN/notes/verify_n8_2/chain_referee'
rc = json.load(open('n8/n8_restricted_counts.json')); gc = json.load(open('n8/jobs/graph_counts.json'))
rows = [l.rstrip('\n').split('\t') for l in open('n8/jobs/jobs.tsv')][1:]
def stats(p):
    m = re.search(r'total (\d+) domfree (\d+) survivors (\d+) boundary (\d+) overflow (\d+)(?: \(charpoly (\d+) shift (\d+) hurwitz (\d+)\))?', open(p).read())
    return [int(x) if x is not None else None for x in m.groups()]
problems = []
agg = [0] * 5; stage_ovf = [0, 0, 0]; nostage = 0
shard_graphs = {k: [l.strip() for l in open(f'{W}/shard{k}.g6') if l.strip()] for k in range(9)}
# shards
for k in range(9):
    s = stats(f'n8/parts/part{k}.stats')
    if s[0] != rc['per_shard'][str(k)]['digraphs']: problems.append(f'shard {k} total {s[0]} != {rc["per_shard"][str(k)]["digraphs"]}')
    for i in range(5): agg[i] += s[i]
    if s[5] is None: nostage += 1
    else:
        for i in range(3): stage_ovf[i] += s[5 + i]
# jobs
split_sum = collections.defaultdict(int); decomp_sum = 0; orbit_sum = collections.defaultdict(int)
job_graphs = {}
for jid, ex, sp, gs in rows:
    s = stats(f'n8/jobs/job{jid}.stats')
    for i in range(5): agg[i] += s[i]
    if s[5] is None: nostage += 1
    else:
        for i in range(3): stage_ovf[i] += s[5 + i]
    gl = gs.split(); job_graphs[jid] = gl
    g6 = [l.strip() for l in open(f'n8/jobs/job{jid}.g6') if l.strip()]
    if g6 != gl: problems.append(f'job {jid}: job.g6 input list differs from jobs.tsv')
    if sp == '-':
        want = sum(gc['graph_counts'][g] for g in gl)
        if s[0] != want: problems.append(f'job {jid}: total {s[0]} != exact {want}')
    elif sp.startswith('decomp:'): decomp_sum += s[0]
    elif sp.startswith('orbit:'): orbit_sum[gl[0]] += s[0]
    else: split_sum[gl[0]] += s[0]
for g, t in split_sum.items():
    if t != gc['graph_counts'][g]: problems.append(f'split graph {g}: parts sum {t} != exact {gc["graph_counts"][g]}')
for g, t in orbit_sum.items():
    if t != gc['graph_counts'][g]: problems.append(f'orbit graph {g}: parts sum {t} != exact {gc["graph_counts"][g]}')
if decomp_sum != gc['graph_counts']['G~~~~{']: problems.append(f'K8 decomp sum {decomp_sum} != {gc["graph_counts"]["G~~~~{"]}')
tot_file = stats('n8/n8_prefilter_TOTAL.txt')
print('my aggregate (total, domfree, survivors, boundary, overflow):', agg, '; n8_prefilter_TOTAL.txt:', tot_file[:5], '; per-stage overflow (charpoly, shift, hurwitz):', stage_ovf, '; stats files without stage breakdown:', nostage)
if agg != tot_file[:5]: problems.append('aggregate differs from n8_prefilter_TOTAL.txt')
if agg[0] != rc['digraphs_total']: problems.append('grand total != Burnside space total')
# my own Burnside for every graph with >= 24 edges (where all strict survivors live) incl. K8
dense = subprocess.run(['n7/geng', '-q', '8', '24:28'], capture_output=True, text=True).stdout.split()
dense_res = {}
for g in dense:
    c = count(g, 29)[0]
    rec = gc['graph_counts'].get(g)
    if rec is None:
        ks = [k for k in range(9) if g in shard_graphs[k]]
        rec = ('shard', ks)
    dense_res[g] = (len(g6_edges(g)[1]), c, rec)
    if isinstance(rec, int) and rec != c: problems.append(f'dense graph {g}: my count {c} != recorded {rec}')
print('dense graphs (edges, my count, recorded):', dense_res)
# content: every kept line's underlying graph must be in the job's/shard's graph list; kept S set == results set
def underlying_g6(d6):
    n, A = decode_d6(d6)
    bits = []
    for j in range(1, n):
        for i in range(j): bits.append(1 if (A[i][j] or A[j][i]) else 0)
    while len(bits) % 6: bits.append(0)
    return chr(63 + n) + ''.join(chr(63 + int(''.join(map(str, bits[k:k+6])), 2)) for k in range(0, len(bits), 6))
kept = []; owner = []
for k in range(9):
    for l in open(f'n8/parts/part{k}.txt'):
        p = l.split()
        if len(p) >= 3: kept.append(p); owner.append(('shard', k))
for jid in job_graphs:
    for l in open(f'n8/jobs/job{jid}.txt'):
        p = l.split()
        if len(p) >= 3: kept.append(p); owner.append(('job', jid))
flags = collections.Counter(p[2] for p in kept)
print('kept lines by flag:', dict(flags), 'from shards:', sum(1 for o in owner if o[0] == 'shard'))
ug = [underlying_g6(p[0]) for p in kept]
open(f'{W}/kept_underlying.g6', 'w').write('\n'.join(ug) + '\n')
canon = subprocess.run(['n7/labelg', '-q', f'{W}/kept_underlying.g6'], capture_output=True, text=True).stdout.split()
assert len(canon) == len(kept)
bad_owner = 0; arcs_bad = 0
for p, o, cg in zip(kept, owner, canon):
    lst = shard_graphs[o[1]] if o[0] == 'shard' else job_graphs[o[1]]
    if cg not in lst: bad_owner += 1
    n, A = decode_d6(p[0])
    if sum(map(sum, A)) != int(p[1]) or int(p[1]) < 29: arcs_bad += 1
if bad_owner: problems.append(f'{bad_owner} kept lines whose underlying graph is not in the producing job/shard list')
if arcs_bad: problems.append(f'{arcs_bad} kept lines with wrong arc field or < 29 arcs')
res = {}
for l in open('n8/n8_limit_results.jsonl'):
    r = json.loads(l); res[r['d6']] = r
kept_set = {p[0] for p in kept}
if len(kept_set) != len(kept): problems.append(f'duplicate kept d6 lines: {len(kept) - len(kept_set)}')
if set(res) != kept_set: problems.append(f'limit results d6 set != kept S/O set: {len(set(res) ^ kept_set)} differ')
hits = [json.loads(l) for l in open('n8/n8_limit_HITS.jsonl')]
hs = {h['d6'] for h in hits}
cand = {d for d, r in res.items() if r['verdict'] == 'stable_somewhere' and not r['complete']}
if hs != cand or len(hits) != len(hs): problems.append('HITS set != non-complete stable_somewhere set')
esc = sum(1 for r in res.values() if r['verdict'] not in ('never', 'stable_somewhere'))
summ = json.load(open('certificates/n8_limit_summary.json'))
print('results', len(res), 'hits', len(hits), 'escalate/other verdicts', esc, 'complete graphs in results', sum(r['complete'] for r in res.values()))
print('summary keys', list(summ.keys())[:20])
print('PROBLEMS:', problems if problems else 'none')
json.dump({'aggregate': agg, 'stage_overflow': stage_ovf, 'dense': dense_res, 'kept_flags': dict(flags), 'problems': problems}, open(f'{W}/accounting2_result.json', 'w'), indent=1, default=str)
