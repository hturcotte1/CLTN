import subprocess, collections, json, sys
W = '/home/user/CLTN/notes/verify_n8_2/chain_referee'
rows = [l.rstrip('\n').split('\t') for l in open('n8/jobs/jobs.tsv')][1:]
job_graphs = {jid: gs.split() for jid, ex, sp, gs in rows}
shard_graphs = {k: [l.strip() for l in open(f'{W}/shard{k}.g6') if l.strip()] for k in range(9)}
allg = sorted({g for v in job_graphs.values() for g in v} | {g for v in shard_graphs.values() for g in v})
open(f'{W}/all_job_graphs.g6', 'w').write('\n'.join(allg) + '\n')
can = subprocess.run(['n7/labelg', '-q', f'{W}/all_job_graphs.g6'], capture_output=True, text=True).stdout.split()
cmap = dict(zip(allg, can))
print('graphs whose geng labeling is not labelg-canonical:', sum(1 for g in allg if cmap[g] != g), 'of', len(allg), '; canonical forms distinct:', len(set(can)))
kept = []; owner = []
for k in range(9):
    for l in open(f'n8/parts/part{k}.txt'):
        p = l.split()
        if len(p) >= 3: kept.append(p[0]); owner.append(('shard', k))
for jid in job_graphs:
    for l in open(f'n8/jobs/job{jid}.txt'):
        p = l.split()
        if len(p) >= 3: kept.append(p[0]); owner.append(('job', jid))
canon = open(f'{W}/kept_underlying.g6').read().split()
canon = subprocess.run(['n7/labelg', '-q', f'{W}/kept_underlying.g6'], capture_output=True, text=True).stdout.split()
assert len(canon) == len(kept)
bad = collections.Counter(); per_graph = collections.Counter()
for d6, o, cg in zip(kept, owner, canon):
    lst = shard_graphs[o[1]] if o[0] == 'shard' else job_graphs[o[1]]
    per_graph[cg] += 1
    if cg not in {cmap[g] for g in lst}: bad[(o, cg)] += 1
print('kept lines whose underlying graph is not in the producing job list (after canonicalising both sides):', sum(bad.values()), list(bad.items())[:5])
print('kept S lines per underlying graph (canonical):', sorted(per_graph.items(), key=lambda x: -x[1]))
