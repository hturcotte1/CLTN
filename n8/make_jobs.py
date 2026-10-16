#!/usr/bin/env python3
"""Build the (underlying graph, directg part) job list for the n = 8 prefilter.

geng's res/mod sharding is useless here: with the 15:28 edge window, 3,966 of 4,000 shards are empty and
one shard holds 19.5% of all digraphs (n8/n8_restricted_counts.json).  Shards 0..8 of that scheme are
finished and verified (filter8 'total' == exact Burnside count, every shard); the remaining underlying
graphs (geng 8 15:28 minus the graphs of shards 0..8) are packed here into jobs of at most TARGET
digraphs, using the exact per-graph counts from count_restricted.py.  A graph with more than TARGET
digraphs is split with  directg -s p/P  (P = ceil(count/TARGET)); directg splits per input graph, each
orientation lands in exactly one part, so the sum over a graph's parts is exact.

Writes n8/jobs/jobs.tsv (job_id, expected_digraphs, split, graphs), largest jobs first, and
n8/jobs/graph_counts.json (exact per-graph counts + the job table) for check_accounting.py.
Integers only; no model in the loop.
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from count_restricted import pair_index, perm_tables, graph6_edges, count_graph, geng

TARGET = 1_500_000_000
DONE_SHARDS = list(range(9)); K = 4000; N = 8; EDGES = '15:28'; LO = 29
SPACE_TOTAL = 801_558_826_511          # n8/n8_restricted_counts.json

def pack(graphs, cnt, target):
    """pack graphs (in order) into jobs of <= target digraphs; a graph above target is split into
    ceil(count/target) directg parts.  Returns job dicts without ids."""
    jobs = []; cur, cur_sum = [], 0
    for g in graphs:
        c = cnt[g]
        if c > target:
            P = math.ceil(c / target)
            for p in range(P):
                jobs.append({'expected': c // P, 'split': f'{p}/{P}', 'graphs': [g], 'exact': False})
            continue
        if cur and cur_sum + c > target:
            jobs.append({'expected': cur_sum, 'split': '-', 'graphs': cur, 'exact': True}); cur, cur_sum = [], 0
        cur.append(g); cur_sum += c
    if cur: jobs.append({'expected': cur_sum, 'split': '-', 'graphs': cur, 'exact': True})
    return jobs

def write_tables(jobs, gc_extra):
    with open('n8/jobs/jobs.tsv', 'w') as fh:
        fh.write('job_id\texpected_digraphs\tsplit\tgraphs\n')
        for j in jobs:
            fh.write(f"{j['job_id']}\t{j['expected']}\t{j['split']}\t{' '.join(j['graphs'])}\n")
    gc_extra['jobs'] = jobs
    json.dump(gc_extra, open('n8/jobs/graph_counts.json', 'w'), indent=0)

def repack(target):
    """Re-pack the UNFINISHED work into smaller jobs (container restarts lose in-flight jobs).
    Finished jobs keep their ids and rows.  A split graph with some parts finished keeps its P and the
    unfinished parts become jobs of their own (a part cannot be sub-split); every other unfinished graph
    is packed afresh with the new target.  Unfinished old ids disappear (their partial files are removed);
    new jobs get ids above every old id, so an id never changes meaning."""
    gc = json.load(open('n8/jobs/graph_counts.json'))
    old = gc['jobs']; cnt = gc['graph_counts']; order = list(cnt.keys())      # geng order
    done = {j['job_id'] for j in old if os.path.exists(f"n8/jobs/job{j['job_id']}.done")}
    kept = [j for j in old if j['job_id'] in done]
    covered = set()
    for j in kept:
        if j['exact']: covered |= set(j['graphs'])
    next_id = max(j['job_id'] for j in old) + 1
    extra = []
    parts_of = {}
    for j in old:
        if not j['exact']: parts_of.setdefault(j['graphs'][0], []).append(j)
    for g, parts in parts_of.items():
        nd = sum(p['job_id'] in done for p in parts)
        if nd == len(parts): covered.add(g)
        elif nd > 0:
            for p in parts:
                if p['job_id'] not in done:
                    extra.append(dict(p, job_id=next_id)); next_id += 1
            covered.add(g)
    new = pack([g for g in order if g not in covered], cnt, target)
    new.sort(key=lambda j: -j['expected'])
    for j in new: j['job_id'] = next_id; next_id += 1
    for j in old:
        if j['job_id'] not in done:
            for ext in ('.stats', '.txt', '.txt.tmp', '.g6', '.bad'):
                try: os.remove(f"n8/jobs/job{j['job_id']}{ext}")
                except FileNotFoundError: pass
    jobs = kept + extra + new
    split_present = {j['graphs'][0] for j in jobs if not j['exact']}
    assert sum(cnt.values()) == sum(j['expected'] for j in jobs if j['exact']) + sum(cnt[g] for g in split_present), 'repack does not partition the remaining space'
    assert len({j['job_id'] for j in jobs}) == len(jobs), 'duplicate job id'
    total_new = sum(cnt[g] for g in order if g not in covered)
    gc.update({'target': target, 'repacked': gc.get('repacked', []) + [{'kept_done_jobs': len(kept), 'kept_split_parts': len(extra), 'new_jobs': len(new), 'new_digraphs': total_new}]})
    write_tables(jobs, gc)
    print(f'repack: kept {len(kept)} finished jobs, {len(extra)} unfinished parts of partly finished split graphs, '
          f'{len(new)} new jobs (target {target:,}) over {total_new:,} digraphs; table now {len(jobs)} jobs')

def main():
    if len(sys.argv) >= 3 and sys.argv[1] == '--repack':
        repack(int(sys.argv[2])); return
    full = geng(N, EDGES)
    done = set()
    for k in DONE_SHARDS:
        done |= set(geng(N, EDGES, f'{k}/{K}'))
    assert done <= set(full), 'shard graphs not in the full list'
    rem = [g for g in full if g not in done]
    idx = pair_index(N); EP, EF = perm_tables(N, idx)
    cnt = {g: count_graph(graph6_edges(g, N)[0], EP, EF, LO)[0] for g in rem}
    shard_counts = json.load(open('n8/n8_restricted_counts.json'))['per_shard']
    done_total = sum(shard_counts[str(k)]['digraphs'] for k in DONE_SHARDS)
    assert done_total + sum(cnt.values()) == SPACE_TOTAL, 'partition does not add up'
    jobs = []            # dicts: expected, split, graphs, graph (for split parts)
    cur, cur_sum = [], 0
    for g in rem:
        c = cnt[g]
        if c > TARGET:
            P = math.ceil(c / TARGET)
            for p in range(P):
                jobs.append({'expected': c // P, 'split': f'{p}/{P}', 'graphs': [g], 'exact': False})
            continue
        if cur and cur_sum + c > TARGET:
            jobs.append({'expected': cur_sum, 'split': '-', 'graphs': cur, 'exact': True}); cur, cur_sum = [], 0
        cur.append(g); cur_sum += c
    if cur: jobs.append({'expected': cur_sum, 'split': '-', 'graphs': cur, 'exact': True})
    jobs.sort(key=lambda j: -j['expected'])
    for i, j in enumerate(jobs): j['job_id'] = i
    os.makedirs('n8/jobs', exist_ok=True)
    with open('n8/jobs/jobs.tsv', 'w') as fh:
        fh.write('job_id\texpected_digraphs\tsplit\tgraphs\n')
        for j in jobs:
            fh.write(f"{j['job_id']}\t{j['expected']}\t{j['split']}\t{' '.join(j['graphs'])}\n")
    json.dump({'n': N, 'geng_edges': EDGES, 'min_arcs': LO, 'target': TARGET, 'done_shards': DONE_SHARDS, 'K': K,
               'done_shards_digraphs': done_total, 'remaining_graphs': len(rem), 'remaining_digraphs': sum(cnt.values()),
               'space_total': SPACE_TOTAL, 'graph_counts': cnt, 'jobs': jobs},
              open('n8/jobs/graph_counts.json', 'w'), indent=0)
    nsplit = sum(not j['exact'] for j in jobs)
    print(f'full graphs {len(full)}, done in shards 0..8: {len(done)}, remaining {len(rem)} graphs = {sum(cnt.values()):,} digraphs')
    print(f'jobs {len(jobs)} ({nsplit} split parts of {len(set(j["graphs"][0] for j in jobs if not j["exact"]))} big graphs); '
          f'largest {jobs[0]["expected"]:,}, median {jobs[len(jobs)//2]["expected"]:,}, smallest {jobs[-1]["expected"]:,}')

if __name__ == '__main__':
    main()
