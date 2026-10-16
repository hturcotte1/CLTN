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

def main():
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
