#!/usr/bin/env python3
"""check_orbit.py -- post-run ledger for the orbitenum parts deployed by deploy_orbit.py (records in
graph_counts.json['orbitenum']).  For every deployed part: job<j>.done exists; the part's genlog (orbitenum's stderr)
has the B lines of exactly the dry-run representatives with index = p (mod P), with the same counts, and its final
line's total equals the sum of those counts, the filter8 'total' in job<j>.stats and the expected value in jobs.tsv;
over the P parts of a graph every dry-run B appears exactly once and the grand total equals the exact count in
graph_counts.json.  Prints one line per graph and 'orbitenum ledger ok' (exit 0) or the problems (exit 1).
Set N8_JOBS_DIR to check a copy."""
import json, os, re, sys
J = os.environ.get('N8_JOBS_DIR', 'n8/jobs'); problems = []
gc = json.load(open(f'{J}/graph_counts.json'))
rows = {int(l.split('\t')[0]): l.split('\t') for l in open(f'{J}/jobs.tsv').read().splitlines()[1:]}
def blines(path):
    L = open(path).read().splitlines()
    return [(int(l.split()[1]), int(l.split()[2])) for l in L if l.startswith('B ')], [l for l in L if l.startswith('orbitenum ')]
for rec in gc.get('orbitenum', []):
    for t in rec['targets']:
        g, P, exact = t['graph'], t['P'], gc['graph_counts'][t['graph']]
        dry, fin = blines(t['dryrun']); assert len(fin) == 1 and int(fin[0].split()[-1]) == exact == sum(c for _, c in dry), (g, 'dry run')
        seen = []; grand = 0
        for p, j in enumerate(t['added']):
            if j not in rows or rows[j][2] != f'orbit:{p}/{P}' or rows[j][3] != g: problems.append(f'{g}: row {j} missing or not orbit:{p}/{P}'); continue
            if not os.path.exists(f'{J}/job{j}.done'): problems.append(f'{g}: part {j} not done'); continue
            want = [dry[i] for i in range(len(dry)) if i % P == p]
            if not (os.path.exists(f'{J}/job{j}.genlog') and os.path.exists(f'{J}/job{j}.stats')): problems.append(f'{g}: part {j} genlog or stats file missing'); continue
            got, gfin = blines(f'{J}/job{j}.genlog')
            if got != want: problems.append(f'{g}: part {j} genlog B lines differ from the dry run ({len(got)} vs {len(want)})')
            tot = sum(c for _, c in got)
            if len(gfin) != 1 or gfin[0].split()[5] != f'{p}/{P}' or int(gfin[0].split()[-1]) != tot: problems.append(f'{g}: part {j} final line {gfin}')
            st = open(f'{J}/job{j}.stats').read().split(); ftot = int(st[st.index('total') + 1])
            if not (tot == ftot == int(rows[j][1]) == t['expected'][p]): problems.append(f'{g}: part {j} totals genlog {tot} filter8 {ftot} jobs.tsv {rows[j][1]} record {t["expected"][p]}')
            seen += [b for b, _ in got]; grand += tot
        if sorted(seen) != [b for b, _ in dry]: problems.append(f'{g}: union of parts != dry-run representatives ({len(seen)} vs {len(dry)})')
        if grand != exact: problems.append(f'{g}: grand total {grand} != exact {exact}')
        print(f'{g} ({t["name"]}): {P} parts, {len(seen)} representatives, total {grand:,} (exact {exact:,}) {"ok" if grand == exact and sorted(seen) == [b for b, _ in dry] else "PROBLEM"}')
if problems: print('\n'.join(problems)); sys.exit(1)
print('orbitenum ledger ok')
