#!/usr/bin/env python3
"""deploy_orbit.py [--apply] -- replace the directg parts of the high-symmetry graphs in TARGETS by orbitenum parts
(n8/orbitenum/orbitenum G 29 p/P: bidirected subgraphs B up to Aut(G), orientations up to the stabilizer; see
notes/k8_decomposition.md, "Relative form").  Each new row: job_id 13000+, expected = the part's exact count (the sum
of the "B <mask> <count>" lines with index = p (mod P) in the count-only dry run n8/orbitenum/dryrun/<name>.log, whose
total equals graph_counts.json), split = orbit:p/P, exact False (check_accounting.py then sums the parts of each
graph against its exact count; the driver's MISMATCH check applies to every part).  The new rows take the place of
the first old row of their graph, so the file order (largest first) is kept.
Refuses if any row of a target graph is done or running (the two partitions of a graph cannot be mixed) and, with
--apply, if the driver loop is alive.  ORDER MATTERS (the live loop has the whole table buffered):
  bash n8/stop_driver.sh;  python3 n8/orbitenum/deploy_orbit.py --apply;  cp <patched driver> n8/run_n8_jobs.sh;
  bash n8/relaunch.sh.
Set N8_JOBS_DIR to run against a copy of the two tables (for testing)."""
import json, os, re, shutil, subprocess, sys, time
apply = '--apply' in sys.argv
J = os.environ.get('N8_JOBS_DIR', 'n8/jobs'); D = 'n8/orbitenum/dryrun'; BASE = 13000
TARGETS = [('GF~~~{', 'K8-K3', 8), ('GTm~~{', 'K8-K15', 4), ('G]zn^[', 'K8-K2+K15', 2)]   # (graph6, dry-run name, P)
def driver_alive():
    try: pid = int(open(f'{J}/driver.pid').read().strip())
    except Exception: return False
    try: return 'run_n8_jobs' in open(f'/proc/{pid}/cmdline', 'rb').read().replace(b'\0', b' ').decode()
    except Exception: return False
def dryrun(name, g, exact):
    lines = open(f'{D}/{name}.log').read().splitlines()
    B = [(int(l.split()[1]), int(l.split()[2])) for l in lines if l.startswith('B ')]
    fin = [l for l in lines if l.startswith('orbitenum ')]
    assert len(fin) == 1 and fin[0].split()[1] == g and fin[0].split()[3] == '29' and fin[0].split()[5] == '0/1', fin
    assert int(fin[0].split()[-1]) == sum(c for _, c in B) == exact, (fin[0], sum(c for _, c in B), exact)
    assert int(fin[0].split()[9]) == len(B) and [b for b, _ in B] == sorted(set(b for b, _ in B)), 'dry run B lines'
    return B
rows = open(f'{J}/jobs.tsv').read().splitlines(); hdr, body = rows[0], rows[1:]
gc = json.load(open(f'{J}/graph_counts.json'))
running = subprocess.run(['ps', '-eo', 'args'], capture_output=True, text=True).stdout
used = {int(l.split('\t')[0]) for l in body}; assert not any(BASE <= u < BASE + 1000 for u in used), 'ids 13000+ in use'
plan = []; nid = BASE
for g, name, P in TARGETS:
    exact = gc['graph_counts'][g]; B = dryrun(name, g, exact)
    old = [l for l in body if l.split('\t')[3] == g]
    assert old and all(re.fullmatch(r'\d+/\d+', l.split('\t')[2]) for l in old), (g, old[:2])
    for l in old:
        j = l.split('\t')[0]
        assert not os.path.exists(f'{J}/job{j}.done'), f'{g} part {j} already done: cannot mix schemes'
        assert not re.search(rf'{re.escape(J)}/job{j}\.g6\b', running), f'{g} part {j} is running: cannot mix schemes'
    exp = [sum(c for i, (_, c) in enumerate(B) if i % P == p) for p in range(P)]; assert sum(exp) == exact
    ids = list(range(nid, nid + P)); nid += P
    plan.append({'graph': g, 'name': name, 'P': P, 'exact': exact, 'superseded': [int(l.split('\t')[0]) for l in old],
                 'added': ids, 'expected': exp, 'dryrun': f'{D}/{name}.log', 'dryrun_B_lines': len(B)})
    print(f'{g} ({name}): exact {exact:,}; replacing {len(old)} directg parts ({old[0].split(chr(9))[2].split("/")[1]} parts, ids {plan[-1]["superseded"][0]}..{plan[-1]["superseded"][-1]}) by {P} orbit parts (ids {ids[0]}..{ids[-1]}, expected {min(exp):,}..{max(exp):,})')
sup = {j for t in plan for j in t['superseded']}
new_rows = {t['superseded'][0]: [f'{t["added"][p]}\t{t["expected"][p]}\torbit:{p}/{t["P"]}\t{t["graph"]}' for p in range(t['P'])] for t in plan}
out = []
for l in body:
    j = int(l.split('\t')[0])
    if j in new_rows: out += new_rows[j]
    if j not in sup: out.append(l)
print(f'rows {len(body)} -> {len(out)}')
if not apply: sys.exit(0)
assert not driver_alive(), 'the driver loop is alive: run bash n8/stop_driver.sh first'
ts = time.strftime('%Y%m%dT%H%M%S'); shutil.copy(f'{J}/jobs.tsv', f'{J}/jobs.tsv.bak.{ts}'); shutil.copy(f'{J}/graph_counts.json', f'{J}/graph_counts.json.bak.{ts}')
open(f'{J}/jobs.tsv', 'w').write('\n'.join([hdr] + out) + '\n')
gc['jobs'] = [x for x in gc['jobs'] if x['job_id'] not in sup] + [{'expected': t['expected'][p], 'split': f'orbit:{p}/{t["P"]}', 'graphs': [t['graph']], 'exact': False, 'job_id': t['added'][p]} for t in plan for p in range(t['P'])]
gc.setdefault('orbitenum', []).append({'when': ts, 'min_arcs': 29, 'targets': plan})
json.dump(gc, open(f'{J}/graph_counts.json', 'w'))
import glob
for j in sup:
    for f in glob.glob(f'{J}/job{j}.*'): os.remove(f)
print('applied; backups', f'{J}/jobs.tsv.bak.{ts}', f'{J}/graph_counts.json.bak.{ts}')
