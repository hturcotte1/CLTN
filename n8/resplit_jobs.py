#!/usr/bin/env python3
"""Re-split the unfinished, not-running split jobs of the n = 8 run m times finer (directg -s nesting verified
2026-10-22: part p/P == union of parts (p+kP)/(mP), k < m, by set equality on two 17-edge graphs, P = 5 and 89).
Usage: python3 resplit_jobs.py [m] [--apply]   (default m = 3; without --apply only prints the plan)
After --apply: bash n8/stop_driver.sh; bash n8/relaunch.sh   (in-flight jobs keep running; the driver skips them)."""
import json, os, re, subprocess, sys, shutil, time
m = 3; apply = False
for a in sys.argv[1:]:
    if a == '--apply': apply = True
    else: m = int(a)
J = 'n8/jobs'; rows = open(f'{J}/jobs.tsv').read().splitlines(); hdr, body = rows[0], rows[1:]
running = subprocess.run(['pgrep', '-f', '-a', 'directg'], capture_output=True, text=True).stdout
new_rows = []; superseded = []; added = []; next_id = 10000
for l in body:
    j, ex, sp, gs = l.split('\t')
    if os.path.exists(f'{J}/job{j}.done') or not re.fullmatch(r'\d+/\d+', sp) or re.search(rf'{J}/job{j}\.g6\b', running):   # decomp:/orbit: rows pass through unchanged
        new_rows.append(l); continue
    p, P = map(int, sp.split('/'))
    superseded.append(j)
    for k in range(m):
        nid = str(next_id); next_id += 1
        new_rows.append('\t'.join([nid, str(int(ex) // m), f'{p + k * P}/{P * m}', gs])); added.append((nid, int(ex) // m, f'{p + k * P}/{P * m}', gs))
print(f'superseding {len(superseded)} unfinished split jobs by {len(added)} sub-part jobs (m = {m}); rows {len(body)} -> {len(new_rows)}')
if not apply: sys.exit(0)
ts = time.strftime('%Y%m%dT%H%M%S'); shutil.copy(f'{J}/jobs.tsv', f'{J}/jobs.tsv.bak.{ts}'); shutil.copy(f'{J}/graph_counts.json', f'{J}/graph_counts.json.bak.{ts}')
open(f'{J}/jobs.tsv', 'w').write('\n'.join([hdr] + new_rows) + '\n')
gc = json.load(open(f'{J}/graph_counts.json')); sup = set(superseded)
gc['jobs'] = [x for x in gc['jobs'] if str(x['job_id']) not in sup] + [{'expected': e, 'split': s, 'graphs': g.split(), 'exact': False, 'job_id': int(i)} for i, e, s, g in added]
gc.setdefault('resplit', []).append({'when': ts, 'm': m, 'superseded': superseded, 'added': [a[0] for a in added]})
json.dump(gc, open(f'{J}/graph_counts.json', 'w'))
print('applied; backups', f'{J}/jobs.tsv.bak.{ts}', f'{J}/graph_counts.json.bak.{ts}')
