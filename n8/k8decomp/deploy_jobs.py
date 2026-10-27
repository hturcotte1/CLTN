#!/usr/bin/env python3
"""Replace the 348 directg parts of K8 (graph G~~~~{, split p/348) in n8/jobs/jobs.tsv and graph_counts.json by
the K bidirected-subgraph batches of n8/k8decomp/batches (see notes/k8_decomposition.md).
   python3 n8/k8decomp/deploy_jobs.py [--apply]
Refuses if any K8 part is done or running (the decomposition replaces ALL of K8; a mixed scheme would double
count).  New rows get ids 12000.., expected = the batch's exact count from the dry run (manifest batch_expected; the
driver checks it), split = decomp:<batchfile>, exact False, so check_accounting.py also sums them against
the exact K8 count like any split graph.  Backups are kept.
Then: bash n8/stop_driver.sh; cp <patched driver>; bash n8/relaunch.sh."""
import json, os, re, subprocess, sys, shutil, time, glob
apply = '--apply' in sys.argv
J = 'n8/jobs'; K8 = 'G~~~~{'
rows = open(f'{J}/jobs.tsv').read().splitlines(); hdr, body = rows[0], rows[1:]
running = subprocess.run(['ps', '-eo', 'args'], capture_output=True, text=True).stdout
batches = sorted(glob.glob('n8/k8decomp/batches/batch*.txt'), key=lambda p: int(re.search(r'batch(\d+)', p).group(1)))
man = json.load(open('n8/k8decomp/batches/manifest.json')); K = man['K']; assert len(batches) == K and man['n'] == 8 and man['edges'] == '1:28'
gc = json.load(open(f'{J}/graph_counts.json')); exact = gc['graph_counts'][K8]
old = [l for l in body if l.split('\t')[3] == K8]; keep = [l for l in body if l.split('\t')[3] != K8]
assert len(old) == 348 and all('/348' in l.split('\t')[2] for l in old), (len(old), old[:2])
for l in old:
    j = l.split('\t')[0]
    assert not os.path.exists(f'{J}/job{j}.done'), f'K8 part {j} already done: cannot mix schemes'
    assert not re.search(rf'{J}/job{j}\.g6\b', running), f'K8 part {j} is running: stop it first'
exp = man['batch_expected']; assert sum(exp) == exact, (sum(exp), exact)   # exact per-batch counts from the dry run
new = [f'{12000 + k}\t{exp[k]}\tdecomp:{batches[k]}\t{K8}' for k in range(K)]
print(f'K8 exact {exact}; replacing {len(old)} directg parts (ids {old[0].split(chr(9))[0]}..{old[-1].split(chr(9))[0]}) by {K} batches (ids 12000..{12000 + K - 1}); rows {len(body)} -> {len(keep) + K}')
if not apply: sys.exit(0)
ts = time.strftime('%Y%m%dT%H%M%S'); shutil.copy(f'{J}/jobs.tsv', f'{J}/jobs.tsv.bak.{ts}'); shutil.copy(f'{J}/graph_counts.json', f'{J}/graph_counts.json.bak.{ts}')
open(f'{J}/jobs.tsv', 'w').write('\n'.join([hdr] + keep + new) + '\n')
sup = {int(l.split('\t')[0]) for l in old}
gc['jobs'] = [x for x in gc['jobs'] if x['job_id'] not in sup] + [{'expected': exp[k], 'split': f'decomp:{batches[k]}', 'graphs': [K8], 'exact': False, 'job_id': 12000 + k} for k in range(K)]
gc.setdefault('k8decomp', []).append({'when': ts, 'K': K, 'superseded': sorted(sup), 'added': list(range(12000, 12000 + K)), 'manifest': man})
json.dump(gc, open(f'{J}/graph_counts.json', 'w'))
for j in sup:
    for f in glob.glob(f'{J}/job{j}.*'): os.remove(f)
print('applied; backups', f'{J}/jobs.tsv.bak.{ts}', f'{J}/graph_counts.json.bak.{ts}')
