import numpy as np, sys, gzip, json
from first_order import analyse
from scan_n7pos import from_d6
rows = []
for line in open('/home/user/CLTN/n7/n6_survivors.txt'):
    s, e, f = line.split(); rows.append((s, f))
for line in gzip.open('/home/user/CLTN/n7/n7_survivors_all.txt.gz', 'rt'):
    s, e, f = line.split()
    if f == 'S': rows.append((s, f))
verdict = {}
for line in open('/home/user/CLTN/n7/n7_limit_results.jsonl'):
    r = json.loads(line); verdict[r['d6']] = r['verdict']
st = {}; ex = []
for s, f in rows:
    A = from_d6(s)
    r = analyse(A)
    if r is None: continue
    n = len(A)
    if abs(r['S']) > 1e-7: key = (n, 'fixed-eig', r['pos'])
    elif abs(r['k0'].imag) < 1e-9: key = (n, 'real', r['pos'])
    else: key = (n, 'cplx', 'pos' if r['pos'] else 'nonpos', 'good' if r['ReSp'] <= 0 else 'BAD')
    st[key] = st.get(key, 0) + 1
    if key[1] == 'cplx' and r['pos']:
        ex.append((s, n, verdict.get(s), round(r['beta_c'], 6), round(r['ReSp'], 6)))
for k in sorted(st, key=str): print(k, st[k])
for e in ex: print(e)
