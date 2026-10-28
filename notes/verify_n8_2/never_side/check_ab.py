"""(a) own exact alpha class vs filter8 flag; (b) own domination-freeness vs filter8 printing."""
import json, sys, time
sys.path.insert(0, '.')
from engine import parse_d6, is_domination_free, ones_arcs
from alpha import classify_alpha
res = {}; rows = []
t0 = time.time()
for f, expect_df in [('sample_S.jsonl', True), ('sample_B.jsonl', True), ('sample_NPDF.jsonl', True), ('sample_NDF.jsonl', False), ('sample_stride_domfree.jsonl', True)]:
    cnt = {}
    for line in open(f):
        r = json.loads(line); A = parse_d6(r['d6'])
        df = is_domination_free(A)
        c, det, g = classify_alpha(A) if df else (None, 'not domination-free (alpha not needed)', None)
        key = (r['cls'], df, c); cnt[key] = cnt.get(key, 0) + 1
        rows.append({'d6': r['d6'], 'cls': r['cls'], 'arcs': ones_arcs(A), 'domfree': df, 'alpha_vs_minus_half': c, 'detail': det, 'g': g, 'filter8_flag': r.get('filter8_flag')})
    res[f] = {str(k): v for k, v in cnt.items()}
    print(f, cnt, flush=True)
json.dump({'tables': res, 'seconds': round(time.time() - t0, 1)}, open('check_ab_summary.json', 'w'), indent=1)
with open('check_ab_rows.jsonl', 'w') as fh:
    for r in rows: fh.write(json.dumps(r) + '\n')
