"""Pass 2: every stored pointwise witness in the limit file, plus first LEGAL witness per graph.
Reuses decode_d6/check/is_complete_bidirected from verify.py (own code). Exact only.
Usage: python3 -I all_witnesses.py <repo_root> <out_dir>"""
import sys, json, types
src = open(__file__.replace('all_witnesses.py', 'verify.py')).read().replace('\nmain()\n', '\n')
mod = types.ModuleType('v'); mod.__dict__['__file__'] = 'verify.py'
sys_argv = sys.argv
exec(compile(src, 'verify.py', 'exec'), mod.__dict__)
ROOT, OUT = sys_argv[1], sys_argv[2]
lim = [json.loads(l) for l in open(f'{ROOT}/certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl')]
rows = []; mism = []; per_graph = {}
for r in lim:
    A = mod.decode_d6(r['d6']); assert A == r['A']
    cb = mod.is_complete_bidirected(A)
    for i, w in enumerate(r['pointwise_witnesses']):
        res = mod.check(A, w['eps'], w['delta'])
        stored = {'legal': w['legal'], 'positive': w['is_fixed_point'], 'stable': w['stable']}
        mine = {'legal': res['legal'], 'positive': res['positive'], 'stable': res['stable']}
        x_match = (res['x'] == w.get('x')) if w.get('x') is not None else None
        agree = all(stored[k] == mine[k] for k in stored) and x_match in (True, None)
        cex = res['legal'] is True and res['positive'] is True and res['stable'] is True and not cb
        row = {'d6': r['d6'], 'idx': i, 'eps': res['eps'], 'delta': res['delta'], 'q': res['q'], **mine,
               'clique': cb, 'stored': stored, 'x_match': x_match, 'agrees_with_stored': agree,
               'routh_agrees': res['routh_agrees'], 'counterexample_here': cex}
        rows.append(row)
        if not agree: mism.append(row)
        per_graph.setdefault(r['d6'], []).append(cex)
    print(r['d6'], per_graph[r['d6']], flush=True)
no_legal_witness = [d for d, v in per_graph.items() if not any(v)]
summ = {'n_graphs': len(lim), 'n_points': len(rows), 'n_mismatch_with_stored_flags': len(mism),
        'mismatches': mism, 'graphs_with_no_confirming_witness': no_legal_witness,
        'n_points_counterexample': sum(r['counterexample_here'] for r in rows),
        'n_undecided': sum(1 for r in rows if r['stable'] == 'UNDECIDED')}
json.dump({'summary': summ, 'rows': rows}, open(f'{OUT}/v_sympy_all_witnesses.json', 'w'), indent=1)
print('SUMMARY', json.dumps(summ))
