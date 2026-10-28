"""Exact legal-ray grid test with the OWN engine (engine.py).
   python3 grid.py <sample.jsonl> <out.jsonl> <J>      grid beta = 2 + j/40, j = 1..J, t in {10000, 100000},
plus every record's extra_betas (critical-interval endpoints and midpoints) at both t.  Illegal points
(q <= 1, i.e. beta <= 2 + 1/t) are skipped and counted.  Resumable (d6 already in out are skipped).
Per graph: counts, and the full detail of every stable+positive point (FATAL) and every undecided point."""
import sys, json, os, time
from fractions import Fraction as F
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import parse_d6, decide_point, is_complete, ones_arcs
TS = (10000, 100000)

def run_graph(args):
    rec, J = args
    A = parse_d6(rec['d6'])
    pts = [(40 * 2 + j, 40) for j in range(1, J + 1)] + [tuple(p) for p in rec.get('extra_betas', [])]
    out = {'d6': rec['d6'], 'cls': rec['cls'], 'arcs': ones_arcs(A), 'complete': is_complete(A), 'n_grid': J,
           'n_extra': len(rec.get('extra_betas', [])), 'n_points': 0, 'n_illegal_skipped': 0, 'n_stable': 0,
           'n_positive': 0, 'n_stable_positive': 0, 'n_stability_undecided': 0, 'n_singular': 0,
           'stable_positive_points': [], 'undecided_points': [], 'stable_points': [], 'first_fail_minor_hist': {}}
    t0 = time.time()
    for t in TS:
        for bn, bd in pts:
            beta = F(bn, bd)
            if not (beta > 2 + F(1, t)):
                out['n_illegal_skipped'] += 1; continue
            r = decide_point(A, t, bn, bd)
            assert r['legal']
            out['n_points'] += 1
            key = r['stable_detail'].split(' ')[0] + (r['stable_detail'].split(' ')[2] if r['stable_detail'].startswith('Hurwitz') else '')
            out['first_fail_minor_hist'][key] = out['first_fail_minor_hist'].get(key, 0) + 1
            if r['stable'] is None:
                out['n_stability_undecided'] += 1
                out['undecided_points'].append({'t': t, 'beta': f"{bn}/{bd}", 'detail': r['stable_detail'], 'positive': r['positive']})
            if r['positive'] is None: out['n_singular'] += 1
            if r['stable'] is True:
                out['n_stable'] += 1
                out['stable_points'].append({'t': t, 'beta': f"{bn}/{bd}", 'positive': r['positive'],
                                             'x': [str(v) for v in r['x']] if r['x'] else None})
            if r['positive'] is True: out['n_positive'] += 1
            if r['stable'] is True and r['positive'] is True:
                out['n_stable_positive'] += 1
                out['stable_positive_points'].append({'t': t, 'beta': f"{bn}/{bd}", 'eps': str(r['eps']), 'delta': str(r['delta']),
                                                      'q': str(r['q']), 'x': [str(v) for v in r['x']], 'charpoly_of_scaled_M': r['charpoly']})
    out['seconds'] = round(time.time() - t0, 2)
    return out

if __name__ == '__main__':
    sample, outp, J = sys.argv[1], sys.argv[2], int(sys.argv[3])
    done = set()
    if os.path.exists(outp):
        for line in open(outp):
            try: done.add(json.loads(line)['d6'])
            except Exception: pass
    recs = [json.loads(l) for l in open(sample)]
    todo = [(r, J) for r in recs if r['d6'] not in done]
    print(f"{sample}: {len(recs)} graphs, {len(done)} done, {len(todo)} to do", flush=True)
    t0 = time.time(); k = 0
    with Pool(2) as pool, open(outp, 'a') as fh:
        for out in pool.imap_unordered(run_graph, todo, chunksize=1):
            fh.write(json.dumps(out) + '\n'); fh.flush(); k += 1
            if out['n_stable_positive']:
                print("*** FATAL: stable+positive legal point for", out['cls'], out['d6'], out['stable_positive_points'][0], flush=True)
            if k % 20 == 0:
                print(f"  {k}/{len(todo)} done, {time.time()-t0:.0f}s", flush=True)
    print("DONE", f"{time.time()-t0:.0f}s", flush=True)
