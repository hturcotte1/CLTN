"""
EXACT sharded driver: Lemma-R curve test (fast path) over the n=7 positivity survivors.

  python3 run_curve_n7.py --input n7_pos_all_PU.txt.gz --shard K --nshards S --out DIR

For every line "<d6> <edges> <flag> <twin>" with flag in {P, U, R} and shard = line_index mod S == K:
  edges <  22  -> 'never' (proved: legal stability forces B > Nn, i.e. edges > m(m-1)/2 = 21; notes §5)
  otherwise    -> curve1d.curve_test_fast; 'never' is a complete certificate that the full support is stable at
                  NO legal (eps, delta) (Lemma R), hence never an SFP.
Every non-'never' graph is written to DIR/curve_hits_K.jsonl with an exact alpha test (omega(.;2) Hurwitz), to be
sent to region2d. Counts go to DIR/curve_stats_K.json. A region2d SFP hit with alpha >= -1/2 is a finite-eps-only
counterexample candidate: STOP and hand-verify.
Resumable: if DIR/curve_done_K exists the shard is skipped.
"""
import sys, os, json, gzip, time, argparse
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/home/user/CLTN')
import ctln
from curve1d import curve_test_fast
from certify_M_violation import from_d6
from certify_graph_level_pure import omega_at

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True)
ap.add_argument('--shard', type=int, required=True)
ap.add_argument('--nshards', type=int, required=True)
ap.add_argument('--out', required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
done = os.path.join(a.out, f'curve_done_{a.shard}')
if os.path.exists(done):
    print('shard already done'); sys.exit(0)
op = gzip.open if a.input.endswith('.gz') else open
stats = {"lines": 0, "edges_lt_22": 0, "never": 0, "stable_somewhere": 0, "escalate": 0, "alpha_lt_minus_half": 0,
         "alpha_ge_minus_half": 0}
t0 = time.time()
with op(a.input, 'rt') as fh, open(os.path.join(a.out, f'curve_hits_{a.shard}.jsonl'), 'w') as hits:
    for i, line in enumerate(fh):
        if i % a.nshards != a.shard:
            continue
        d6, edges = line.split()[:2]
        stats["lines"] += 1
        if int(edges) < 22:
            stats["edges_lt_22"] += 1
            continue
        A = from_d6(d6)
        r = curve_test_fast(A)
        v = r['verdict']
        stats[v] = stats.get(v, 0) + 1
        if v != 'never':
            rec = {"d6": d6, "verdict": v, "stable_cells_t": [(c['lo'], c['hi'], c['sample_t']) for c in r.get('cells', []) if c['stable']]}
            if v == 'stable_somewhere':
                hs = ctln.hurwitz_stable(omega_at(A, F(2)))[0]
                rec["alpha_lt_minus_half"] = hs
                stats["alpha_lt_minus_half" if hs is True else "alpha_ge_minus_half"] += 1
            hits.write(json.dumps(rec) + "\n"); hits.flush()
        if stats["lines"] % 50000 == 0:
            print(a.shard, stats, round(time.time() - t0), flush=True)
stats["seconds"] = round(time.time() - t0)
json.dump(stats, open(os.path.join(a.out, f'curve_stats_{a.shard}.json'), 'w'), indent=1)
open(done, 'w').write('ok\n')
print('DONE', a.shard, stats, flush=True)
