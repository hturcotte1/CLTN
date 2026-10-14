"""300 'never' graphs from b_random_shard*.jsonl (all 276 decided by the full curve test, i.e. coeff_never False,
then the densest coeff_never ones up to 300), 20 random EXACT legal points each, ctln.analyze_support 'stable'
field must never be True (positivity irrelevant). Points: 10 near the legality curve (eta in [1e-9, 1e-3]),
10 broad (eta in [1e-3, 10]); t log-uniform in [1e-2, 1e4]; all rationals. Output d_never_points.json"""
import sys, json, glob, random
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
D = '/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_lemmaR/'
recs = [json.loads(l) for f in sorted(glob.glob(D + 'b_random_shard*.jsonl')) for l in open(f)]
nev = [r for r in recs if r['verdict'] == 'never']
full = sorted([r for r in nev if not r['coeff_never']], key=lambda r: r['idx'])
rest = sorted([r for r in nev if r['coeff_never']], key=lambda r: (-r['edges'], r['idx']))
sel = (full + rest)[:300]
rng = random.Random(4242)
def rat_log(lo_exp, hi_exp):           # rational, log-uniform-ish in [10^lo, 10^hi]
    e = rng.randint(lo_exp * 1000, hi_exp * 1000)
    mant = F(rng.randint(1000, 9999), 1000)
    return mant * (F(10) ** (e // 1000)) / 10
n_pts = 0; stable_true = []; none_ct = 0; summary = []
for r in sel:
    A = [[int(ch) for ch in row] for row in r['A']]
    for k in range(20):
        t = rat_log(-2, 4)
        eta = rat_log(-9, -3) if k < 10 else rat_log(-3, 1)
        beta = 2 + 1 / t + eta; eps = 1 / (t + 1); delta = (beta - 1) * eps
        assert ctln.legal(eps, delta)
        a = ctln.analyze_support(A, list(range(7)), eps, delta)
        n_pts += 1
        if a['stable'] is True: stable_true.append({"idx": r['idx'], "eps": str(eps), "delta": str(delta)})
        if a['stable'] is None: none_ct += 1
out = {"graphs": len(sel), "from_full_test": min(len(full), 300), "points": n_pts, "stable_true": stable_true, "stable_none": none_ct,
       "edges_range": [min(r['edges'] for r in sel), max(r['edges'] for r in sel)]}
json.dump(out, open(D + 'd_never_points.json', 'w'), indent=1); print(json.dumps(out))
