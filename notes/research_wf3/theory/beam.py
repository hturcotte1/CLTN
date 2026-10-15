"""HEURISTIC beam search (floats only steer; nothing here is a certificate).
level n -> n+1: all 2^(2n) one-node extensions of each beam graph; score q_lim (heur_extend.qlim on a beta grid);
dedupe by nauty labelg canonical digraph6; keep top K.  Output JSON per level.  Certification is done by certify_beam.py."""
import numpy as np, itertools, json, sys, subprocess, time
sys.path.insert(0, '.')
import heur_extend as he, tlib
LABELG = '/home/user/CLTN/n7/labelg'
def canon(d6s):
    p = subprocess.run([LABELG, '-q'], input='\n'.join(d6s) + '\n', capture_output=True, text=True)
    return p.stdout.split()
def ext_chunks(A, chunk=1 << 15):
    A = np.array(A, dtype=np.int8); m = len(A); n = m + 1
    total = 1 << (2 * m)
    for start in range(0, total, chunk):
        idx = np.arange(start, min(total, start + chunk), dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(2 * m - 1, -1, -1)) & 1).astype(np.int8)
        B = np.zeros((len(idx), n, n), dtype=np.int8)
        B[:, :m, :m] = A
        B[:, :m, m] = bits[:, :m]; B[:, m, :m] = bits[:, m:]
        yield B
def run_level(beam, K, betas, keep_per_base=200):
    cand = []
    t0 = time.time()
    for bi, A in enumerate(beam):
        local = []
        for B in ext_chunks(A):
            ev = he.alphas(B)
            a = ev.real.max(axis=1)
            sel = np.where(a < -0.5 - 1e-9)[0]
            if len(sel) == 0: continue
            a2, q = he.qlim(B[sel], ev[sel], betas)
            for j in np.where(q > 0)[0]:
                local.append((float(q[j]), B[sel[j]].tolist(), float(a2[j])))
        local.sort(key=lambda x: -x[0])
        cand += local[:keep_per_base]
        print('  base %d/%d  local hits %d  %.0fs' % (bi + 1, len(beam), len(local), time.time() - t0), flush=True)
    cand.sort(key=lambda x: -x[0])
    d6s = [tlib.to_d6(c[1]) for c in cand]
    cs = canon(d6s) if d6s else []
    out, seen = [], set()
    for c, cf in zip(cand, cs):
        if cf in seen: continue
        seen.add(cf); out.append({'q_heur': c[0], 'd6': tlib.to_d6(c[1]), 'canon': cf, 'alpha_heur': c[2]})
        if len(out) >= K: break
    return out
if __name__ == '__main__':
    start = json.load(open(sys.argv[1]))      # list of d6
    levels = int(sys.argv[2]); K = int(sys.argv[3]); tag = sys.argv[4]
    betas = np.concatenate([np.linspace(2.0005, 2.6, 240), np.linspace(2.605, 4.5, 380)])
    beam = [tlib.from_d6(s) for s in start]
    for L in range(levels):
        n = len(beam[0]) + 1
        print('level n=%d, beam %d' % (n, len(beam)), flush=True)
        out = run_level(beam, K, betas)
        json.dump(out, open('beam_%s_n%d.json' % (tag, n), 'w'), indent=0)
        print('n=%d top q_heur:' % n, [round(o['q_heur'], 4) for o in out[:10]], flush=True)
        beam = [tlib.from_d6(o['d6']) for o in out]
