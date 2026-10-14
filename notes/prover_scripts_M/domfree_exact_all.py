"""EXACT (no floats in any decision). For every domination-free digraph on n nodes (n <= 6 from the repo key caches):
run the Lemma-R curve test (legal-stable somewhere, full support, pure stability) and, for each graph that is
legal-stable somewhere, decide exactly whether alpha(A|perp) < -1/2 (omega(.;2) Hurwitz, Lemma C'/I).
Graph-level statement tested:  domination-free & legal-stable somewhere  =>  alpha < -1/2."""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from enum_cache import key_from_int
from graphs import from_key, is_domination_free
from curve1d import curve_test
from certify_graph_level_pure import omega_at
from to_d6 import to_d6
n = int(sys.argv[1]); lo = int(sys.argv[2]); hi = int(sys.argv[3]); outp = sys.argv[4]
if n <= 4:
    from graphs import enumerate_digraphs
    keys = list(enumerate_digraphs(n)); As = [from_key(k, n) for k in keys]
else:
    As = [from_key(key_from_int(int(l), n), n) for l in open(f'/home/user/CLTN/certificates/digraphs_n{n}.keys')]
As = As[lo:hi]
tally = {}; t0 = time.time()
with open(outp, 'w') as fh:
    for idx, A in enumerate(As):
        if not is_domination_free(A):
            tally['dominated'] = tally.get('dominated', 0) + 1; continue
        r = curve_test(A)
        v = r['verdict']
        if v == 'stable_somewhere':
            hs = ctln.hurwitz_stable(omega_at(A, F(2)))[0]
            v = v + ('|alpha<-1/2' if hs is True else ('|alpha>=-1/2 !!' if hs is False else '|alpha boundary'))
            fh.write(json.dumps({"idx": lo + idx, "d6": to_d6(A), "verdict": v,
                                 "stable_cells_t": [(c['lo'], c['hi'], c['sample_t']) for c in r['cells'] if c['stable']]}) + "\n")
        elif v == 'escalate':
            fh.write(json.dumps({"idx": lo + idx, "d6": to_d6(A), "verdict": v}) + "\n")
        tally[v] = tally.get(v, 0) + 1
        if idx % 20000 == 0: print(lo + idx, tally, round(time.time() - t0), flush=True)
print('DONE', n, lo, hi, tally, round(time.time() - t0), flush=True)
