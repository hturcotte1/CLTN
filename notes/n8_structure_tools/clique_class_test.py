"""Exhaustive test of the bidirected-clique bound omega_bi <= n - 3 on the class of digraphs that contain a
bidirected clique on n - 2 vertices, for given n (validation at n = 6, 7, 8 against the complete maps; the
new test is n = 9).  Every such digraph is: K = {0..n-3} bidirected clique, x = n-2, y = n-1; each clique node
has one of 16 patterns (i->x, x->i, i->y, y->i); the pair (x, y) has one of 4 relations.  Up to relabelling
inside K the graph is a multiset of n - 2 patterns plus the (x, y) relation, so the enumeration is complete.
Pipeline (decisions exact): graphs.is_domination_free -> strict spectral test (all roots of
psi_D(z) = chi_D - chi_{D+J} in Re z > -1/2, ctln.hurwitz_stable on psi_D(-x - 1/2); a vanishing minor means
alpha >= -1/2, which is never limit-stable by Lemma I.3(c2)) -> nauty dedup -> fastregion.region(A, V, INF,
replay=True), the production limit classifier.
Usage: nice -n 10 python3 clique_class_test.py n      Output: clique_class_test_n{n}.json"""
import sys, os, json, time
from itertools import combinations_with_replacement
from fractions import Fraction as F
from multiprocessing import Pool
from collections import Counter
sys.path.insert(0, '/home/user/CLTN')
import lib
from graphs import is_domination_free, is_complete_bidirected
from ctln import hurwitz_stable

HERE = os.path.dirname(os.path.abspath(__file__))
PAT = [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]
XY = [(0, 0), (1, 0), (0, 1), (1, 1)]   # (x->y, y->x)


def build(n, pats, xy):
    A = [[0] * n for _ in range(n)]
    x, y = n - 2, n - 1
    for i in range(n - 2):
        for j in range(n - 2):
            if i != j: A[i][j] = 1
        a, b, c, d = PAT[pats[i]]
        A[i][x], A[x][i], A[i][y], A[y][i] = a, b, c, d
    A[x][y], A[y][x] = xy
    return A


def compose_shift(p):
    """coefficients (hi->lo) of p(-x - 1/2), made monic-positive."""
    deg = len(p) - 1
    # Horner in Fractions on polynomials in x
    res = [F(0)]
    lin = [F(-1), F(-1, 2)]              # -x - 1/2
    for c in p:
        # res = res * lin + c
        out = [F(0)] * (len(res) + 1)
        for i, r in enumerate(res):
            out[i] += r * lin[0]; out[i + 1] += r * lin[1]
        out[-1] += c
        res = out
    while len(res) > 1 and res[0] == 0:
        res = res[1:]
    if res[0] < 0:
        res = [-r for r in res]
    return res


def strict(A):
    n = len(A)
    D = lib.nonedge(A)
    DJ = [[D[i][j] + 1 for j in range(n)] for i in range(n)]
    c1, c2 = lib.charpoly_int(D), lib.charpoly_int(DJ)
    psi = [a - b for a, b in zip(c1, c2)][1:]
    assert psi[0] == n
    v, _ = hurwitz_stable(compose_shift(psi))
    return v is True


def job(args):
    n, pats, xy = args
    A = build(n, pats, xy)
    if not is_domination_free(A):
        return ('dominated', None)
    if not strict(A):
        return ('spectral', None)
    return ('strict', lib.to_d6(A))


def main(n):
    t0 = time.time()
    tasks = [(n, pats, xy) for pats in combinations_with_replacement(range(16), n - 2) for xy in XY]
    with Pool(2) as p:
        res = p.map(job, tasks, chunksize=2000)
    c = Counter(r[0] for r in res)
    surv = sorted({r[1] for r in res if r[0] == 'strict'})
    canon = lib.labelg_canon(surv) if surv else []
    uniq = {}
    for d, cn in zip(surv, canon):
        uniq.setdefault(cn, d)
    from fastregion import region, beta_data, INF
    verdicts = Counter(); hits = []
    for cn, d in sorted(uniq.items()):
        A = lib.from_d6(d)
        r = region(A, list(range(n)), INF, bd=beta_data(A, list(range(n))), replay=True)
        comp = is_complete_bidirected(A)
        key = r['verdict'] + (' (complete)' if comp else '')
        verdicts[key] += 1
        if r['verdict'] != 'never' and not comp:
            hits.append(dict(d6=d, verdict=r['verdict']))
    out = dict(n=n, configurations=len(tasks), stage_counts=dict(c), strict_survivors_labelled=len(surv),
               strict_survivors_up_to_isomorphism=len(uniq), limit_verdicts=dict(verdicts),
               noncomplete_limit_stable_or_undecided=hits, seconds=round(time.time() - t0, 1))
    json.dump(out, open(os.path.join(HERE, f'clique_class_test_n{n}.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(int(sys.argv[1]))
