"""Per-candidate exact structural features of all 9,133 n = 8 limit counterexamples.
Output: features.jsonl (one record per line of n8/n8_limit_HITS.jsonl, same order).
Run:  nice -n 10 python3 features.py   (2 worker processes)
"""
import sys, json, os
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

HERE = os.path.dirname(os.path.abspath(__file__))


def gfeat(A):
    """Graph-only exact features (any n)."""
    d6 = lib.to_d6(A)
    n = len(A)
    D = lib.nonedge(A)
    b, s, z = lib.pair_types(A)
    e = lib.arcs(A)
    assert 2 * b + s == e and b + s + z == n * (n - 1) // 2
    Dsinks = [v for v in range(n) if lib.outdeg(D, v) == 0]      # A: v sends to all others
    Dsources = [v for v in range(n) if lib.indeg(D, v) == 0]     # A: v receives from all others
    Diso = [v for v in Dsinks if v in Dsources]
    Asources = [v for v in range(n) if lib.indeg(A, v) == 0]
    Asinks = [v for v in range(n) if lib.outdeg(A, v) == 0]
    und = lib.underlying_edges(D)          # pairs that are NOT bidirected in A
    shape = lib.undirected_shape(n, und)
    # oriented part of D (= reversed single arcs of A)
    Dor = [[D[i][j] if not D[j][i] else 0 for j in range(n)] for i in range(n)]
    omega, Kmax = lib.bidirected_clique_number(A)
    tw = lib.twin_pairs(A)
    ctc = lib.closed_twin_classes(A)
    # sink classes: group D-sinks by their in-set from the non-sink part (Lemma X1'/X2' classes)
    core = [v for v in range(n) if v not in Dsinks]
    sink_insets = {}
    for v in Dsinks:
        S = tuple(sorted(i for i in core if A[i][v]))
        sink_insets.setdefault(S, []).append(v)
    # iterated peeling of D-sinks
    alive = list(range(n)); layers = []
    while True:
        sub = lib.induced(A, alive)
        sk = [alive[i] for i in range(len(alive)) if all(sub[i][j] for j in range(len(alive)) if j != i)]
        if not sk or len(sk) == len(alive):
            layers.append(len(sk)) if sk else None
            break
        layers.append(len(sk))
        alive = [v for v in alive if v not in sk]
    # D restricted to the core (non-sink nodes): strongly connected?
    Dcore = lib.induced(D, core) if core else None
    Dcore_strong = lib.strongly_connected(Dcore) if core and len(core) > 1 else None
    # one-way arcs of A: do they contain a directed 3-cycle?
    one = [[1 if (A[i][j] and not A[j][i]) else 0 for j in range(n)] for i in range(n)]
    ow3 = sum(1 for i in range(n) for j in range(n) for k in range(n)
              if i < j and i < k and one[i][j] and one[j][k] and one[k][i])
    return dict(
        d6=d6, n=n, arcs=e, missing=n * (n - 1) - e, Dcore_strong=Dcore_strong, oneway_3cycles=ow3,
        oneway_has_cycle=lib.has_directed_cycle(one), bi=b, single=s, nonadj=z,
        outdeg=[lib.outdeg(A, i) for i in range(n)], indeg=[lib.indeg(A, i) for i in range(n)],
        Dsinks=Dsinks, Dsources=Dsources, Diso=Diso, Asources=Asources, Asinks=Asinks,
        D_has_cycle=lib.has_directed_cycle(D), Dor_has_cycle=lib.has_directed_cycle(Dor),
        D_shape=shape, D_strong=lib.strongly_connected(D), A_strong=lib.strongly_connected(A),
        omega_bi=omega, Kmax=Kmax,
        n_bi_in_twins=len(tw['bi_in']), n_closed_twins=len(tw['closed']), n_nonadj_in_twins=len(tw['nonadj_in']),
        n_bi_out_twins=len(tw['bi_out']), closed_twin_class_sizes=[len(c) for c in ctc],
        dominated=lib.dominated_nodes(A),
        sink_classes=sorted([[list(S), vs] for S, vs in sink_insets.items()], key=lambda x: -len(x[1])),
        core_size=len(core), core_d6=lib.to_d6(lib.induced(A, core)) if core else None,
        peel_layers=layers, peel_residual=len(alive),
        aut=lib.automorphism_count(A),
        Dclass=[list(x) for x in lib.D_class_key(A)[0:1]] + [[list(v) for v in lib.D_class_key(A)[1]]],
        charpoly_A=lib.charpoly_int(A), charpoly_D=lib.charpoly_int(D),
    )


def feat(rec):
    A = lib.from_d6(rec['d6'])
    r = gfeat(A)
    assert r['d6'] == rec['d6'] and r['arcs'] == rec['edges']
    st = rec['cert']['stable_intervals']
    assert len(st) == 1
    lo, hi = st[0]['lo'], st[0]['hi']
    r.update(
        lo_kind=lo['kind'], lo_beta=lo.get('beta'), lo_vanishing=lo.get('vanishing'),
        hi_kind=hi['kind'], hi_vanishing=hi.get('vanishing'), hi_beta_interval=hi.get('beta_interval'),
        hi_q_interval=hi.get('q_interval'), hi_q_minpoly=hi.get('q_minpoly_hi_to_lo'),
    )
    return r


if __name__ == '__main__':
    hits = lib.load_hits()
    with Pool(2) as p:
        out = p.map(feat, hits, chunksize=50)
    with open(os.path.join(HERE, 'features.jsonl'), 'w') as f:
        for r in out:
            f.write(json.dumps(r) + '\n')
    print('records', len(out))
