"""
validate_fastregion.py -- cross-validate fastregion.region against the
independent symbolic path symbolic.q_region (sympy symbolic matrices + real_roots)
and against ctln.analyze_support at every sample point.

For each (graph, support, eps): the stable q-intervals must agree exactly:
same number of intervals, and each endpoint of q_region (a sympy exact number)
must be a root of fastregion's minimal polynomial and lie in its isolating
interval. Any discrepancy is printed and counted. Run:
   python3 validate_fastregion.py --n 4            (all 218 graphs, all supports of size>=2)
   python3 validate_fastregion.py --n 5 --sample 300
"""
import argparse, json, random, sys, time
from fractions import Fraction as F
from itertools import combinations
from multiprocessing import Pool

import sympy as sp
from graphs import enumerate_digraphs, from_key
from fastregion import region, beta_data, INF
from symbolic import q_region

def check_case(args):
    A, sigma, eps = args
    t = (1 - eps) / eps
    try:
        c = region(A, sigma, t, crosscheck=True)
    except Exception as e:
        return ("fastregion_error", A, sigma, str(eps), repr(e))
    if c["verdict"] == "escalate":
        return ("escalate", A, sigma, str(eps), None)
    try:
        s = q_region(A, sigma, sp.Rational(eps.numerator, eps.denominator))
    except Exception as e:
        return ("symbolic_error", A, sigma, str(eps), repr(e))
    ours = c.get("stable_intervals", [])
    theirs = s["intervals"]
    if len(ours) != len(theirs):
        return ("count_mismatch", A, sigma, str(eps), (len(ours), len(theirs), [str(x) for x in theirs]))
    q = sp.Symbol("q")
    for o, (lo, hi) in zip(ours, theirs):
        for e, v in ((o["lo"], lo), (o["hi"], hi)):
            if e["kind"] == "beta_min":
                if v != 1: return ("endpoint_mismatch", A, sigma, str(eps), ("lo should be 1", str(v)))
            elif e["kind"] == "inf":
                if v is not None: return ("endpoint_mismatch", A, sigma, str(eps), ("hi should be inf", str(v)))
            else:
                mp = sp.Poly(e["q_minpoly_hi_to_lo"], q)
                val = sp.nsimplify(v)
                if sp.simplify(mp.as_expr().subs(q, val)) != 0:
                    return ("endpoint_mismatch", A, sigma, str(eps), ("not a root", str(v), e["q_minpoly_hi_to_lo"]))
                lo_i, hi_i = sp.Rational(e["q_interval"][0]), sp.Rational(e["q_interval"][1])
                if not (lo_i <= val <= hi_i):
                    return ("endpoint_mismatch", A, sigma, str(eps), ("outside interval", str(v), e["q_interval"]))
    return ("ok", A, sigma, str(eps), len(ours))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--sample", type=int, default=None)
    ap.add_argument("--eps", nargs="+", default=["1/2", "1/4", "1/10"])
    ap.add_argument("--procs", type=int, default=2)
    args = ap.parse_args()
    n = args.n
    keys = sorted(enumerate_digraphs(n))
    rng = random.Random(7)
    graphs = [from_key(k, n) for k in keys]
    if args.sample:
        graphs = rng.sample(graphs, args.sample)
    cases = []
    for A in graphs:
        sups = [list(range(n))] if args.sample else [list(s) for size in range(2, n + 1) for s in combinations(range(n), size)]
        for sigma in sups:
            for e in args.eps:
                cases.append((A, sigma, F(e)))
    print(f"{len(cases)} cases", flush=True)
    t0 = time.time()
    tally = {}
    bad = []
    with Pool(args.procs) as pool:
        for i, r in enumerate(pool.imap_unordered(check_case, cases, chunksize=4)):
            tally[r[0]] = tally.get(r[0], 0) + 1
            if r[0] != "ok":
                bad.append(r); print("  ", r, flush=True)
            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(cases)}  {tally}  {time.time()-t0:.0f}s", flush=True)
    print("DONE", tally, f"{time.time()-t0:.0f}s")
    json.dump({"n": n, "cases": len(cases), "tally": tally, "bad": [list(map(str, b)) for b in bad]},
              open(f"certificates/validate_fastregion_n{n}.json", "w"), indent=1)

if __name__ == "__main__":
    main()
