"""Exact Python cross-check of filter8 output lines (n = 8).

For each selected line:
  (1) independent big-int F2: m*chi_perp = chi_A - chi_{A+J} (fastregion.icharpoly),
      g(mu) = 2^d * p(mu - 1/2) via Fractions, Hurwitz minors via fastregion.idet;
      verdict S / B / N must match the C flag (an 'O' line may be anything).
  (2) fastregion.region(A, range(8), INF): limit-stage verdict.
Usage: python crosscheck8.py <filter8 output> <n_random_S> [seed]
"""
import sys, random, json, time
from fractions import Fraction as F
sys.path.insert(0, "/home/user/CLTN"); sys.path.insert(0, "/home/user/CLTN/n7")
import fastregion
from curve_test import from_d6
from math import comb

def f2_python(A):
    m = len(A)
    cA = fastregion.icharpoly(A)
    cAJ = fastregion.icharpoly([[a + 1 for a in row] for row in A])
    p = [cA[k + 1] - cAJ[k + 1] for k in range(m)]          # degree d = m-1, leading first
    d = m - 1
    assert p[0] == m
    # g(mu) = 2^d p(mu - 1/2) exactly, via Fractions then integer check
    g = [F(0)] * (d + 1)                                     # g[idx] = coeff of mu^(d-idx)
    for k in range(d + 1):
        e = d - k
        for j in range(e + 1):
            g[d - j] += F(p[k]) * comb(e, j) * F(-1, 2) ** (e - j) * 2 ** d
    assert all(x.denominator == 1 for x in g)
    g = [int(x) for x in g]
    if g[0] <= 0: return "N", g
    if any(x < 0 for x in g[1:]): return "N", g
    if any(x == 0 for x in g[1:]): return "B", g
    for k in range(1, d + 1):
        H = [[g[2 * (i + 1) - (j + 1)] if 0 <= 2 * (i + 1) - (j + 1) <= d else 0 for j in range(k)] for i in range(k)]
        mnr = fastregion.idet(H)
        if mnr < 0: return "N", g
        if mnr == 0: return "B", g
    return "S", g

def main():
    out, nS = sys.argv[1], int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    lines = list({tuple(l.split()) for l in open(out) if l.strip()}); lines.sort()
    S = [l for l in lines if l[2] == "S"]; O = [l for l in lines if l[2] == "O"]; B = [l for l in lines if l[2] == "B"]
    random.seed(seed)
    sel = random.sample(S, min(nS, len(S))) + O + random.sample(B, min(20, len(B)))
    print(f"lines {len(lines)}: S {len(S)} B {len(B)} O {len(O)}; checking {len(sel)}", flush=True)
    bad = 0; res = []; maxg = 0
    for d6, edges, flag in sel:
        A = from_d6(d6); assert sum(map(sum, A)) == int(edges)
        v, g = f2_python(A); maxg = max(maxg, max(abs(x) for x in g))
        t0 = time.time(); cert = fastregion.region(A, list(range(len(A))), fastregion.INF); dt = time.time() - t0
        ok = (flag == "O") or (v == flag)
        if not ok: bad += 1
        rec = {"d6": d6, "edges": int(edges), "c_flag": flag, "py_F2": v, "g": g,
               "region_verdict": cert["verdict"], "region_sec": round(dt, 2)}
        if cert["verdict"] == "stable_somewhere":
            rec["stable_intervals"] = cert["stable_intervals"]
        res.append(rec)
        print(f"{d6} {edges} C={flag} py={v} region={cert['verdict']} ({dt:.1f}s) {'OK' if ok else 'MISMATCH'}", flush=True)
    print(f"mismatches {bad}; max |g| coefficient {maxg} (~2^{maxg.bit_length()})")
    json.dump(res, open(out + ".crosscheck.json", "w"), indent=1)

main()
