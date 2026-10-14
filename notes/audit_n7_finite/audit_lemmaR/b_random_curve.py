"""Audit of n7/curve_test.py on random 7-node digraphs (exact; no floats in any decision).
usage: python3 b_random_curve.py SHARD NSHARDS [NGRAPHS=3000]  -> b_random_shard{SHARD}.jsonl
Per graph:
  * coefficient polynomials c_k(t) from curve_test.coeff_polys: degree bound deg c_k <= k+1, and
    values at 3 random NON-integer rationals vs an independent sympy QQ charpoly of -L(t0) (L built here);
  * curve_test vs curve_test_fast verdicts; whether coeff_never decided it;
  * for every stable cell: ctln.analyze_support at beta = 2 + 1/t + 10^-6 (eps = 1/(t+1), delta = (beta-1) eps)
    must give stable == True (if not, retried at 10^-9, 10^-12 and recorded)."""
import sys, os, json, random, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp
from n7.curve_test import coeff_polys, curve_test, curve_test_fast, coeff_never, coeffs_at
import ctln

def gen(idx, seed=777):
    rng = random.Random(seed * 1000003 + idx)
    kind = idx % 3
    if kind == 0:
        p = rng.choice([0.2, 0.35, 0.5, 0.65, 0.8]); A = [[0 if i == j else int(rng.random() < p) for j in range(7)] for i in range(7)]
    elif kind == 1:
        p = rng.choice([0.85, 0.9, 0.95]); A = [[0 if i == j else int(rng.random() < p) for j in range(7)] for i in range(7)]
    else:
        k = rng.randint(1, 21); arcs = [(i, j) for i in range(7) for j in range(7) if i != j]
        A = [[0 if i == j else 1 for j in range(7)] for i in range(7)]
        for (i, j) in rng.sample(arcs, k): A[i][j] = 0
    return A

def indep_coeffs(A, t0):
    m = len(A)
    N = [[1 if (i != j and A[j][i] == 0) else 0 for j in range(m)] for i in range(m)]
    L = sp.Matrix(m, m, lambda i, j: sp.Rational(t0.numerator, t0.denominator) * (1 if i == j else 0)
                  + (2 * sp.Rational(t0.numerator, t0.denominator) + 1) * N[i][j] + sp.Rational(t0.numerator, t0.denominator) ** 2)
    x = sp.Symbol('x')
    cp = (-L).charpoly(x).all_coeffs()          # det(x I + L)
    return [F(int(sp.fraction(c)[0]), int(sp.fraction(c)[1])) for c in cp]

def witness(A, t, eta):
    beta = 2 + 1 / t + eta; eps = 1 / (t + 1); delta = (beta - 1) * eps
    assert ctln.legal(eps, delta)
    return ctln.analyze_support(A, list(range(len(A))), eps, delta)["stable"]

def main():
    shard, nsh = int(sys.argv[1]), int(sys.argv[2]); ng = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
    out = open(f'/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_lemmaR/b_random_shard{shard}.jsonl', 'w')
    rng = random.Random(99 + shard)
    for idx in range(shard, ng, nsh):
        A = gen(idx); t0 = time.time()
        P = coeff_polys(A)
        degok = all(len(P[k]) <= k + 2 for k in range(len(P)))
        valok = True
        for _ in range(3):
            tq = F(rng.randint(1, 10 ** 6), rng.randint(2, 10 ** 4))
            if tq.denominator == 1: tq += F(1, 3)
            if indep_coeffs(A, tq) != coeffs_at(P, tq): valok = False
        r = curve_test(A, P); rf = curve_test_fast(A); cn = coeff_never(A, P)
        wit = []
        if r["verdict"] == "stable_somewhere":
            for c in r["cells"]:
                if c["stable"]:
                    t = F(c["sample_t"]); res = {"t": str(t), "1e-6": witness(A, t, F(1, 10 ** 6))}
                    if res["1e-6"] is not True:
                        res["1e-9"] = witness(A, t, F(1, 10 ** 9)); res["1e-12"] = witness(A, t, F(1, 10 ** 12))
                    wit.append(res)
        rec = {"idx": idx, "A": ["".join(map(str, row)) for row in A], "edges": sum(map(sum, A)),
               "deg_ok": degok, "values_ok": valok, "verdict": r["verdict"], "fast_verdict": rf["verdict"], "coeff_never": cn,
               "n_critical": r["n_critical"], "cells": [(c["lo"], c["hi"], c["stable"]) for c in r["cells"]],
               "witness": wit, "sec": round(time.time() - t0, 3)}
        out.write(json.dumps(rec) + "\n"); out.flush()

main()
