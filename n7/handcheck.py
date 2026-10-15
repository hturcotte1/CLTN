"""
handcheck.py -- print everything a human needs to verify one candidate by hand, exactly.

   python3 n7/handcheck.py '<digraph6>' <eps> <delta>          e.g. '&FNFL~[udo?' 1/49 33/1568

Prints: the adjacency rows (row i = out-neighbours of i), the parameters (t, beta, q, legality),
the matrix I - W with every entry as a fraction (W[i][j] = -1+eps if A[j][i] = 1, else -1-delta,
0 on the diagonal; CLAUDE.md section 2), the solution x of (I - W) x = 1 (exact), the characteristic
polynomial of -I + W (exact, via ctln), and all Hurwitz minors (exact), plus whether the support is a
clique and whether the graph is domination-free. Independent of the certificate files. No floats.
"""
import sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ctln import weight_matrix, analyze_support, is_clique, charpoly as ctln_charpoly
from graphs import is_domination_free
from n7.curve_test import from_d6


def hurwitz_minors(coeffs):
    """coeffs highest degree first, monic; returns the leading principal minors of the Hurwitz matrix."""
    d = len(coeffs) - 1
    H = [[(coeffs[2 * (i + 1) - (j + 1)] if 0 <= 2 * (i + 1) - (j + 1) <= d else F(0)) for j in range(d)] for i in range(d)]
    out = []
    for k in range(1, d + 1):
        M = [row[:k] for row in H[:k]]
        det = F(1)
        for c in range(k):
            p = next((r for r in range(c, k) if M[r][c] != 0), None)
            if p is None:
                det = F(0); break
            if p != c:
                M[c], M[p] = M[p], M[c]; det = -det
            det *= M[c][c]
            for r in range(c + 1, k):
                f = M[r][c] / M[c][c]
                for cc in range(c, k):
                    M[r][cc] -= f * M[c][cc]
        out.append(det)
    return out


def main():
    d6, eps, delta = sys.argv[1], F(sys.argv[2]), F(sys.argv[3])
    A = from_d6(d6); n = len(A); V = list(range(n))
    t = (1 - eps) / eps; beta = 1 + delta / eps; q = delta * (1 - eps) / eps
    print(f"graph {d6}: {n} nodes, {sum(map(sum, A))} arcs; clique: {is_clique(A, V)}; domination-free: {is_domination_free(A)}")
    for i, row in enumerate(A):
        print(f"  node {i} -> {[j for j in range(n) if row[j]]}")
    print(f"eps = {eps}, delta = {delta}, t = (1-eps)/eps = {t}, beta = 1+delta/eps = {beta}, q = delta(1-eps)/eps = {q}, legal (q > 1): {q > 1}")
    W = weight_matrix(A, eps, delta)
    print("I - W (rows):")
    for i in range(n):
        print("  [" + ", ".join(str((1 if i == j else 0) - W[i][j]) for j in range(n)) + "]")
    r = analyze_support(A, V, eps, delta)
    print("solution x of (I - W) x = 1:", [str(v) for v in r["x"]] if r["x"] is not None else None)
    print("all x_i > 0:", r["x"] is not None and all(F(str(v)) > 0 for v in r["x"]))
    cp = [F(str(c)) for c in r["charpoly"]]
    print("characteristic polynomial of -I + W (highest degree first):", [str(c) for c in cp])
    H = hurwitz_minors(cp)
    print("Hurwitz minors Delta_1..Delta_%d:" % len(H), [str(h) for h in H])
    print("all minors > 0:", all(h > 0 for h in H), "| ctln verdict:", r["stable"], "-", r.get("stability_detail"))
    print("stable positive fixed point on the full support:", bool(r.get("stable_fixed_point")), "| target-free clique:", bool(r.get("target_free")))


if __name__ == "__main__":
    main()
