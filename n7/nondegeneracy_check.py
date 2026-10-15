"""
nondegeneracy_check.py -- exact nondegeneracy certificates at witness parameters.

Geneson's conjecture (and Curto-Morrison's) is stated for NONDEGENERATE CTLNs. At a parameter point
(eps, delta) we certify nondegeneracy by checking, for EVERY support sigma (all 2^n - 1 nonempty
subsets): det(I - W_sigma) != 0, the solution x_sigma of (I - W_sigma) x = 1 has no zero entry, and
for every k outside sigma the drive 1 + sum_j W[k][j] x_j is nonzero. All in Fractions.
Inputs: the three signed-off n = 6 graphs (certificates/COUNTEREXAMPLE_CANDIDATE_n6_*.json and
Geneson's graph) at their witness points, and the 124 n = 7 candidates at one rational point inside
their first stable cell (from certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl).
Output: certificates/nondegeneracy_witnesses.json. No floating point anywhere.
"""
import sys, os, json, itertools
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ctln import weight_matrix, analyze_support
from n7.curve_test import from_d6
from fastregion import idet


def det_frac(M):
    n = len(M); M = [row[:] for row in M]; d = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                for k in range(c, n):
                    M[r][k] -= f * M[c][k]
    return d


def solve(M, b):
    n = len(M); A = [row[:] + [b[i]] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]
                for k in range(c, n + 1):
                    A[r][k] -= f * A[c][k]
    return [A[i][n] / A[i][i] for i in range(n)]


def nondegenerate(A, eps, delta):
    n = len(A); W = weight_matrix(A, eps, delta)
    zero_det = []; zero_x = []; zero_drive = []
    for m in range(1, n + 1):
        for sigma in itertools.combinations(range(n), m):
            M = [[(1 if i == j else 0) - W[i][j] for j in sigma] for i in sigma]
            d = det_frac(M)
            if d == 0:
                zero_det.append(sigma); continue
            x = solve(M, [F(1)] * m)
            if any(v == 0 for v in x):
                zero_x.append(sigma)
            for k in range(n):
                if k in sigma:
                    continue
                if 1 + sum(W[k][j] * x[i] for i, j in enumerate(sigma)) == 0:
                    zero_drive.append((sigma, k))
    return {"nondegenerate": not (zero_det or zero_x or zero_drive), "singular_supports": zero_det,
            "zero_coordinate_supports": zero_x, "zero_drive": zero_drive, "supports_checked": 2 ** n - 1}


def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(base)
    out = {"definition": "det(I - W_sigma) != 0, x_sigma has no zero entry, no off-support drive is exactly 0, for all nonempty sigma; exact Fractions", "graphs": []}
    # n = 6 signed-off graphs: Geneson (paper parameters) and the two candidate files
    pts6 = [("Geneson &EZjtwWE (paper parameters)", from_d6(open("n7/geneson.d6").read().split()[0]), F(1, 100), F(29, 2500))]
    for f in ["certificates/COUNTEREXAMPLE_CANDIDATE_n6_idx1387528.json", "certificates/COUNTEREXAMPLE_CANDIDATE_n6_idx1424629.json"]:
        c = json.load(open(f))
        A = c["A_ours"]; w0 = c["pointwise_witnesses_ctln_analyze_support"][0]
        pts6.append((os.path.basename(f), A, F(w0["eps"]), F(w0["delta"])))
    for name, A, eps, delta in pts6:
        r = nondegenerate(A, eps, delta); r.update(graph=name, eps=str(eps), delta=str(delta))
        s = analyze_support(A, list(range(len(A))), eps, delta); r["full_support_stable"] = s["stable"]
        out["graphs"].append(r); print(name, eps, delta, "nondegenerate:", r["nondegenerate"], "stable:", s["stable"], flush=True)
    # n = 7: one rational point inside the first stable cell of each hit
    bad = 0
    for line in open("certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl"):
        h = json.loads(line); A = from_d6(h["d6"])
        lo, hi, bs, tlo, thi = h["stable_cells"][0]
        beta = F(bs); t = F(tlo) + 1 if thi is None else (F(tlo) + F(thi)) / 2
        eps = 1 / (t + 1); delta = (beta - 1) * eps
        s = analyze_support(A, list(range(7)), eps, delta)
        r = nondegenerate(A, eps, delta); r.update(graph=h["d6"], eps=str(eps), delta=str(delta), full_support_stable=s["stable"], full_support_positive=bool(s["x"] and all(F(str(v)) > 0 for v in s["x"])))
        out["graphs"].append(r)
        if not (r["nondegenerate"] and s["stable"] is True): bad += 1; print("ATTENTION", h["d6"], r["nondegenerate"], s["stable"], flush=True)
    out["n7_count"] = len(out["graphs"]) - len(pts6); out["n7_all_nondegenerate_and_stable"] = bad == 0
    json.dump(out, open("certificates/nondegeneracy_witnesses.json", "w"), indent=1)
    print(json.dumps({"n6": [(g["graph"], g["nondegenerate"]) for g in out["graphs"][:3]], "n7_checked": out["n7_count"], "n7_all_ok": bad == 0}))


if __name__ == "__main__":
    main()
