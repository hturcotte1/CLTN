"""Rigorous worst-case magnitude bounds (Hadamard) for every __int128 intermediate in filter7p, all integer arithmetic."""
from math import isqrt, factorial
LIM = 2**127 - 1
def ceil_sqrt(x): r = isqrt(x); return r if r*r == x else r + 1
def had(colsq, k):
    # bound on any k x k minor: product of the k largest column norms (ceil of sqrt of each squared norm, squared bound via product of squares)
    cs = sorted(colsq, reverse=True)[:k]
    prod = 1
    for c in cs: prod *= c
    return ceil_sqrt(prod)          # sqrt of product of squared norms >= product of norms
def bareiss_max(n, p, q):
    # matrix: column i = q*ones (sq norm n q^2), other columns: diag q, off-diag p or 0 (sq norm <= q^2 + (n-1) p^2)
    colsq = [n*q*q] + [q*q + (n-1)*p*p]*(n-1)
    worst = 0
    for c in range(n - 1):
        Hc = had(colsq, c + 1)                # entries at start of step c are (c+1)-minors
        worst = max(worst, 2*Hc*Hc, had(colsq, c + 2))
    return worst
SAMPLES = [(1,1),(1,2),(1,4),(1,8),(2,1),(1,16),(4,1),(1,32),(8,1),(16,1),(48,1),(1,64),(3,8),(3,4)]
for n in (6, 7, 8):
    w = max(bareiss_max(n, 2*q+p, q) for p, q in SAMPLES)
    wd = max(bareiss_max(n, 2+k, 1) for k in range(n + 2))
    # Descartes path
    m = n
    fmax = max(had([n] + [1 + (n-1)*(2+k)**2]*(n-1), n) for k in range(m + 2))
    Dmax = [2**k * fmax for k in range(m + 1)]
    # unsigned Stirling numbers of the first kind |s(k,j)|
    S = [[0]*(m+2) for _ in range(m+2)]; S[0][0] = 1
    for k in range(1, m+1):
        for j in range(1, k+1): S[k][j] = S[k-1][j-1] + (k-1)*S[k-1][j]
    scale = [Dmax[k] * (factorial(m)//factorial(k)) for k in range(m + 1)]
    termmax = max(scale[k]*S[k][j] for k in range(m+1) for j in range(k+1))
    coefmax = [sum(scale[k]*S[k][j] for k in range(j, m+1)) for j in range(m+1)]
    horner = sum(coefmax[j]*(m+1)**j for j in range(m+1))
    rhs = factorial(m) * fmax
    allmax = max(w, wd, termmax, max(coefmax), horner, rhs)
    print(f"n={n}: samples Bareiss {w:.3e}  descartes Bareiss {wd:.3e}  m!*coef {max(coefmax):.3e}  horner {horner:.3e}  ratio-to-2^127 {allmax/LIM:.3e}  SAFE={allmax < LIM}")
# largest p at q=64 that is still safe for n=7, 8
for n in (7, 8):
    p = 129
    while bareiss_max(n, p + 1, 64) < LIM: p += 1
    print(f"n={n}: q=64 Bareiss-safe for p <= {p} (beta <= {p/64:.2f})")
