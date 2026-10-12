"""
ctln.py -- exact-arithmetic verifier for CTLN fixed points.

NO FLOATING POINT APPEARS ANYWHERE IN A DECISION PATH.
Every comparison is on exact rationals (fractions.Fraction) or exact integers.

--------------------------------------------------------------------------
CONVENTIONS  (pin these down; drift here silently corrupts everything)
--------------------------------------------------------------------------
n           number of neurons, indexed 0..n-1
A           adjacency matrix, A[i][j] == 1  iff  there is an edge  i -> j
            (so row i lists who i SENDS to)
eps, delta  the two CTLN knobs, rationals with 0 < eps, 0 < delta
theta       constant external drive, rational > 0 (default 1)

Weight matrix:
    W[i][j] = 0            if i == j
    W[i][j] = -1 + eps     if A[j][i] == 1     (j sends to i -> WEAK inhibition)
    W[i][j] = -1 - delta   otherwise           (no j->i edge -> STRONG inhibition)

NOTE the transpose: W[i][j] is the influence OF j ON i, so it reads A[j][i].

Dynamics:
    dx_i/dt = -x_i + [ sum_j W[i][j] x_j + theta ]_+

Legal parameter range: q := delta*(1-eps)/eps  must satisfy  q > 1,
equivalently 0 < eps < delta/(1+delta).

--------------------------------------------------------------------------
FIXED POINT WITH SUPPORT sigma
--------------------------------------------------------------------------
Three conditions, all checked exactly:

  (1) ON-SUPPORT SOLVE   (I - W_sigma) x_sigma = theta * 1
      and x_sigma > 0 strictly, componentwise.

      Useful exact identity (verified in tests below):
          I - W_sigma = -delta*I + (1+delta)*J - (eps+delta)*A_sigma
      where J is all-ones and A_sigma is the adjacency restricted to sigma.
      Everything is a polynomial in (eps, delta) with integer coefficients.

  (2) OFF-SUPPORT        for every k not in sigma:
          sum_{j in sigma} W[k][j] * x_j + theta  <  0
      (strict, for nondegeneracy) so that x_k = 0 is consistent.

  (3) STABILITY          all eigenvalues of  M = -I + W_sigma  have Re < 0.
      Decided by the HURWITZ DETERMINANT criterion on the exact characteristic
      polynomial -- never by a numerical eigenvalue solver. Near the q -> 1
      boundary (where the whole open problem lives) a float solver cannot tell
      Re = -1e-12 from Re = +1e-12.

--------------------------------------------------------------------------
GRAPH NOTIONS
--------------------------------------------------------------------------
clique(sigma)       every distinct i,j in sigma have BOTH i->j and j->i
target(sigma)       a node k not in sigma with j->k for EVERY j in sigma
target_free(sigma)  sigma is a clique with no target

THE CONJECTURE UNDER TEST (Curto-Morrison Open Question 10):
    stable fixed point supports  ==  target-free cliques
Geneson (arXiv:2607.21396, Jul 2026) disproved this with a 6-node graph at
q ~ 1.16, and proved it holds for q >= n - 2 - (n-3)*eps/2. The open region
is the gap between q slightly above 1 and q around n-2.
"""

from fractions import Fraction as F
from itertools import combinations

# ---------------------------------------------------------------- linear algebra
# All routines take and return Fractions. No floats, ever.


def zeros(r, c):
    return [[F(0) for _ in range(c)] for _ in range(r)]


def eye(m):
    return [[F(1) if i == j else F(0) for j in range(m)] for i in range(m)]


def matmul(X, Y):
    n, k, m = len(X), len(Y), len(Y[0])
    out = zeros(n, m)
    for i in range(n):
        Xi = X[i]
        for t in range(k):
            x = Xi[t]
            if x:
                Yt = Y[t]
                outi = out[i]
                for j in range(m):
                    outi[j] += x * Yt[j]
    return out


def trace(X):
    return sum((X[i][i] for i in range(len(X))), F(0))


def solve_exact(Amat, b):
    """Solve Amat x = b over the rationals by Gaussian elimination with
    exact pivoting. Returns None if singular."""
    m = len(Amat)
    M = [row[:] + [b[i]] for i, row in enumerate(Amat)]
    for col in range(m):
        piv = None
        for r in range(col, m):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(m):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * c for a, c in zip(M[r], M[col])]
    return [M[i][m] for i in range(m)]


def charpoly(M):
    """Exact characteristic polynomial of M via Faddeev-LeVerrier.
    Returns [1, c1, c2, ..., cm] meaning  lam^m + c1 lam^(m-1) + ... + cm."""
    m = len(M)
    coeffs = [F(1)]
    Mk = eye(m)
    for k in range(1, m + 1):
        if k > 1:
            Mk = matmul(M, Mk)
            for i in range(m):
                Mk[i][i] += coeffs[-1]
        AM = matmul(M, Mk) if k > 1 else M
        ck = -trace(AM) / k
        coeffs.append(ck)
    return coeffs


def _charpoly_direct(M):
    """Independent char-poly via exact interpolation, used only to cross-check
    charpoly() in the test suite. det(lam I - M) evaluated at m+1 points."""
    m = len(M)
    pts = list(range(-(m // 2) - 1, -(m // 2) - 1 + m + 1))
    vals = []
    for t in pts:
        B = [[(F(t) if i == j else F(0)) - M[i][j] for j in range(m)] for i in range(m)]
        vals.append(det_exact(B))
    # Lagrange interpolate the degree-m monic polynomial, return coeff list.
    coeffs = [F(0)] * (m + 1)
    for idx, t in enumerate(pts):
        num = [F(1)]
        den = F(1)
        for jdx, s in enumerate(pts):
            if jdx == idx:
                continue
            num = _polymul(num, [F(1), F(-s)])
            den *= F(t - s)
        scale = vals[idx] / den
        for p in range(len(num)):
            coeffs[p + (m + 1 - len(num))] += num[p] * scale
    return coeffs


def _polymul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def det_exact(M):
    m = len(M)
    if m == 0:
        return F(1)
    M = [row[:] for row in M]
    det = F(1)
    for col in range(m):
        piv = None
        for r in range(col, m):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return F(0)
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
            det = -det
        det *= M[col][col]
        pv = M[col][col]
        for r in range(col + 1, m):
            if M[r][col] != 0:
                f = M[r][col] / pv
                M[r] = [a - f * c for a, c in zip(M[r], M[col])]
    return det


def hurwitz_stable(coeffs):
    """Exact Hurwitz test. coeffs = [a0, a1, ..., am] for
    a0 lam^m + a1 lam^(m-1) + ... + am, with a0 > 0.

    Returns (verdict, detail) where verdict is True / False / None.
    None means a leading principal minor vanished -- the case sits exactly on
    the stability boundary. NEVER silently treat None as stable; such cases are
    logged and escalated for symbolic treatment.
    """
    m = len(coeffs) - 1
    if m == 0:
        return True, "degree 0"
    a = coeffs
    if a[0] <= 0:
        return False, "leading coefficient not positive"
    # Necessary condition: all coefficients strictly positive.
    if any(c <= 0 for c in a):
        if any(c == 0 for c in a):
            return None, "zero coefficient -- boundary case"
        return False, "negative coefficient"

    def aidx(k):
        return a[k] if 0 <= k <= m else F(0)

    H = [[aidx(2 * (i + 1) - (j + 1)) for j in range(m)] for i in range(m)]
    for k in range(1, m + 1):
        minor = det_exact([row[:k] for row in H[:k]])
        if minor == 0:
            return None, f"Hurwitz minor {k} vanishes -- boundary case"
        if minor < 0:
            return False, f"Hurwitz minor {k} negative"
    return True, "all Hurwitz minors positive"


# ---------------------------------------------------------------- CTLN core


def q_of(eps, delta):
    """The derived parameter q = delta(1-eps)/eps. Legal CTLN requires q > 1."""
    return delta * (1 - eps) / eps


def legal(eps, delta):
    return eps > 0 and delta > 0 and q_of(eps, delta) > 1


def weight_matrix(A, eps, delta):
    """W[i][j] = influence of j on i. Reads A[j][i]."""
    n = len(A)
    W = zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i == j:
                W[i][j] = F(0)
            elif A[j][i] == 1:
                W[i][j] = F(-1) + eps
            else:
                W[i][j] = F(-1) - delta
    return W


def submatrix(W, sigma):
    return [[W[i][j] for j in sigma] for i in sigma]


def identity_check_IminusW(A, sigma, eps, delta):
    """Cross-check the closed form  I - W_sigma = -d I + (1+d) J - (e+d) A_sigma
    against the directly constructed matrix. Returns True iff they agree."""
    W = weight_matrix(A, eps, delta)
    Ws = submatrix(W, sigma)
    m = len(sigma)
    direct = [[(F(1) if i == j else F(0)) - Ws[i][j] for j in range(m)] for i in range(m)]
    closed = zeros(m, m)
    for i, si in enumerate(sigma):
        for j, sj in enumerate(sigma):
            # A_sigma[i][j] uses the same transposed reading as W: edge sj -> si.
            a_ij = F(1) if (si != sj and A[sj][si] == 1) else F(0)
            kron = F(1) if i == j else F(0)
            closed[i][j] = -delta * kron + (1 + delta) * F(1) - (eps + delta) * a_ij
    return direct == closed


def is_clique(A, sigma):
    for i, j in combinations(sigma, 2):
        if A[i][j] != 1 or A[j][i] != 1:
            return False
    return True


def targets(A, sigma):
    """Nodes outside sigma that receive an edge from every member of sigma."""
    n = len(A)
    out = []
    for k in range(n):
        if k in sigma:
            continue
        if all(A[j][k] == 1 for j in sigma):
            out.append(k)
    return out


def is_target_free_clique(A, sigma):
    return is_clique(A, sigma) and not targets(A, sigma)


def analyze_support(A, sigma, eps, delta, theta=F(1)):
    """Full exact analysis of one (graph, support, parameter) case.

    Returns a dict -- the certificate record. Keys:
        is_fixed_point   True/False
        x                the on-support solution (list of Fractions) or None
        positive         bool
        off_support_ok   bool
        stable           True/False/None   (None == exact boundary case)
        stability_detail str
        charpoly         exact coefficients
        clique           bool
        targets          list
        target_free      bool
        predicted        what the conjecture says (stable iff target-free clique)
        agrees           whether reality matches the conjecture
    """
    sigma = sorted(sigma)
    n = len(A)
    W = weight_matrix(A, eps, delta)
    m = len(sigma)
    Ws = submatrix(W, sigma)

    IW = [[(F(1) if i == j else F(0)) - Ws[i][j] for j in range(m)] for i in range(m)]
    x = solve_exact(IW, [theta] * m)

    rec = {
        "sigma": sigma,
        "eps": str(eps),
        "delta": str(delta),
        "q": str(q_of(eps, delta)),
        "theta": str(theta),
        "clique": is_clique(A, sigma),
        "targets": targets(A, sigma),
    }
    rec["target_free"] = rec["clique"] and not rec["targets"]

    if x is None:
        # det(I - W_sigma) == 0 at this exact parameter value. The fixed-point
        # equation is degenerate here. Treat as UNDECIDED and escalate -- do not
        # silently classify it as "not a fixed point", which would hide any
        # counterexample sitting on this locus.
        rec.update(is_fixed_point=False, x=None, positive=False,
                   off_support_ok=False, stable=None,
                   stability_detail="singular I - W_sigma -- degenerate parameter",
                   charpoly=None, stable_fixed_point=False)
        rec["predicted"] = rec["target_free"]
        rec["BOUNDARY"] = True
        rec["agrees"] = None
        return rec

    positive = all(v > 0 for v in x)

    off_ok = True
    for k in range(n):
        if k in sigma:
            continue
        drive = theta + sum((W[k][sigma[j]] * x[j] for j in range(m)), F(0))
        if drive >= 0:
            off_ok = False
            break

    M = [[(F(-1) if i == j else F(0)) + Ws[i][j] for j in range(m)] for i in range(m)]
    cp = charpoly(M)
    stable, detail = hurwitz_stable(cp)

    is_fp = positive and off_ok
    actually_stable = bool(is_fp and stable is True)

    rec.update(
        is_fixed_point=is_fp,
        x=[str(v) for v in x],
        positive=positive,
        off_support_ok=off_ok,
        stable=stable,
        stability_detail=detail,
        charpoly=[str(c) for c in cp],
        stable_fixed_point=actually_stable,
    )
    rec["predicted"] = rec["target_free"]
    rec["BOUNDARY"] = (stable is None)
    # THREE-VALUED on purpose. An undecided boundary case is NOT a counterexample.
    # Collapsing None into False here would make the search chase phantoms.
    if stable is None:
        rec["agrees"] = None          # undecided; escalate to symbolic treatment
    else:
        rec["agrees"] = (actually_stable == rec["target_free"])
    return rec


def scan_supports(A, eps, delta, theta=F(1), min_size=1, max_size=None):
    """Every nonempty support, returning only the stable fixed points found
    plus any boundary cases that need escalation."""
    n = len(A)
    max_size = max_size or n
    stable_fps, boundary, violations = [], [], []
    for size in range(min_size, max_size + 1):
        for sigma in combinations(range(n), size):
            rec = analyze_support(A, list(sigma), eps, delta, theta)
            if rec.get("BOUNDARY"):
                boundary.append(rec)
            if rec.get("stable_fixed_point"):
                stable_fps.append(rec)
            if rec["agrees"] is False:      # strictly False, never None
                violations.append(rec)
    return {"stable_fixed_points": stable_fps,
            "boundary_cases": boundary,
            "conjecture_violations": violations}
