"""
symbolic.py -- exact parameter-region analysis for CTLN supports.

WHY THIS MODULE IS THE CENTRE OF THE PROJECT
--------------------------------------------
A sweep of 1,495 (graph, q) pairs over mixed n=6 graphs at sampled rational q
found ZERO decided counterexamples but 84 UNDECIDED boundary cases. The
counterexamples live in the undecided pile. Sampling q at points is the wrong
instrument: a case that is stable on a thin q-interval will be missed by every
sample that misses the interval, and a case sitting exactly on a Hurwitz
boundary at the sampled point tells you nothing about its neighbourhood.

So: keep the parameters SYMBOLIC and solve for the exact region where each
condition holds. That is what turns "we found another counterexample" into
"we characterized the gap", which is the result worth having.

TWO MODES
---------
full 2-parameter   conditions_symbolic(A, sigma)
                   returns polynomial conditions in (eps, delta). Exact and
                   complete, but determinants of 6x6 matrices in two symbols
                   get heavy. Use for the cases that matter.

1-parameter slice  q_region(A, sigma, eps_value)
                   fixes eps to an exact rational and keeps q symbolic via
                        delta = q*eps/(1-eps)
                   This is the practical workhorse: still fully exact, far
                   faster, and it answers the question the conjecture is
                   actually phrased in ("for which q does this flip?").
                   Run it at several eps values to reconstruct the 2-D picture.

CONVENTIONS are identical to ctln.py. Do not re-derive them here.
"""

import sympy as sp
from sympy import Rational, Symbol, Poly, QQ

eps_s = Symbol("eps", positive=True)
delta_s = Symbol("delta", positive=True)
q_s = Symbol("q", positive=True)
theta_s = Rational(1)


# ---------------------------------------------------------------- construction


def W_symbolic(A, eps=eps_s, delta=delta_s):
    """W[i][j] = influence of j on i; reads A[j][i]. Same transpose as ctln.py."""
    n = len(A)
    return sp.Matrix(n, n, lambda i, j:
                     sp.Integer(0) if i == j
                     else (-1 + eps if A[j][i] == 1 else -1 - delta))


def I_minus_W_sigma(A, sigma, eps=eps_s, delta=delta_s):
    """Built two ways and cross-checked against the closed form
          I - W_s = -delta*I + (1+delta)*J - (eps+delta)*A_s
    Raises if they disagree -- that would mean a convention drift."""
    W = W_symbolic(A, eps, delta)
    m = len(sigma)
    direct = sp.eye(m) - sp.Matrix(m, m, lambda i, j: W[sigma[i], sigma[j]])
    As = sp.Matrix(m, m, lambda i, j:
                   sp.Integer(1) if (i != j and A[sigma[j]][sigma[i]] == 1)
                   else sp.Integer(0))
    closed = -delta * sp.eye(m) + (1 + delta) * sp.ones(m, m) - (eps + delta) * As
    if sp.simplify(direct - closed) != sp.zeros(m, m):
        raise AssertionError("closed form disagrees with direct construction "
                             "-- convention drift, STOP")
    return direct


# ---------------------------------------------------------------- conditions


def hurwitz_minors(poly_coeffs):
    """Leading principal minors of the Hurwitz matrix, as symbolic expressions.
    coeffs = [a0, a1, ..., am] for a0*lam^m + ... + am, a0 > 0.
    All roots have Re < 0  iff  every minor is > 0."""
    m = len(poly_coeffs) - 1
    a = poly_coeffs

    def ai(k):
        return a[k] if 0 <= k <= m else sp.Integer(0)

    H = sp.Matrix(m, m, lambda i, j: ai(2 * (i + 1) - (j + 1)))
    return [sp.expand(H[:k, :k].det()) for k in range(1, m + 1)]


def conditions_symbolic(A, sigma, eps=eps_s, delta=delta_s):
    """All conditions for 'sigma supports a stable fixed point', symbolically.

    Returns dict with keys:
      positivity   list of expressions that must each be > 0  (x_i > 0)
      off_support  list of expressions that must each be < 0  (drive at k)
      stability    list of expressions that must each be > 0  (Hurwitz minors)
      x            the symbolic fixed point
      denom        the common denominator det(I - W_sigma); its sign matters
    """
    sigma = sorted(sigma)
    n, m = len(A), len(sigma)
    IW = I_minus_W_sigma(A, sigma, eps, delta)
    detIW = sp.expand(IW.det())
    adj = IW.adjugate()
    # x = (I-W)^-1 theta*1 = adj * theta*1 / det. Keep numerators separate so
    # the sign conditions stay polynomial rather than rational.
    num = [sp.expand(sum(adj[i, j] for j in range(m)) * theta_s) for i in range(m)]

    W = W_symbolic(A, eps, delta)
    off = []
    for k in range(n):
        if k in sigma:
            continue
        # drive_k = theta + sum_j W[k][sigma_j] * x_j ; multiply through by det
        expr = theta_s * detIW + sum(W[k, sigma[j]] * num[j] for j in range(m))
        off.append(sp.expand(expr))

    M = sp.Matrix(m, m, lambda i, j:
                  (-1 if i == j else 0) + W[sigma[i], sigma[j]])
    lam = Symbol("lam")
    cp = sp.Poly(sp.expand((lam * sp.eye(m) - M).det()), lam)
    coeffs = cp.all_coeffs()

    return {
        "sigma": sigma,
        "x_numerators": num,
        "denom": detIW,
        "positivity": num,          # each num[i]/detIW > 0
        "off_support": off,         # each off[k]/detIW < 0
        "stability": hurwitz_minors(coeffs),
        "charpoly": coeffs,
    }


# ---------------------------------------------------------------- q slice


def q_region(A, sigma, eps_value, verbose=False):
    """Fix eps to an exact rational, keep q symbolic, and return the exact set
    of q > 1 for which sigma supports a stable fixed point.

    Substitution:  delta = q*eps/(1-eps)   (so that q = delta(1-eps)/eps)

    Returns dict:
      intervals   list of (lo, hi) exact endpoints where ALL conditions hold,
                  intersected with the legal region q > 1. Endpoints are
                  sympy Rationals or algebraic numbers.
      critical_q  sorted exact roots of every boundary polynomial in q > 1 --
                  these are the only places the verdict can change.
      conditions  the substituted polynomial conditions, for the certificate
    """
    e = Rational(eps_value)
    if not (0 < e < 1):
        raise ValueError("eps must lie strictly between 0 and 1")
    d = q_s * e / (1 - e)

    cond = conditions_symbolic(A, sigma, eps=e, delta=d)

    def as_poly(expr):
        return sp.Poly(sp.together(sp.expand(expr)).as_numer_denom()[0], q_s)

    den_poly = as_poly(cond["denom"])
    pos_polys = [as_poly(p) for p in cond["positivity"]]
    off_polys = [as_poly(p) for p in cond["off_support"]]
    stab_polys = [as_poly(p) for p in cond["stability"]]

    # Every sign change happens at a real root of one of these polynomials.
    # Roots are kept EXACT (radicals or CRootOf). No nsimplify: it is a
    # numerical heuristic and produced nonsense closed forms (ANOMALIES.md
    # 2026-10-13). Probe points are suggested numerically but then VERIFIED
    # exactly to lie strictly inside the cell before they are used.
    crit = set()
    for p in [den_poly] + pos_polys + off_polys + stab_polys:
        if p.degree() < 1:
            continue
        for r in sp.real_roots(p):
            if r.is_real and (r - 1).is_positive:
                crit.add(r)
    crit = sorted(crit, key=lambda r: r.evalf(50))
    # distinct exact roots may compare equal only if they are the same number
    dedup = []
    for r in crit:
        if not dedup or not sp.simplify(r - dedup[-1]) == 0:
            dedup.append(r)
    crit = dedup

    def probe_between(lo, hi):
        """A rational strictly between exact numbers lo < hi (hi may be None)."""
        if hi is None:
            cand = Rational(sp.floor(lo) + 2)
            assert (cand - lo).is_positive
            return cand
        prec = 30
        while True:
            mid = (sp.N(lo, prec) + sp.N(hi, prec)) / 2
            cand = Rational(str(mid))
            if (cand - lo).is_positive and (hi - cand).is_positive:
                return cand
            prec *= 2
            if prec > 2000:
                raise ArithmeticError("could not separate consecutive critical values exactly")

    bounds = [Rational(1)] + list(crit) + [None]
    intervals = []
    for i in range(len(bounds) - 1):
        lo, hi = bounds[i], bounds[i + 1]
        probe = probe_between(lo, hi)
        dsign = den_poly.eval(probe)
        if dsign == 0:
            raise ArithmeticError("probe hit a root of the denominator -- impossible between critical values")
        ok = all(p.eval(probe) / dsign > 0 for p in pos_polys)
        if ok:
            ok = all(p.eval(probe) / dsign < 0 for p in off_polys)
        if ok:
            ok = all(p.eval(probe) > 0 for p in stab_polys)
        if ok:
            intervals.append((lo, hi))
        if verbose:
            print(f"   q in ({lo}, {hi}): {'STABLE FP' if ok else '-'}")

    return {"intervals": intervals, "critical_q": crit,
            "eps": e, "sigma": sorted(sigma)}
