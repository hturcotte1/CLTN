"""
curve_test.py -- EXACT Lemma-R test: is the full-support matrix K = I + beta N + tJ positive stable
somewhere on the LEGAL BOUNDARY CURVE t(beta - 2) = 1 (beta = 2 + 1/t, t > 0)?

Lemma R (notes/lemma_M_refuted.md section 3, PROVED; re-derived in ANOMALIES.md 2026-10-14):
  K(lambda beta, lambda t) = (1 - lambda) I + lambda K(beta, t)   (Remark F), so stability at a legal
  (beta, t) propagates down the ray to lambda in (0, 1]; f(lambda) = lambda t (lambda beta - 2) - 1 has
  f(2/beta) = -1 < 0 < f(1), so the ray meets the curve at a stable point with beta' > 2.  Conversely a
  stable curve point has stable legal neighbours (openness).  Hence
      the full support is stable at SOME legal (eps, delta)  <=>  K(2 + 1/t, t) is stable for SOME t > 0.
  'never' therefore certifies: no legal (eps, delta) gives a stable full-support fixed point, positivity
  irrelevant; so the graph is never a counterexample (stable positive non-clique full support).

Along the curve t K(2 + 1/t, t) = L(t) := t I + (2t + 1) N + t^2 J, an integer polynomial matrix with the
same stability as K for t > 0.  p(mu, t) = det(mu I + L(t)) = sum_k c_k(t) mu^(m-k), deg_t c_k <= k + 1
(c_k = sum of k x k principal minors; matrix determinant lemma on each block).  The coefficients are
found by exact integer evaluation at m + 3 points, interpolation and one check point.  Critical set:
positive real roots of c_m(t) = det L(t) and of the Hurwitz minor Delta_{m-1}(t) (Lemma B with t as the
parameter: a root of p crosses the imaginary axis only where c_m = 0 or, by Orlando, Delta_{m-1} = 0).
Stability is constant on each open cell between consecutive roots and is decided at an exact rational
sample by ctln.hurwitz_stable (three-valued; a cell undecided at three samples is 'escalate', never
'never').  A coefficient c_k identically zero in t is 'never' (Hurwitz => all coefficients > 0, (b1)).
Exact fast path (curve_test_fast): if no t > 0 has all c_1..c_m > 0 the curve is never stable (b1);
decided from the positive roots of the low-degree c_k.  Same verdicts, ~4x faster.

Validation (see n7/validate_curve_test.py): n = 5 all 9,608 digraphs; n = 6 all 28,350 strict positivity
survivors (every graph stable somewhere in the exact 2-D map must be stable_somewhere here).
No floating point anywhere.
"""
import sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ctln
from fastregion import icharpoly, interp, to_int, peval, hurwitz_det_poly, _cells, _hi_to_lo, _frac
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.domains import ZZ


def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


def nonedge(A):
    """N[i][j] = 1 iff j -/-> i (i != j): receiver-first nonedge matrix, as fastregion.beta_data."""
    m = len(A)
    return [[1 if (i != j and A[j][i] == 0) else 0 for j in range(m)] for i in range(m)]


def coeff_polys(A):
    """c_k(t), low->high integer coefficient lists, k = 0..m, for det(mu I + L(t))."""
    m = len(A); N = nonedge(A)
    deg = m + 1
    xs = list(range(deg + 2))
    vals = []
    for x in xs:
        L = [[(x if i == j else 0) + (2 * x + 1) * N[i][j] + x * x for j in range(m)] for i in range(m)]
        vals.append(icharpoly([[-v for v in row] for row in L]))      # det(mu I + L)
    P = []
    for k in range(m + 1):
        ys = [v[k] for v in vals]
        p = to_int(interp(xs[:deg + 1], ys[:deg + 1]))
        if peval(p, xs[deg + 1]) != ys[deg + 1]:
            raise ArithmeticError("interpolation check failed -- degree bound violated")
        P.append(p)
    return P


def coeffs_at(P, t):
    return [peval(p, t) for p in P]


def curve_test(A, P=None):
    m = len(A)
    P = P or coeff_polys(A)
    assert P[0] == [1], "p is not monic"
    zc = [k for k in range(1, m + 1) if not P[k]]
    if zc:
        return {"verdict": "never", "reason": f"coefficient(s) {zc} identically zero along the curve (b1)",
                "n_critical": 0, "deg_Delta": None, "cells": [], "critical": []}
    Dp = P[m]
    Delta = hurwitz_det_poly(P, m - 1) if m >= 2 else [1]
    if not Dp or not Delta:
        return {"verdict": "never", "reason": "identically zero D or Delta (b1)", "n_critical": 0, "deg_Delta": None, "cells": [], "critical": []}
    iso = dup_isolate_real_roots_list([_hi_to_lo(Dp), _hi_to_lo(Delta)], ZZ, inf=0, strict=True, basis=True)
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == 0)]
    cells = _cells(iso, F(0))
    out_cells = []; stable_any = False; escal = False
    for lo_i, hi_i, lo, hi, sample in cells:
        tries = [sample] + ([lo + (hi - lo) / 3, lo + 2 * (hi - lo) / 3] if hi is not None else [lo + 2, lo + 7])
        v = None
        for s0 in tries:
            v, det = ctln.hurwitz_stable(coeffs_at(P, s0))
            if v is not None:
                sample = s0; break
        out_cells.append({"lo": str(lo), "hi": None if hi is None else str(hi), "sample_t": str(sample), "stable": v})
        if v is None: escal = True
        if v: stable_any = True
    verdict = "stable_somewhere" if stable_any else ("escalate" if escal else "never")
    return {"verdict": verdict, "n_critical": len(iso), "deg_Delta": len(Delta) - 1, "cells": out_cells,
            "critical": [{"interval": [str(_frac(a)), str(_frac(b))], "minpoly_hi_to_lo": [int(c) for c in f]} for (a, b), idx, f in iso]}


def coeff_never(A, P=None):
    """True iff NO t > 0 has all coefficients c_1..c_m of det(mu I + L(t)) positive (then never stable, (b1))."""
    m = len(A)
    P = P or coeff_polys(A)
    if any(not P[k] for k in range(1, m + 1)):
        return True
    nonconst = [_hi_to_lo(P[k]) for k in range(1, m + 1) if len(P[k]) > 1]
    iso = dup_isolate_real_roots_list(nonconst, ZZ, inf=0, strict=True, basis=True) if nonconst else []
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == 0)]
    for _, _, lo, hi, s in _cells(iso, F(0)):
        if all(peval(P[k], s) > 0 for k in range(1, m + 1)):
            return False
    return True


def curve_test_fast(A):
    P = coeff_polys(A)
    if coeff_never(A, P):
        return {"verdict": "never", "reason": "no t > 0 with all coefficients of det(mu I + L(t)) positive (b1)",
                "n_critical": None, "deg_Delta": None, "cells": [], "critical": []}
    return curve_test(A, P)


if __name__ == "__main__":
    for s in sys.argv[1:]:
        r = curve_test_fast(from_d6(s))
        print(s, r["verdict"], r.get("reason", ""), [(c["lo"], c["hi"], c["stable"]) for c in r["cells"] if c["stable"]])
