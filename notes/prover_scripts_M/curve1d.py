"""
EXACT (no floating point). Lemma R test: is the full-support matrix K stable somewhere on the LEGAL BOUNDARY CURVE
    t (beta - 2) = 1,   i.e.  beta = 2 + 1/t,  t in (0, inf)?
By Lemma R (notes_lemmaM.md) this holds  <=>  -I + W_sigma (full support) is Hurwitz at SOME legal (eps, delta).
Hence 'never' here certifies: no legal (eps, delta) gives a stable full-support fixed point (positivity irrelevant).

Along the curve, t*K = L(t) := t I + (2t+1) N + t^2 J  (integer polynomial matrix; same stability as K for t > 0).
p(mu, t) = det(mu I + L(t)) = sum_k c_k(t) mu^(m-k),  deg_t c_k <= k+1.
Critical set: positive real roots of c_m(t) = det L(t) and of the Hurwitz minor Delta_{m-1}(t) (Lemma B, t as parameter).
Stability is constant on each open cell; decided exactly at a rational sample by ctln.hurwitz_stable.
"""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from fastregion import idet, icharpoly, interp, to_int, peval, hurwitz_det_poly, _cells, _hi_to_lo, _frac
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.domains import ZZ

def nonedge(A):
    m = len(A)
    return [[1 if (i != j and A[j][i] == 0) else 0 for j in range(m)] for i in range(m)]

def coeff_polys(A):
    """c_k(t) (low->high in t), k = 0..m, for det(mu I + L(t))."""
    m = len(A); N = nonedge(A)
    deg = m + 1
    xs = list(range(deg + 2))
    vals = []
    for x in xs:
        L = [[(x if i == j else 0) + (2 * x + 1) * N[i][j] + x * x for j in range(m)] for i in range(m)]
        cp = icharpoly([[-v for v in row] for row in L])   # det(mu I + L)
        vals.append(cp)
    P = []
    for k in range(m + 1):
        ys = [v[k] for v in vals]
        p = to_int(interp(xs[:deg + 1], ys[:deg + 1]))
        assert peval(p, xs[deg + 1]) == ys[deg + 1], "interpolation check failed"
        P.append(p)
    return P

def coeffs_at(P, t):
    return [peval(p, t) for p in P]

def curve_test(A):
    m = len(A)
    P = coeff_polys(A)
    assert P[0] == [1]
    zc = [k for k in range(1, m + 1) if not P[k]]
    if zc:
        # A coefficient of p(mu, t) vanishes IDENTICALLY in t: by (b1) of notes/lemmas.md (Hurwitz => all coefficients
        # > 0) p is Hurwitz at NO t. This is an identity over the whole curve, not a sampled boundary value.
        return {"verdict": "never", "reason": f"coefficient(s) {zc} identically zero along the curve", "n_critical": 0,
                "deg_Delta": None, "cells": [], "critical": []}
    Dp = P[m]
    Delta = hurwitz_det_poly(P, m - 1) if m >= 2 else [1]
    if not Dp or not Delta:
        return {"verdict": "never", "reason": "identically zero D or Delta"}
    iso = dup_isolate_real_roots_list([_hi_to_lo(Dp), _hi_to_lo(Delta)], ZZ, inf=0, strict=True, basis=True)
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == 0)]
    cells = _cells(iso, F(0))
    out_cells = []
    stable_any = False; escal = False
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

if __name__ == '__main__':
    from certify_M_violation import from_d6
    for s in sys.argv[1:]:
        r = curve_test(from_d6(s))
        print(s, r['verdict'], 'n_crit', r['n_critical'], 'deg Delta', r['deg_Delta'],
              [(c['lo'][:10], (c['hi'] or 'inf')[:10], c['stable']) for c in r['cells']])


def coeff_never(A, P=None):
    """EXACT early exit. True iff there is NO t > 0 at which all coefficients c_1..c_m of det(mu I + L(t)) are > 0.
    Then by (b1) (Hurwitz => all coefficients > 0) the curve is never stable. Cells between the positive real roots of
    all c_k have constant coefficient signs; roots themselves have a zero coefficient."""
    m = len(A)
    if P is None:
        P = coeff_polys(A)
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
    """curve_test with the exact coefficient-positivity early exit (same verdicts; about 4x faster on n=7 data)."""
    if coeff_never(A):
        return {"verdict": "never", "reason": "no t > 0 with all coefficients of det(mu I + L(t)) positive (b1)",
                "n_critical": None, "deg_Delta": None, "cells": [], "critical": []}
    return curve_test(A)
