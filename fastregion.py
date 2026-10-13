"""
fastregion.py -- exact parameter regions for ONE (graph, support), fast enough
to run over all 1,540,944 digraphs on 6 nodes.

Same mathematics as symbolic.q_region, different engine:

  * NO symbolic matrix algebra. Every boundary polynomial is obtained by
    evaluating integer matrices at integer parameter points (fraction-free
    determinants, integer Faddeev-LeVerrier) and interpolating exactly.
    Every interpolation is checked at one extra point.
  * Real roots are isolated exactly by sympy's VAS continued-fraction isolator
    (sympy.polys.rootisolation.dup_isolate_real_roots_list), which also returns
    the irreducible integer polynomial of each root. An endpoint is therefore
    serialised EXACTLY as (irreducible polynomial, isolating interval).
  * Each open cell between consecutive critical values is decided at a rational
    sample point with exact arithmetic; stability through ctln.hurwitz_stable
    (the validated three-valued routine). A vanishing minor at a sample point
    means re-sample; if it persists the cell is ESCALATE, never a verdict.

THE STRUCTURE THIS EXPLOITS (verified exactly, see ANOMALIES.md 2026-10-13)
--------------------------------------------------------------------------
    I - W_sigma = eps * K,      K = C + t J,      C = I + beta N_sigma
    t = (1-eps)/eps,   beta = 1 + delta/eps,   q = (beta-1) t/(t+1) = (beta-1)(1-eps)
    legal  q > 1  <=>  t (beta - 2) > 1  <=>  beta > 2 + 1/t
where N_sigma is the receiver-first NONEDGE matrix of the induced subgraph:
    N[i][j] = 1  iff  i != j and sigma_j -/-> sigma_i    (i.e. A[sigma_j][sigma_i] == 0).
(Same transpose as W: entry (i,j) is about the influence of j on i.)

Because J has rank one:
    adj(K) 1   = adj(C) 1                       =: a(beta)      (t-free)
    det K      = det C + t * sum(a)             =: d + t*s      (linear in t)
    det(mu I + K) = det(mu I + C) + t * 1^T adj(mu I + C) 1
               =: phi(mu) + t * omega(mu)                       (linear in t)
    off-support drive at k (times det K / theta) = d - beta * sum_j N[k][j] a_j =: o_k(beta)
So for a fixed t everything is a polynomial in beta alone, and the eps -> 0
limit (t -> infinity) is the one-parameter problem
    permitted:  all a_i(beta) nonzero with a common sign
    off:        o_k(beta) * s(beta) < 0
    stable:     omega(mu; beta) (degree m-1, leading coeff m) is Hurwitz in mu.

Conditions for "sigma supports a stable fixed point" at (beta, t):
    D := d + t s != 0;   a_i / D > 0  (all i in sigma);   o_k / D < 0  (all k outside);
    p(mu) := phi(mu) + t omega(mu) = det(mu I + K) is Hurwitz   (<=> -K Hurwitz <=> M = -I + W_sigma Hurwitz).

Boundary set in beta (where any verdict can change): roots of a_i, D, o_k, and
of the Hurwitz determinant Delta_{deg-1} of p (Orlando: a pair lambda_i + lambda_j = 0)
-- the constant term of p is D itself, so real-axis crossings are already in.
"""
from fractions import Fraction as F
from sympy.polys.rootisolation import dup_isolate_real_roots_list, dup_refine_real_root
from sympy.polys.domains import ZZ

import ctln

INF = None   # marker for t = infinity (eps -> 0 limit)

# ------------------------------------------------------------- integer algebra

def idet(M):
    """Bareiss fraction-free determinant of an integer matrix."""
    n = len(M)
    if n == 0:
        return 1
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        pk = M[k][k]
        Mk = M[k]
        for i in range(k + 1, n):
            Mi = M[i]
            mik = Mi[k]
            for j in range(k + 1, n):
                Mi[j] = (Mi[j] * pk - mik * Mk[j]) // prev
        prev = pk
    return sign * M[n - 1][n - 1]


def imatmul(X, Y):
    n = len(X)
    return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def icharpoly(M):
    """[1, c1, ..., cn] with det(lam I - M) = lam^n + c1 lam^(n-1) + ... + cn.
    Faddeev-LeVerrier; every division is exact over the integers."""
    n = len(M)
    coeffs = [1]
    Mk = None
    for k in range(1, n + 1):
        if k == 1:
            Mk = [row[:] for row in M]
        else:
            for i in range(n):
                Mk[i][i] += coeffs[-1]
            Mk = imatmul(M, Mk)
        tr = sum(Mk[i][i] for i in range(n))
        assert tr % k == 0, "Faddeev-LeVerrier division not exact -- bug"
        coeffs.append(-tr // k)
    return coeffs


def cramer_ones(M):
    """(adj(M) 1)_i = det(M with column i replaced by the all-ones vector)."""
    n = len(M)
    out = []
    for i in range(n):
        Mi = [row[:] for row in M]
        for r in range(n):
            Mi[r][i] = 1
        out.append(idet(Mi))
    return out

# ------------------------------------------------------------- polynomials (low -> high)

def interp(xs, ys):
    """Exact Newton interpolation; returns monomial coefficients low->high (Fractions)."""
    n = len(xs)
    c = [F(y) for y in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j])
    poly = [F(0)] * n
    basis = [F(1)]
    for i in range(n):
        for k, b in enumerate(basis):
            poly[k] += c[i] * b
        nb = [F(0)] * (len(basis) + 1)
        for k, b in enumerate(basis):
            nb[k] -= xs[i] * b
            nb[k + 1] += b
        basis = nb
    return poly


def to_int(poly):
    out = []
    for v in poly:
        v = F(v)
        if v.denominator != 1:
            raise ArithmeticError("non-integer coefficient where an integer was required")
        out.append(int(v))
    while out and out[-1] == 0:
        out.pop()
    return out


def peval(p, x):
    r = F(0)
    for c in reversed(p):
        r = r * x + c
    return r


def padd(p, q):
    n = max(len(p), len(q))
    out = [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]
    while out and out[-1] == 0:
        out.pop()
    return out


def pscale(p, c):
    out = [c * v for v in p]
    while out and out[-1] == 0:
        out.pop()
    return out


def pshift(p):
    """multiply by the variable"""
    return [0] + list(p) if p else []


def pdeg(p):
    return len(p) - 1


def pcompose_affine(p, c0, c1):
    """p(c0 + c1*x) as a polynomial in x, exact Fractions, low->high."""
    out = [F(0)]
    for c in reversed(p):
        # out = out*(c0 + c1 x) + c
        new = [F(0)] * (len(out) + 1)
        for k, v in enumerate(out):
            new[k] += v * c0
            new[k + 1] += v * c1
        new[0] += c
        out = new
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def primitive_int(p):
    """Scale a rational polynomial to a primitive integer polynomial (same roots)."""
    from math import gcd, lcm
    den = 1
    for v in p:
        den = lcm(den, F(v).denominator)
    ints = [int(F(v) * den) for v in p]
    g = 0
    for v in ints:
        g = gcd(g, abs(v))
    if g > 1:
        ints = [v // g for v in ints]
    while ints and ints[-1] == 0:
        ints.pop()
    return ints

# ------------------------------------------------------------- beta data

def nonedge_matrix(A, sigma):
    m = len(sigma)
    return [[0 if i == j else 1 - A[sigma[j]][sigma[i]] for j in range(m)] for i in range(m)]


def beta_data(A, sigma):
    """All t-free polynomials in beta for (A, sigma). Integer coefficient lists low->high."""
    sigma = sorted(sigma)
    m, n = len(sigma), len(A)
    Ns = nonedge_matrix(A, sigma)
    pts = list(range(m + 2))                      # m+1 for interpolation, 1 check
    dv, av, phiv, phiJv = [], [], [], []
    for b in pts:
        C = [[(1 if i == j else 0) + b * Ns[i][j] for j in range(m)] for i in range(m)]
        dv.append(idet(C))
        av.append(cramer_ones(C))
        phiv.append(icharpoly([[-C[i][j] for j in range(m)] for i in range(m)]))          # det(mu I + C)
        phiJv.append(icharpoly([[-C[i][j] - 1 for j in range(m)] for i in range(m)]))    # det(mu I + C + J)
    X = pts[:m + 1]
    chk = pts[m + 1]

    def fit(vals):
        p = to_int(interp(X, vals[:m + 1]))
        if peval(p, chk) != vals[m + 1]:
            raise ArithmeticError("interpolation check failed -- degree bound violated")
        return p

    d = fit(dv)
    a = [fit([av[j][i] for j in range(m + 2)]) for i in range(m)]
    phi = [fit([phiv[j][k] for j in range(m + 2)]) for k in range(m + 1)]
    phiJ = [fit([phiJv[j][k] for j in range(m + 2)]) for k in range(m + 1)]
    omega = [padd(phiJ[k], pscale(phi[k], -1)) for k in range(m + 1)]
    assert omega[0] == [] and (m == 0 or omega[1] == [m]), "omega leading structure wrong"
    s = a[0]
    for i in range(1, m):
        s = padd(s, a[i])
    o = {}
    for k in range(n):
        if k in sigma:
            continue
        acc = []
        for j in range(m):
            if A[sigma[j]][k] == 0:            # sigma_j -/-> k  : N[k][sigma_j] = 1
                acc = padd(acc, a[j])
        o[k] = padd(d, pscale(pshift(acc), -1))   # d - beta * sum_j N_kj a_j
    return {"sigma": sigma, "m": m, "n": n, "N": Ns, "d": d, "a": a, "s": s,
            "phi": phi, "omega": omega, "o": o}

# ------------------------------------------------------------- per-t polynomials

def stability_coeffs(bd, t):
    """Coefficient polynomials (in beta) of the stability polynomial, highest degree
    first, scaled to integers. t = INF gives omega (degree m-1, leading coeff m)."""
    m = bd["m"]
    if t is INF:
        return [list(bd["omega"][k]) for k in range(1, m + 1)]
    t = F(t)
    num, den = t.numerator, t.denominator
    return [padd(pscale(bd["phi"][k], den), pscale(bd["omega"][k], num)) for k in range(m + 1)]


def hurwitz_det_poly(P, k):
    """k-th leading principal Hurwitz minor of the polynomial with coefficient
    polynomials P (highest degree first), as an integer polynomial in beta,
    via exact evaluation + interpolation with a check point."""
    deg = len(P) - 1

    def ent(i, j):                       # 1-indexed Hurwitz matrix entry a_{2i-j}
        idx = 2 * i - j
        return P[idx] if 0 <= idx <= deg else []

    bound = 0
    for i in range(1, k + 1):
        bound += max((pdeg(ent(i, j)) for j in range(1, k + 1) if ent(i, j)), default=0)
    xs = list(range(bound + 2))
    vals = []
    for x in xs:
        H = [[int(peval(ent(i, j), x)) if ent(i, j) else 0 for j in range(1, k + 1)] for i in range(1, k + 1)]
        vals.append(idet(H))
    p = to_int(interp(xs[:bound + 1], vals[:bound + 1]))
    if peval(p, xs[bound + 1]) != vals[bound + 1]:
        raise ArithmeticError("Hurwitz interpolation check failed")
    return p


def D_poly(bd, t):
    if t is INF:
        return list(bd["s"])
    t = F(t)
    return padd(pscale(bd["d"], t.denominator), pscale(bd["s"], t.numerator))

# ------------------------------------------------------------- pointwise decision

def decide_point(bd, t, b0):
    """Exact verdict at beta = b0 (Fraction). Returns dict with positive/off_ok/stable."""
    m = bd["m"]
    if t is INF:
        D0 = peval(bd["s"], b0)
    else:
        D0 = peval(bd["d"], b0) + F(t) * peval(bd["s"], b0)
    if D0 == 0:
        return {"positive": None, "off_ok": None, "stable": None, "detail": "singular (D = 0)"}
    positive = all(peval(ai, b0) / D0 > 0 for ai in bd["a"])
    off_ok = all(peval(ok, b0) / D0 < 0 for ok in bd["o"].values())
    coeffs = [peval(p, b0) for p in stability_coeffs(bd, t)]
    stable, detail = ctln.hurwitz_stable(coeffs)
    return {"positive": positive, "off_ok": off_ok, "stable": stable, "detail": detail}

# ------------------------------------------------------------- region

def _hi_to_lo(p):
    return [int(c) for c in reversed(p)]


def region(A, sigma, t, bd=None, crosscheck=False, positivity_only=False, replay=False):
    """Exact set of beta (hence q) in the legal range for which sigma supports a
    stable fixed point, at fixed t = (1-eps)/eps (Fraction) or t = INF.

    Returns a certificate dict. Key fields:
      verdict        'never' | 'stable_somewhere' | 'escalate'
      stable_beta    list of cells [lo, hi] (endpoint descriptors) where all hold
      cells          every cell with its sample and verdict
      critical       every critical value: (interval, irreducible poly, which polys)
    """
    sigma = sorted(sigma)
    if bd is None:
        bd = beta_data(A, sigma)
    m = bd["m"]
    b_min = F(2) if t is INF else F(2) + 1 / F(t)
    cert = {"sigma": sigma, "t": None if t is INF else str(F(t)),
            "eps": "0 (limit)" if t is INF else str(1 / (F(t) + 1)),
            "beta_min": str(b_min)}

    names, polys = [], []
    for i, ai in enumerate(bd["a"]):
        names.append(f"a{i}"); polys.append(ai)
    for k, ok in bd["o"].items():
        names.append(f"o{k}"); polys.append(ok)
    zero = [nm for nm, p in zip(names, polys) if not p]
    if zero:
        cert.update(verdict="never", reason=f"identically zero: {zero}", anomaly=True)
        return cert

    # Stage 1: positivity (t-free) and off-support cells only.
    iso1 = dup_isolate_real_roots_list([_hi_to_lo(p) for p in polys], ZZ, inf=b_min, strict=True, basis=True)
    iso1 = [iv for iv in iso1 if not (iv[0][0] == iv[0][1] == b_min)]
    cells1 = _cells(iso1, b_min)
    # Lemma A.2 (notes/lemmas.md): a stable K has det K > 0, so x_i = a_i/D > 0 forces
    # a_i > 0 and the off-support drive o_k/D < 0 forces o_k < 0.
    sign_ok = False
    for _li, _hi, lo, hi, sample in cells1:
        if all(peval(ai, sample) > 0 for ai in bd["a"]) and all(peval(ok, sample) < 0 for ok in bd["o"].values()):
            sign_ok = True
            break
    if replay:
        _replay_isolation(iso1, [_hi_to_lo(p) for p in polys], b_min)
    if not sign_ok:
        cert.update(verdict="never", reason="no legal beta with all a_i > 0 and all o_k < 0",
                    n_critical_stage1=len(iso1))
        return cert
    if positivity_only:
        cert.update(verdict="positivity_possible", n_critical_stage1=len(iso1))
        return cert

    # Stage 2: full boundary set.
    P = stability_coeffs(bd, t)
    Dp = D_poly(bd, t)
    deg = len(P) - 1
    Delta = hurwitz_det_poly(P, deg - 1) if deg >= 2 else [1]
    names2 = names + ["D", "Delta"]
    polys2 = polys + [Dp, Delta]
    zero = [nm for nm, p in zip(names2, polys2) if not p]
    if zero:
        cert.update(verdict="never", reason=f"identically zero: {zero}", anomaly=True)
        return cert
    iso = dup_isolate_real_roots_list([_hi_to_lo(p) for p in polys2], ZZ, inf=b_min, strict=True, basis=True)
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == b_min)]
    if replay:
        _replay_isolation(iso, [_hi_to_lo(p) for p in polys2], b_min)
    cells = _cells(iso, b_min)

    crit = []
    for (a_, b_), idx, f in iso:
        crit.append({"interval": [str(_frac(a_)), str(_frac(b_))], "minpoly_hi_to_lo": [int(c) for c in f],
                     "vanishing": sorted(names2[i] for i in idx)})

    def endpoint(ci):
        if ci == -1:
            return {"kind": "beta_min", "beta": str(b_min), "q": "1"}
        if ci is None:
            return {"kind": "inf"}
        c = crit[ci]
        lo, hi = F(c["interval"][0]), F(c["interval"][1])
        e = {"kind": "root", "index": ci, "beta_interval": c["interval"], "beta_minpoly_hi_to_lo": c["minpoly_hi_to_lo"],
             "vanishing": c["vanishing"]}
        # exact image in q:  q = (beta-1) t/(t+1)  (t finite)  or  q = beta - 1 (limit)
        if t is INF:
            c0, c1 = F(1), F(1)            # beta = 1 + q
        else:
            c0, c1 = F(1), (F(t) + 1) / F(t)   # beta = 1 + q (t+1)/t
        qpoly = primitive_int(pcompose_affine([F(v) for v in reversed(c["minpoly_hi_to_lo"])], c0, c1))
        e["q_interval"] = [_q_of_beta(lo, t), _q_of_beta(hi, t)]
        e["q_minpoly_hi_to_lo"] = list(reversed(qpoly))
        if lo == hi:
            e["beta"] = str(lo); e["q"] = _q_of_beta(lo, t)
        return e

    out_cells, stable_cells, escal = [], [], False
    for lo_i, hi_i, lo, hi, sample in cells:
        v = None
        tries = [sample]
        if hi is not None:
            tries += [lo + (hi - lo) / 3, lo + 2 * (hi - lo) / 3]
        else:
            tries += [lo + 2, lo + 3]
        for s0 in tries:
            v = decide_point(bd, t, s0)
            if v["stable"] is not None and v["positive"] is not None:
                sample = s0
                break
        rec = {"lo": lo_i, "hi": hi_i, "sample_beta": str(sample), "sample_q": _q_of_beta(sample, t), **v}
        if crosscheck and t is not INF:
            eps = 1 / (F(t) + 1)
            delta = (sample - 1) * eps
            r = ctln.analyze_support(A, sigma, eps, delta)
            agree = (r["positive"] == v["positive"] and r["off_support_ok"] == v["off_ok"] and r["stable"] == v["stable"])
            rec["crosscheck_analyze_support"] = agree
            if not agree:
                raise AssertionError(f"fastregion disagrees with ctln.analyze_support at beta={sample}: {v} vs "
                                     f"{(r['positive'], r['off_support_ok'], r['stable'])}")
        out_cells.append(rec)
        if v["stable"] is None or v["positive"] is None:
            escal = True
        elif v["positive"] and v["off_ok"] and v["stable"]:
            stable_cells.append(rec)

    cert.update(critical=crit, cells=out_cells, Delta_degree=pdeg(Delta), escalated_cells=escal)
    if stable_cells:
        # a decided stable cell is a hit regardless of other cells escalating
        cert["verdict"] = "stable_somewhere"
        cert["stable_intervals"] = [{"lo": endpoint(c["lo"]), "hi": endpoint(c["hi"])} for c in stable_cells]
    elif escal:
        cert["verdict"] = "escalate"
    else:
        cert["verdict"] = "never"
    return cert


def _frac(x):
    return F(int(x.numerator), int(x.denominator))


def _q_of_beta(b, t):
    b = F(b)
    return str(b - 1) if t is INF else str((b - 1) * F(t) / (F(t) + 1))


def _cells(iso, b_min):
    """Open cells (lo_idx, hi_idx, lo, hi, sample) between consecutive isolating
    intervals; idx -1 is beta_min, None is +infinity. Refines the first interval
    until it lies strictly above b_min (its root does, since exact-b_min roots
    were dropped)."""
    bounds = []
    for (a_, b_), idx, f in iso:
        bounds.append([_frac(a_), _frac(b_), f])
    if bounds and bounds[0][0] <= b_min:
        a_, b_, f = bounds[0]
        while a_ <= b_min:
            a_, b_ = dup_refine_real_root(f, a_, b_, ZZ, steps=1)
            a_, b_ = _frac(a_), _frac(b_)
        bounds[0][0], bounds[0][1] = a_, b_
    cells = []
    prev, prev_i = b_min, -1
    for ci, (a_, b_, f) in enumerate(bounds):
        assert prev < a_, "isolating intervals not strictly separated"
        cells.append((prev_i, ci, prev, a_, (prev + a_) / 2))
        prev, prev_i = b_, ci
    cells.append((prev_i, None, prev, None, prev + 1))
    return cells


def stable_q_intervals(cert):
    """Human-readable exact q-intervals from a certificate."""
    out = []
    for iv in cert.get("stable_intervals", []):
        def show(e):
            if e["kind"] == "beta_min": return "1"
            if e["kind"] == "inf": return "inf"
            if "q" in e: return e["q"]
            return f"root of {e['q_minpoly_hi_to_lo']} in ({e['q_interval'][0]}, {e['q_interval'][1]})"
        out.append((show(iv["lo"]), show(iv["hi"])))
    return out


# ------------------------------------------------------------- independent replay of root isolation

def _pdivmod(a, b):
    """polynomial division over Fractions, low->high lists; returns remainder."""
    a = [F(x) for x in a]
    b = [F(x) for x in b]
    while b and b[-1] == 0:
        b.pop()
    db = len(b) - 1
    while len(a) - 1 >= db and any(a):
        if a[-1] == 0:
            a.pop(); continue
        coef = a[-1] / b[-1]
        shift = len(a) - 1 - db
        for i, bv in enumerate(b):
            a[shift + i] -= coef * bv
        a.pop()
    while a and a[-1] == 0:
        a.pop()
    return a


def sturm_sequence(p):
    p = [F(x) for x in p]
    while p and p[-1] == 0:
        p.pop()
    dp = [F(i) * p[i] for i in range(1, len(p))]
    seq = [p, dp]
    while seq[-1]:
        r = _pdivmod(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-x for x in r])
    return seq


def _sign_changes(vals):
    s = [v for v in vals if v != 0]
    return sum(1 for i in range(len(s) - 1) if (s[i] > 0) != (s[i + 1] > 0))


def sturm_count(p, a, b=None):
    """Number of DISTINCT real roots of p in (a, b]  (b None = +infinity). Requires p(a) != 0."""
    seq = sturm_sequence(p)
    va = _sign_changes([peval(s, a) for s in seq])
    if b is None:
        vb = _sign_changes([s[-1] for s in seq])          # signs at +infinity: leading coefficients
    else:
        vb = _sign_changes([peval(s, b) for s in seq])
    return va - vb


def _replay_isolation(iso, polys_hi_to_lo, b_min):
    """Independent check of sympy's isolation output:
       * every non-degenerate isolating interval shows a sign change of its irreducible
         factor (so it contains a root), every degenerate one is an exact rational root;
       * for every input polynomial the number of intervals attributed to it equals the
         Sturm count of its distinct real roots in (b_min, +inf).
    Raises AssertionError on any discrepancy."""
    per_poly = {}
    for (a_, c_), idx, f in iso:
        a_, c_ = _frac(a_), _frac(c_)
        flo = [F(int(v)) for v in reversed(f)]
        if a_ == c_:
            assert peval(flo, a_) == 0, "degenerate isolating interval is not a root"
        else:
            assert peval(flo, a_) * peval(flo, c_) < 0, "isolating interval shows no sign change"
        for i in idx:
            per_poly[i] = per_poly.get(i, 0) + 1
    for i, p in enumerate(polys_hi_to_lo):
        plo = [F(int(v)) for v in reversed(p)]
        if len(plo) <= 1:
            continue
        if peval(plo, b_min) == 0:
            # a root exactly at b_min was filtered out; count roots strictly above by a tiny shift
            # (exact): use the deflated polynomial p/(beta - b_min)
            q = plo
            while peval(q, b_min) == 0:
                q = _pdivmod_exact(q, b_min)
            plo = q
        n_sturm = sturm_count(plo, b_min)
        assert n_sturm == per_poly.get(i, 0), f"Sturm count {n_sturm} != {per_poly.get(i, 0)} isolating intervals for polynomial {i}"


def _pdivmod_exact(p, r):
    """exact division of p (low->high) by (beta - r); p(r) must be 0."""
    hi = list(reversed(p))
    out = []
    acc = F(0)
    for c in hi:
        acc = acc * r + c
        out.append(acc)
    assert out[-1] == 0
    return list(reversed(out[:-1]))
