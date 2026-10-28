"""OWN replication of the eps -> 0 limit stage over ALL beta > 2 (what 'never_positivity' / 'never_full' certify).

Derivation (own): with eps = 1/(t+1), delta = (beta-1) eps,  (t+1)(I - W) = K := C + t J,  C = I + beta N,
N[i][j] = 1 iff i != j and A[j][i] == 0.  Fixed point x = (t+1) K^{-1} 1; stable iff -I + W = -K/(t+1) is Hurwitz
iff det(mu I + K) is Hurwitz in mu.  det(mu I + K) = det(mu I + C) + t * omega(mu; beta) with
omega = det(mu I + C + J) - det(mu I + C) (matrix determinant lemma), and K^{-1} 1 = adj(C) 1 / (det C + t * s),
s = 1^T adj(C) 1, a := adj(C) 1 (t-free).  A stable K has det K > 0 (product of eigenvalues with positive real
parts), so positivity forces a_i > 0 for every i, at EVERY t.  The limit stage (t = infinity) declares
   never_positivity : no beta > 2 with all a_i(beta) > 0           (hence no stable positive point at ANY t, any legal beta)
   never_full       : positivity cells exist, but on no open cell of beta > 2 (cells cut by the roots of a_i, s and
                      the Orlando minor Delta_6 of omega) is omega Hurwitz AND all a_i > 0.
Everything exact: integer evaluation + own interpolation; real roots isolated by sympy over ZZ; cells decided at
rational sample points with own Hurwitz."""
import json, sys, time
sys.path.insert(0, '.')
from fractions import Fraction as F
from sympy import Poly, symbols, ZZ, QQ, Rational
from engine import parse_d6, charpoly_int, det_int, hurwitz
beta_sym = symbols('beta')

def interp(xs, ys):
    """Newton interpolation, exact; returns coefficients low->high (Fractions)."""
    n = len(xs); c = [F(y) for y in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j])
    poly = [F(0)] * n; basis = [F(1)]
    for i in range(n):
        for k, b in enumerate(basis): poly[k] += c[i] * b
        nb = [F(0)] * (len(basis) + 1)
        for k, b in enumerate(basis): nb[k] -= xs[i] * b; nb[k + 1] += b
        basis = nb
    while len(poly) > 1 and poly[-1] == 0: poly.pop()
    return poly

def fit_int(vals, deg):
    """vals at x = 0..deg+1 (deg+2 values); fit degree <= deg, check the extra point, return int coeffs low->high."""
    p = interp(list(range(deg + 1)), vals[:deg + 1])
    assert all(v.denominator == 1 for v in p)
    p = [int(v) for v in p]
    assert peval(p, deg + 1) == vals[deg + 1], "interpolation degree bound violated"
    return p

def peval(p, x):
    r = F(0)
    for c in reversed(p): r = r * x + c
    return r

def Cmat(A, b):
    n = len(A)
    return [[(1 if i == j else 0) + (b if (i != j and A[j][i] == 0) else 0) for j in range(n)] for i in range(n)]

def limit_data(A):
    n = len(A); pts = range(n + 2)
    a_v, s_v, om_v = [], [], []
    for b in pts:
        C = Cmat(A, b)
        a = []
        for i in range(n):
            Ci = [row[:] for row in C]
            for r in range(n): Ci[r][i] = 1
            a.append(det_int(Ci))                                   # (adj(C) 1)_i by Cramer
        a_v.append(a); s_v.append(sum(a))
        cpC = charpoly_int([[-v for v in row] for row in C])        # det(mu I + C), highest first
        cpCJ = charpoly_int([[-v - 1 for v in row] for row in C])   # det(mu I + C + J)
        om_v.append([cpCJ[k] - cpC[k] for k in range(n + 1)])       # omega coefficients, highest first (k=0 is 0)
    a = [fit_int([a_v[b][i] for b in pts], n - 1) for i in range(n)]
    s = fit_int(s_v, n - 1)
    omega = [fit_int([om_v[b][k] for b in pts], k) for k in range(n + 1)]   # coefficient of mu^(n-k): degree <= k in beta
    assert omega[0] == [] or omega[0] == [0] or all(c == 0 for c in omega[0])
    omega = omega[1:]                                               # degree n-1 in mu, leading coefficient n
    assert omega[0] == [n]
    assert omega[-1] == s or [c for c in omega[-1]] == s, "constant term of omega must be s"
    return a, s, omega

def hurwitz_minor_poly(omega, k):
    """k-th leading principal Hurwitz minor of omega(mu; beta) as an integer polynomial in beta (own interpolation)."""
    d = len(omega) - 1
    def ent(i, j):
        idx = 2 * i - j
        return omega[idx] if 0 <= idx <= d else []
    bound = sum(max((len(ent(i, j)) - 1 for j in range(1, k + 1) if ent(i, j)), default=0) for i in range(1, k + 1))
    vals = [det_int([[int(peval(ent(i, j), x)) if ent(i, j) else 0 for j in range(1, k + 1)] for i in range(1, k + 1)]) for x in range(bound + 2)]
    return fit_int(vals, bound)

def cells_above_2(polys):
    """Open cells of beta > 2 cut by the real roots (> 2) of the given integer polynomials (low->high):
    returns rational sample points, one inside every open cell, plus the list of root isolating intervals."""
    P = Poly(1, beta_sym, domain=ZZ)
    for p in polys:
        if len(p) > 1: P = P * Poly(list(reversed(p)), beta_sym, domain=ZZ)
    # remove roots exactly at 2
    two = Poly([1, -2], beta_sym, domain=ZZ)
    while P.eval(2) == 0: P = P.quo(two)
    P = P.sqf_part()                              # squarefree part: same roots, isolation and refinement well-defined
    if P.degree() == 0: return [F(3)], []
    iv = P.intervals(inf=2)
    iv = [(Rational(a), Rational(b)) for (a, b), _ in iv]
    iv = [(a, b) for a, b in iv if b > 2]
    if iv and iv[0][0] <= 2:                      # refine the first interval until it lies strictly above 2
        a, b = iv[0]
        while a <= 2:
            a, b = P.refine_root(a, b, steps=1)
        iv[0] = (Rational(a), Rational(b))
    # open every touching pair of isolating intervals (shared endpoint) by exact refinement, so that each open cell
    # between consecutive roots contains a non-degenerate rational gap
    for k in range(1, len(iv)):
        while iv[k - 1][1] >= iv[k][0]:
            assert iv[k - 1][1] == iv[k][0], "isolating intervals overlap"
            a0, b0 = iv[k - 1]; a1, b1 = iv[k]
            if b0 > a0: a0, b0 = (Rational(v) for v in P.refine_root(a0, b0, steps=1))
            if b1 > a1: a1, b1 = (Rational(v) for v in P.refine_root(a1, b1, steps=1))
            assert (b0 > a0) or (b1 > a1), "two rational roots coincide -- impossible for squarefree P"
            iv[k - 1] = (a0, b0); iv[k] = (a1, b1)
    for k in range(1, len(iv)): assert iv[k - 1][1] < iv[k][0], "isolating intervals still touch"
    bounds = [F(2)] + [F(int(x.p), int(x.q)) for a, b in iv for x in (a, b)]
    samples = []; cells = []
    for k in range(0, len(bounds) - 1, 2):
        lo, hi = bounds[k], bounds[k + 1]          # gap between consecutive isolating intervals (may be a single point)
        samples.append((lo + hi) / 2 if hi > lo else lo); cells.append((lo, hi))
    samples.append(bounds[-1] + 1); cells.append((bounds[-1], None))
    cells_above_2.last_cells = cells
    # every sample must lie in an open cell: no boundary polynomial vanishes there, and it is > 2
    Pfull = Poly(1, beta_sym, domain=ZZ)
    for p in polys:
        if len(p) > 1: Pfull = Pfull * Poly(list(reversed(p)), beta_sym, domain=ZZ)
    for x in samples:
        assert x > 2 and Pfull.eval(Rational(x.numerator, x.denominator)) != 0, "sample point is a root -- bug"
    return samples, [(str(a), str(b)) for a, b in iv]

def replay(A):
    n = len(A)
    a, s, omega = limit_data(A)
    if any(len(p) == 0 or all(c == 0 for c in p) for p in a):
        return {'verdict': 'never', 'stage': 'positivity', 'reason': 'some a_i identically zero'}
    samp1, iv1 = cells_above_2(a)
    pos_cells = [str(x) for x in samp1 if all(peval(ai, x) > 0 for ai in a)]
    if not pos_cells:
        return {'verdict': 'never', 'stage': 'positivity', 'n_roots_stage1': len(iv1), 'n_cells': len(samp1)}
    Delta = hurwitz_minor_poly(omega, n - 2)                 # Orlando minor Delta_{d-1}, d = n-1
    samp2, iv2 = cells_above_2(a + [s, Delta])
    cells = []; stable_pos = []; undecided = []
    cell_bounds = cells_above_2.last_cells
    for x, (lo, hi) in zip(samp2, cell_bounds):
        # the verdict is constant on the open cell; a sample where a coefficient or minor vanishes is re-sampled
        # at other points of the SAME cell (the point itself stays undecided; the cell is decided by a clean point)
        tries = [x] + ([lo + (hi - lo) / 3, lo + 2 * (hi - lo) / 3] if hi is not None and hi > lo else [lo + 2, lo + 3] if hi is None else [])
        for x in tries:
            positive = all(peval(ai, x) > 0 for ai in a)
            coeffs = [peval(p, x) for p in omega]
            den = 1
            for c in coeffs: den = den * c.denominator // __import__('math').gcd(den, c.denominator)
            st, det = hurwitz([int(c * den) for c in coeffs])
            if st is not None: break
        cells.append((str(x), positive, st))
        if st is None: undecided.append((str(x), det))
        if positive and st is True: stable_pos.append(str(x))
    v = 'stable_somewhere' if stable_pos else ('escalate' if undecided else 'never')
    return {'verdict': v, 'stage': 'full', 'n_roots_stage2': len(iv2), 'n_cells': len(samp2), 'stable_positive_cells': stable_pos,
            'undecided': undecided, 'positivity_cells_stage1': pos_cells, 'Delta_degree': len(Delta) - 1}

if __name__ == '__main__':
    sample, outp = sys.argv[1], sys.argv[2]
    t0 = time.time(); tally = {}
    with open(outp, 'w') as fh:
        for line in open(sample):
            r = json.loads(line); A = parse_d6(r['d6'])
            res = replay(A); res['d6'] = r['d6']; res['cls'] = r['cls']
            key = (r['cls'], res['verdict'], res['stage']); tally[key] = tally.get(key, 0) + 1
            if res['verdict'] != 'never': print("*** DISAGREEMENT with limit stage:", r['d6'], res, flush=True)
            fh.write(json.dumps(res) + '\n')
    print(sample, tally, f"{time.time()-t0:.0f}s")
