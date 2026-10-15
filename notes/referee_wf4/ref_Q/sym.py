"""Referee ref_Q: re-derive the limit polynomials of Q, Q2, Y symbolically in k, check them against explicit
graphs, verify the thresholds by exact root isolation in k, and the beta>2 positivity. Exact (sympy rationals/ints)."""
import sys, json
import sympy as sp
from fractions import Fraction as F
import fam, exact
z, k, x, b = sp.symbols('z k x beta')

def Dmat(A):
    n = len(A); return [[1 if (i != j and not A[i][j]) else 0 for j in range(n)] for i in range(n)]

def family_Phi(B, classes):
    """Phi(z,k) = z psi_D + p chi_D + 1^T adj(zI-D) Bvec  with p = sum sizes, Bvec = sum size_j b_j (sizes linear in k)."""
    m = len(B); D = sp.Matrix(Dmat(B)); Z = z * sp.eye(m) - D
    adj = Z.adjugate(); one = sp.ones(m, 1)
    psi = sp.expand((one.T * adj * one)[0]); chi = sp.expand(Z.det())
    p = sum(c * k for _, c in classes)
    Bv = sp.zeros(m, 1)
    for S, c in classes:
        Bv += c * k * sp.Matrix([1 - s for s in S])
    Phi = sp.expand(z * psi + p * chi + (one.T * adj * Bv)[0])
    return Phi, sp.expand(z * psi), sp.expand(Phi.coeff(k, 1)), sp.expand(chi)

def explicit_psi(A):
    c = exact.psi_int(Dmat(A)); return sp.Poly(list(reversed(c)), z).as_expr()

def hurwitz_minors_k(Phi, bstar):
    dz = sp.Poly(Phi, z).degree()
    h = sp.Poly(sp.expand((-1) ** dz * Phi.subs(z, -1 / bstar - x)), x)   # leading coeff = LC_z(Phi) > 0
    co = [sp.expand(c) for c in h.all_coeffs()]
    m = len(co) - 1
    def a(i): return co[i] if 0 <= i <= m else 0
    H = sp.Matrix(m, m, lambda i, j: a(2 * (i + 1) - (j + 1)))
    mins = [sp.expand(H[:r, :r].det(method='berkowitz')) for r in range(1, m + 1)]
    return co, mins

def threshold(Phi, bstar, k0, kmax_scan=None):
    co, mins = hurwitz_minors_k(Phi, bstar)
    rec = []; allroots = []
    ok = True
    for mn in [co[0]] + mins:
        P = sp.Poly(mn, k)
        lc = P.LC(); deg = P.degree()
        nge = P.count_roots(k0, None) if deg > 0 else 0   # roots in [k0, inf)
        rr = [r for r in sp.real_roots(P)] if deg > 0 else []
        allroots += rr
        rec.append(dict(deg=deg, LC_pos=bool(lc > 0), roots_ge_k0=int(nge)))
        ok = ok and lc > 0 and nge == 0
    # minimal integer threshold: smallest K such that every minor > 0 for all integers k >= K
    big = max([sp.ceiling(r) for r in allroots if r.is_real] + [1])
    # scan downward from max real root to find exact set of good integers near the top
    def good(kk): return co[0].subs(k, kk) > 0 and all(mn.subs(k, kk) > 0 for mn in mins)
    kmin = None
    K = int(big) + 1
    assert good(K)
    while K - 1 >= 1 and good(K - 1): K -= 1
    return ok, rec, K   # K = least integer such that all k in [K, max_root+1] are good (and beyond max root all good)

def positivity(B, Slist):
    m = len(B); g = 1 - 1 / b
    At = sp.Matrix(B).T
    Mm = (b * (At + g * sp.eye(m))).applyfunc(sp.expand)   # = b A^T + (b-1) I, polynomial
    det = sp.expand(Mm.det()); adj = Mm.adjugate()
    one = sp.ones(m, 1); v = (adj * one).applyfunc(sp.expand)    # w = b * adj(Mm) 1 / det(Mm)
    res = {}
    def nr(p):
        P = sp.Poly(p, b); return int(P.count_roots(2, None)) if P.degree() > 0 else 0
    sgn3 = sp.sign(det.subs(b, 3))
    res['det'] = (str(sp.factor(det)), nr(det), int(det.subs(b, 3)))
    ok = nr(det) == 0
    for i in range(m):
        p = v[i]; res['w%d' % i] = (str(sp.factor(p)), nr(p), int(p.subs(b, 3)))
        ok = ok and nr(p) == 0 and sp.sign(p.subs(b, 3)) == sgn3
    for S in Slist:
        # sigma_S - 1 = (b * sum_{i in S} v_i - det)/det
        num = sp.expand(b * sum(v[i] for i in range(m) if S[i]) - det)
        res['sigma-1 S=%s' % ''.join(map(str, S))] = (str(sp.factor(num)), nr(num), int(num.subs(b, 3)))
        ok = ok and nr(num) == 0 and sp.sign(num.subs(b, 3)) == sgn3
    return ok, res

FAMS = {
  'Q':  dict(B=fam.BQ,  classes=[(fam.SQ, 1)], build=fam.Q, ks=[1, 2, 3, 4, 5], thr=[(sp.Rational(5, 2), 25), (sp.Rational(27, 10), 79), (sp.Rational(14, 5), 286)]),
  'Q2': dict(B=fam.BQ2, classes=[(fam.SQ2a, 1), (fam.SQ2b, 3)], build=fam.Q2, ks=[1, 2, 3], thr=[(sp.Rational(5, 2), 1), (sp.Integer(3), 3), (sp.Rational(7, 2), 10), (sp.Rational(19, 5), 33)]),
  'Y':  dict(B=fam.BY,  classes=[(fam.SY1, 4), (fam.SY2, 1)], build=fam.Y, ks=[1, 2], thr=[(sp.Integer(3), 13), (sp.Rational(7, 2), 24), (sp.Integer(4), 46), (sp.Rational(9, 2), 117), (sp.Rational(24, 5), 363)]),
}
if __name__ == '__main__':
    which = sys.argv[1:] or list(FAMS)
    out = {}
    for name in which:
        f = FAMS[name]; B = f['B']; m = len(B)
        Phi, zpsi, Lk, chi = family_Phi(B, f['classes'])
        p_per_k = sum(c for _, c in f['classes'])
        r = {'Phi(z,k)': str(sp.factor(Phi)), 'z_psi_D': str(sp.factor(zpsi)), 'coeff_of_k': str(sp.factor(Lk)),
             'limit_poly_monic': str(sp.factor(Lk / p_per_k))}
        # explicit graphs: psi_{D'} == z^(p-1) Phi(z, k)
        ex = []
        for kk in f['ks']:
            A = f['build'](kk); p = p_per_k * kk
            lhs = explicit_psi(A); rhs = sp.expand(z ** (p - 1) * Phi.subs(k, kk))
            ex.append((kk, len(A), sp.expand(lhs - rhs) == 0))
        r['explicit_check'] = ex
        ok, pos = positivity(B, [S for S, _ in f['classes']])
        r['positivity_all_beta_gt_2'] = ok; r['positivity_polys'] = pos
        thr = {}
        for bstar, k0 in f['thr']:
            okt, rec, Kmin = threshold(Phi, bstar, k0)
            thr[str(bstar)] = dict(k0_claimed=k0, all_minors_LCpos_and_no_root_ge_k0=bool(okt), least_K_all_good_from_K_on=Kmin, minors=rec)
        r['thresholds'] = thr
        out[name] = r
        print(name, json.dumps({kk: v for kk, v in r.items() if kk != 'positivity_polys'}, indent=0, default=str), flush=True)
    json.dump(out, open('sym_%s.json' % '_'.join(which), 'w'), indent=1, default=str)
