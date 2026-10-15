"""Exact first-order coefficient R* (Lemma AS, notes/lemma_M_refuted.md section 4) for the 124 n=7 finite-eps
candidates.  Pure exact arithmetic (sympy integers/rationals, Fractions); no float enters any decision.

Conventions (CLAUDE.md section 2; lemma_M_refuted.md section 1):
  A[i][j] = 1 iff i -> j.  chi(z) = det(zI + A),  psi(z) = det(zI + A + J) - det(zI + A) = 1^T adj(zI + A) 1,
  det(zI - (cJ - A)) = chi - c psi.  h = gcd(chi, psi), chi = h chi1, psi = h psi1.
  nu_inf = min Re(roots of psi) = -alpha,  beta_c = 1/(1 - nu_inf).
  For a simple root z0 of psi1 with Re z0 = nu_inf:  moving branch z0 + w/c + O(c^-2),  w = chi1(z0)/psi1'(z0),
  R* = min Re w over all such roots (conjugates give conjugate w).
Exact method per graph:
  * factor psi over Q; isolate every root in a rational rectangle (sympy Poly.intervals, exact);
  * select the roots of minimal real part by rectangle comparison; ties between non-conjugate roots are
    certified equal only when their real-part intervals sit inside an interval that isolates a single real root of
    the real-part resultant  P_f(x) = Res_z(f(z), f(2x - z))  (product over the factors involved);
  * w = p(z0) with p = chi1 * (psi1')^{-1} mod f in Q[z];  2 Re w is a root of T(y) = Res_u(C(u), C(y - u)),
    C(u) = Res_z(f(z), u - p(z)) (for real z0: Re w = p(z0) is a root of C); factor, and pick the unique
    irreducible factor with exactly one real root in a rigorous complex-interval enclosure of 2 Re w;
  * sign: minimal polynomial = x  <=>  R* = 0;  otherwise refine the isolating interval until it excludes 0.
  * positivity at large t, beta = beta_c:  y = adj(I + beta N) 1 must be strictly of one sign
    (x = K^{-1} 1 = y / (det (1 + t sigma)), sigma = 1^T (I+beta N)^{-1} 1); evaluated exactly at beta_c.
Usage: python3 rstar_exact.py [jsonl] [out.json]
"""
import sys, json
from fractions import Fraction as F
import sympy as sp

z, u, y, x, b = sp.symbols('z u y x b')
QQ = sp.QQ


def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


# ---------------- exact interval arithmetic (Fractions) ----------------
def imul(a, b_):
    ps = [a[0] * b_[0], a[0] * b_[1], a[1] * b_[0], a[1] * b_[1]]
    return (min(ps), max(ps))
def iadd(a, b_): return (a[0] + b_[0], a[1] + b_[1])
def isub(a, b_): return (a[0] - b_[1], a[1] - b_[0])

def ceval(coeffs_hi_lo, rect):
    """Rigorous enclosure of p(z) for z in rect = ((xl,xh),(yl,yh)); coefficients rational, Horner."""
    zr, zi = rect
    re, im = (F(0), F(0)), (F(0), F(0))
    for c in coeffs_hi_lo:
        c = F(int(sp.numer(c)), int(sp.denom(c)))
        nre = isub(imul(re, zr), imul(im, zi)); nim = iadd(imul(re, zi), imul(im, zr))
        re, im = iadd(nre, (c, c)), nim
    return re, im

def fr(q): q = sp.Rational(q); return F(int(q.p), int(q.q))


def roots_rects(f, eps):
    """All roots of irreducible sympy Poly f as (rect, is_real, im_sign)."""
    re_l, cx_l = f.intervals(all=True, eps=eps)
    out = []
    for (lo, hi), m in re_l:
        out.append((((fr(lo), fr(hi)), (F(0), F(0))), True, 0))
    for (c1, c2), m in cx_l:
        r1, i1 = sp.re(c1), sp.im(c1); r2, i2 = sp.re(c2), sp.im(c2)
        rect = ((fr(r1), fr(r2)), (fr(i1), fr(i2)))
        out.append((rect, False, 1 if rect[1][0] > 0 else -1))
    return out


def realpart_poly(fs):
    P = sp.Integer(1)
    for f in fs:
        fe = f.as_expr()
        P *= sp.resultant(fe, fe.subs(z, 2 * x - z), z)
    return sp.Poly(P, x)


def sqf_count(P, lo, hi):
    Ps = sp.Poly(sp.sqf_part(P.as_expr()), P.gens[0])
    return Ps.count_roots(sp.Rational(lo.numerator, lo.denominator), sp.Rational(hi.numerator, hi.denominator))


def isolate_in(polys_irred, lo, hi):
    """Return the irreducible factor with a root in [lo,hi] if exactly one root among all factors lies there."""
    hits = []
    for M in polys_irred:
        k = M.count_roots(sp.Rational(lo.numerator, lo.denominator), sp.Rational(hi.numerator, hi.denominator))
        if k: hits.append((M, k))
    if len(hits) == 1 and hits[0][1] == 1:
        return hits[0][0]
    return None


def refine_real_root(M, lo, hi, target_excl=None, width=F(1, 10**12)):
    """Bisect [lo,hi] (containing exactly one root of irreducible M, endpoints not roots) to width; exact."""
    lo, hi = F(lo), F(hi)
    def val(q): return M.eval(sp.Rational(q.numerator, q.denominator))
    while hi - lo > width or (target_excl is not None and lo <= target_excl <= hi):
        mid = (lo + hi) / 2
        vm = val(mid)
        if vm == 0:
            return (mid, mid)
        if M.count_roots(sp.Rational(lo.numerator, lo.denominator), sp.Rational(mid.numerator, mid.denominator)) == 1:
            hi = mid
        else:
            lo = mid
    return (lo, hi)


def simplify_interval(M, lo, hi, avoid=None):
    """Outward-rounded decimal isolating interval for the unique root of irreducible M in [lo,hi] (exact checks)."""
    if lo == hi:
        return lo, hi
    import math
    for k in range(2, 80):
        D = 10**k
        lo2 = F(math.floor(lo * D), D); hi2 = F(math.ceil(hi * D), D)
        if avoid is not None and lo2 <= avoid <= hi2 and not (lo <= avoid <= hi):
            continue
        if M.eval(sp.Rational(lo2.numerator, lo2.denominator)) == 0 or M.eval(sp.Rational(hi2.numerator, hi2.denominator)) == 0:
            continue
        if M.count_roots(sp.Rational(lo2.numerator, lo2.denominator), sp.Rational(hi2.numerator, hi2.denominator)) == 1:
            return lo2, hi2
    return lo, hi


def coeffs_list(P):
    return [int(c) for c in sp.Poly(P.as_expr() * sp.lcm([sp.denom(c) for c in P.all_coeffs()]), P.gens[0]).all_coeffs()]


def analyze(A):
    m = len(A)
    MA = sp.Matrix(A); J = sp.ones(m, m)
    chi = sp.Poly((-MA).charpoly(z).as_expr(), z)              # det(zI + A)
    chiJ = sp.Poly((-(MA + J)).charpoly(z).as_expr(), z)      # det(zI + A + J)
    psi = chiJ - chi
    h = sp.gcd(chi, psi); chi1 = sp.quo(chi, h); psi1 = sp.quo(psi, h)
    _, psifac = sp.factor_list(psi.as_expr(), z)
    _, chifac = sp.factor_list(chi.as_expr(), z)
    chimult = {sp.Poly(f, z).monic(): e for f, e in chifac}
    facs = [(sp.Poly(f, z), e) for f, e in psifac if sp.Poly(f, z).degree() > 0]
    rec = {'chi_hi_lo': [int(c) for c in chi.all_coeffs()], 'psi_hi_lo': [int(c) for c in psi.all_coeffs()],
           'gcd_chi_psi_hi_lo': [int(c) for c in sp.Poly(h, z).monic().all_coeffs()] if sp.Poly(h, z).degree() > 0 else [1],
           'psi_factors': [{'f_hi_lo': [int(c) for c in f.all_coeffs()], 'mult_psi': e,
                            'mult_chi': chimult.get(f.monic(), 0)} for f, e in facs]}
    # ---- min real part selection
    eps = sp.Rational(1, 10**30)
    for attempt in range(4):
        allr = []
        for fi, (f, e) in enumerate(facs):
            for rect, isr, sg in roots_rects(f, eps):
                allr.append((fi, rect, isr, sg))
        minhi = min(r[1][0][1] for r in allr)
        cand = [r for r in allr if r[1][0][0] <= minhi]
        # keep one representative per conjugate class (real roots, or Im > 0)
        reps = [r for r in cand if r[2] or r[3] > 0]
        if len(reps) == 1:
            tie = False; break
        # attempt tie certification
        fis = sorted(set(r[0] for r in reps))
        P = realpart_poly([facs[i][0] for i in fis])
        lo = min(r[1][0][0] for r in reps); hi = max(r[1][0][1] for r in reps)
        if lo < hi and sqf_count(P, lo, hi) == 1 and P.eval(sp.Rational(lo.numerator, lo.denominator)) != 0 \
                and P.eval(sp.Rational(hi.numerator, hi.denominator)) != 0:
            tie = True; break
        eps = eps / 10**20
    else:
        raise RuntimeError('min real part selection not resolved')
    # nu_inf exact description
    r0 = reps[0]; f0 = facs[r0[0]][0]
    if r0[2]:
        Mnu_cands = [sp.Poly(g, x) for g, _ in sp.factor_list(f0.as_expr().subs(z, x), x)[1]]
    else:
        Pf = realpart_poly([f0])
        Mnu_cands = [sp.Poly(g, x) for g, _ in sp.factor_list(Pf.as_expr(), x)[1]]
    nlo, nhi = r0[1][0]
    Mnu = isolate_in(Mnu_cands, nlo, nhi)
    if Mnu is None:
        raise RuntimeError('nu_inf minpoly not isolated')
    nlo, nhi = refine_real_root(Mnu, nlo, nhi, width=F(1, 10**30)) if nlo != nhi else (nlo, nhi)
    if Mnu.degree() == 1:
        nv = fr(-Mnu.all_coeffs()[1] / Mnu.all_coeffs()[0]); nlo = nhi = nv
    else:
        nlo, nhi = simplify_interval(Mnu, nlo, nhi, avoid=F(1, 2))
    rec['nu_inf'] = {'minpoly_hi_to_lo': coeffs_list(Mnu), 'interval': [str(nlo), str(nhi)]}
    Malpha = sp.Poly(Mnu.as_expr().subs(x, -x), x)
    Malpha = sp.Poly(Malpha.as_expr() * (1 if Malpha.LC() > 0 else -1), x)
    rec['alpha'] = {'minpoly_hi_to_lo': coeffs_list(Malpha), 'interval': [str(-nhi), str(-nlo)]}
    # beta_c = 1/(1 - nu) ; minpoly in x: Mnu(1 - 1/x) * x^d
    d = Mnu.degree()
    Mb = sp.Poly(sp.expand(Mnu.as_expr().subs(x, 1 - 1 / x) * x**d), x)
    Mb = sp.Poly(sp.factor_list(Mb.as_expr(), x)[1][0][0] if len(sp.factor_list(Mb.as_expr(), x)[1]) == 1 else Mb.as_expr(), x)
    Mb = sp.Poly(Mb.as_expr() * (1 if Mb.LC() > 0 else -1), x)
    blo, bhi = (1 / (1 - nlo), 1 / (1 - nhi)) if nhi < 1 else (None, None)
    rec['beta_c'] = {'minpoly_hi_to_lo': coeffs_list(Mb), 'interval': [str(blo), str(bhi)]}
    half = F(1, 2)
    if Mnu.as_expr() == x - sp.Rational(1, 2) or Mnu.monic().as_expr() == x - sp.Rational(1, 2):
        rec['nu_vs_half'] = 0
    else:
        rec['nu_vs_half'] = 1 if nlo > half else (-1 if nhi < half else None)
    # ---- per-root data
    roots_out = []
    status = 'ok'
    for r in reps:
        fi, rect, isr, sg = r
        f, e_psi = facs[fi]
        e_chi = chimult.get(f.monic(), 0)
        desc = {'factor_hi_lo': [int(c) for c in f.all_coeffs()], 'real': isr,
                'rect': [[str(rect[0][0]), str(rect[0][1])], [str(rect[1][0]), str(rect[1][1])]],
                'mult_psi': e_psi, 'mult_chi': e_chi}
        if e_chi > 0:
            desc['kind'] = 'fixed (root of gcd(chi,psi)): nu_A(c) <= nu_inf for every c'
            status = 'fixed_critical'; roots_out.append(desc); continue
        if e_psi > 1:
            desc['kind'] = 'multiple moving root: Lemma AS does not apply'
            status = 'multiple'; roots_out.append(desc); continue
        desc['kind'] = 'simple moving root'
        fq = sp.Poly(f.as_expr(), z, domain=QQ)
        dpsi1 = sp.Poly(sp.diff(psi1.as_expr(), z), z, domain=QQ)
        inv = sp.invert(sp.rem(dpsi1, fq), fq)
        p = sp.rem(sp.Poly(chi1.as_expr(), z, domain=QQ) * inv, fq)
        pc = p.all_coeffs()
        C = sp.Poly(sp.resultant(fq.as_expr(), u - p.as_expr(), z), u)
        if isr:
            T = sp.Poly(C.as_expr().subs(u, y), y); scale = 1
        else:
            T = sp.Poly(sp.resultant(C.as_expr(), C.as_expr().subs(u, y - u), u), y); scale = 2
        Tf = [sp.Poly(g, y) for g, _ in sp.factor_list(T.as_expr(), y)[1]]
        eps2 = sp.Rational(1, 10**30); rr = rect
        for _ in range(6):
            re_enc, _im = ceval(pc, rr)
            lo, hi = scale * re_enc[0], scale * re_enc[1]
            M = isolate_in(Tf, lo, hi) if lo < hi else None
            if M is not None:
                break
            eps2 = eps2 / 10**30
            # re-isolate the same root more tightly
            cands = [q for q in roots_rects(f, eps2) if q[0][0][0] >= rr[0][0] - F(1, 10**25) and q[0][0][1] <= rr[0][1] + F(1, 10**25)
                     and q[0][1][0] >= rr[1][0] - F(1, 10**25) and q[0][1][1] <= rr[1][1] + F(1, 10**25)]
            assert len(cands) == 1, 'root tracking failed'
            rr = cands[0][0]
        else:
            raise RuntimeError('R* minpoly not isolated')
        # minimal polynomial of R* = y/scale
        MR = sp.Poly(sp.expand(M.as_expr().subs(y, scale * x)), x)
        MR = sp.Poly(MR.as_expr() * (1 if MR.LC() > 0 else -1), x)
        Rlo, Rhi = lo / scale, hi / scale
        if MR.degree() == 1:
            val = -MR.all_coeffs()[1] / MR.all_coeffs()[0]
            desc['Rstar_exact'] = str(val); sign = int(sp.sign(val)); Rlo = Rhi = fr(val)
        else:
            Rlo, Rhi = refine_real_root(MR, Rlo, Rhi, target_excl=F(0), width=F(1, 10**20))
            Rlo, Rhi = simplify_interval(MR, Rlo, Rhi, avoid=F(0))
            sign = 1 if Rlo > 0 else -1
            desc['Rstar_exact'] = None
        desc['Rstar_minpoly_hi_to_lo'] = coeffs_list(MR)
        desc['Rstar_interval'] = [str(Rlo), str(Rhi)]
        desc['Rstar_sign'] = sign
        desc['_Rlo'] = Rlo; desc['_Rhi'] = Rhi
        roots_out.append(desc)
    rec['tie_between_distinct_conjugate_classes'] = tie
    rec['critical_roots'] = roots_out
    rec['status'] = status
    if status == 'ok':
        # R* = min over critical roots; decide which is min exactly when several
        best = None
        for d_ in roots_out:
            if best is None or d_['_Rhi'] < best['_Rlo']:
                best = d_
            elif not (d_['_Rlo'] > best['_Rhi']):
                if d_['Rstar_minpoly_hi_to_lo'] != best['Rstar_minpoly_hi_to_lo'] or d_['Rstar_interval'] != best['Rstar_interval']:
                    raise RuntimeError('cannot order R* values of tied roots')
        rec['Rstar'] = {'minpoly_hi_to_lo': best['Rstar_minpoly_hi_to_lo'], 'interval': best['Rstar_interval'],
                        'exact': best['Rstar_exact'], 'sign': best['Rstar_sign']}
        rec['Rstar_vs_1/8'] = (1 if best['_Rlo'] > F(1, 8) else -1 if best['_Rhi'] < F(1, 8) else
                               (0 if best['Rstar_exact'] == '1/8' else None))
        lo_, hi_ = best['_Rlo'], best['_Rhi']
        if blo is not None:
            rec['beta_c3_Rstar_interval_rounded'] = [float(blo**3 * lo_), float(bhi**3 * hi_)] if lo_ > 0 else [float(bhi**3 * lo_), float(blo**3 * hi_)]
    for d_ in roots_out:
        d_.pop('_Rlo', None); d_.pop('_Rhi', None)
    # ---- positivity at beta_c (large t): adj(I + beta N) 1 strictly one-signed
    N = sp.Matrix(m, m, lambda i, j: 1 if (i != j and A[j][i] == 0) else 0)
    Mb_ = sp.eye(m) + b * N
    v = Mb_.adjugate() * sp.ones(m, 1)
    signs = []
    for i in range(m):
        wi = sp.Poly(sp.expand(sp.together(v[i].subs(b, 1 / (1 - x)) * (1 - x)**m)), x)
        if sp.rem(wi, Mnu).is_zero:
            signs.append(0); continue
        lo_, hi_ = nlo, nhi
        while True:
            # interval eval of wi on [lo_,hi_]
            acc = (F(0), F(0))
            for c in wi.all_coeffs():
                c = fr(c); acc = iadd(imul(acc, (lo_, hi_)), (c, c))
            if acc[0] > 0: signs.append(1); break
            if acc[1] < 0: signs.append(-1); break
            lo_, hi_ = refine_real_root(Mnu, lo_, hi_, width=(hi_ - lo_) / 2**40) if lo_ != hi_ else (lo_, hi_)
    rec['positivity_signs_at_beta_c'] = signs
    rec['positive_at_beta_c'] = (all(s_ == 1 for s_ in signs) or all(s_ == -1 for s_ in signs))
    return rec


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else '/home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl'
    out = sys.argv[2] if len(sys.argv) > 2 else 'rstar_raw.json'
    only = sys.argv[3].split(',') if len(sys.argv) > 3 else None
    L = [json.loads(l) for l in open(src)]
    res = []
    for k, r in enumerate(L):
        if only and r['d6'] not in only: continue
        rec = analyze(from_d6(r['d6']))
        rec['d6'] = r['d6']
        res.append(rec)
        print(k, r['d6'], rec['status'], rec.get('Rstar', {}).get('sign'), rec.get('Rstar', {}).get('exact'),
              rec['positive_at_beta_c'], r['bounded_in_t_somewhere'], flush=True)
    json.dump(res, open(out, 'w'), indent=1)

if __name__ == '__main__':
    main()
