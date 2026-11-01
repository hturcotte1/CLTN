"""Independent referee checks of the fixes for the external correctness report of 2026-11-01.

Everything here is rebuilt from the definitions (CLAUDE.md section 2 and paper Section 2) with sympy exact
arithmetic and fractions.Fraction; the project engines are NOT imported, except that the final cross-check of
part 2 builds W from the CLAUDE.md convention by hand.  Floats appear only in 'display' fields.
Run:  cd /home/user/CLTN && nice -n 10 python3 notes/referee_external/fixes_2026-11-01/referee/referee_check.py
Output: referee_check.json next to this script.
"""
import json, os, sys, subprocess, itertools
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = '/home/user/CLTN'
OUT = {}
b, t, s, u, z, mu, c_, e_, d_ = sp.symbols('beta t s u z mu c eps delta')

def check(name, cond, info=None):
    OUT.setdefault('checks', []).append({'name': name, 'pass': bool(cond), 'info': info})
    print(('PASS ' if cond else 'FAIL ') + name + ('' if info is None else f'  [{info}]'), flush=True)
    return cond

# ---------------------------------------------------------------- exact helpers
def hurwitz_dets(coeffs):
    """coeffs c0..cm of c0 mu^m + ... + cm (sympy expressions); returns [Delta_1..Delta_m]."""
    m = len(coeffs) - 1
    def ent(i, j):
        k = 2 * i - j
        return coeffs[k] if 0 <= k <= m else 0
    out = []
    for k in range(1, m + 1):
        H = sp.Matrix(k, k, lambda i, j: ent(i + 1, j + 1))
        out.append(sp.expand(H.det(method='berkowitz')))
    return out

def poly_frac(p, x):
    """sympy poly in x -> list of Fraction coefficients high->low."""
    P = sp.Poly(p, x)
    return [Fr(int(sp.fraction(cc)[0]), int(sp.fraction(cc)[1])) for cc in P.all_coeffs()]

def peval(c, x):
    v = Fr(0)
    for a in c:
        v = v * x + a
    return v

def pdiv_rem(a, bb):
    a = list(a)
    while a and a[0] == 0: a.pop(0)
    while len(a) >= len(bb) and a:
        q = a[0] / bb[0]
        for i in range(len(bb)):
            a[i] -= q * bb[i]
        a.pop(0)
        while a and a[0] == 0: a.pop(0)
    return a

def deriv(c):
    n = len(c) - 1
    return [c[i] * (n - i) for i in range(n)]

def sturm_count(c, lo, hi):
    """number of distinct real roots in (lo, hi] of the polynomial c (Fractions, high->low), own Sturm sequence."""
    seq = [c, deriv(c)]
    while True:
        r = pdiv_rem(seq[-2], seq[-1])
        if not r: break
        seq.append([-x for x in r])
    def V(x):
        vals = [peval(p, x) for p in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))
    return V(lo) - V(hi)

def mobius_positive(p, x, lo, hi):
    """all coefficients of (1+u)^deg p((lo + hi u)/(1+u)) strictly positive => p > 0 on [lo, hi]."""
    P = sp.Poly(p, x)
    dg = P.degree()
    q = sp.expand(sp.cancel(P.as_expr().subs(x, (lo + hi * u) / (1 + u)) * (1 + u) ** dg))
    Q = sp.Poly(q, u)
    co = Q.all_coeffs()
    return len(co) == dg + 1 and all(cc > 0 for cc in co), dg

# ================================================================ PART 1: G_1 and the vertical section
rows = ['011000', '000011', '000111', '110011', '110101', '101110']
A = sp.Matrix([[int(ch) for ch in r] for r in rows])
n = 6
I6, J6 = sp.eye(n), sp.ones(n, n)
N = J6 - I6 - A.T
C = I6 + b * N
K = C + t * J6
# convention check against CLAUDE.md section 2, built entry by entry
W = sp.Matrix(n, n, lambda i, j: 0 if i == j else (-1 + e_ if A[j, i] == 1 else -1 - d_))
lhs = (I6 - W).subs({e_: 1 / (1 + t), d_: (b - 1) / (1 + t)})
check('P1 convention: I - W = eps*K with eps=1/(1+t), delta=(beta-1)/(1+t)', sp.simplify(lhs - K / (1 + t)) == sp.zeros(n, n))

# same graph as Geneson's labelling in the stage-0 certificate (nauty labelg canonical forms)
def to_d6(M):
    m = M.shape[0]; bits = [int(M[i, j]) for i in range(m) for j in range(m)]
    while len(bits) % 6: bits.append(0)
    return '&' + chr(63 + m) + ''.join(chr(63 + int(''.join(map(str, bits[k:k + 6])), 2)) for k in range(0, len(bits), 6))
st0 = json.load(open(os.path.join(ROOT, 'certificates/stage0_geneson_2607.21396.json')))
Ag = sp.Matrix(st0['A_ours'])
def labelg(d6s):
    r = subprocess.run([os.path.join(ROOT, 'n7/labelg'), '-q'], input='\n'.join(d6s) + '\n', capture_output=True, text=True)
    return r.stdout.split()
cf = labelg([to_d6(A), to_d6(Ag)])
check('P1 Table-1 G_1 is isomorphic to the stage-0 Geneson graph (labelg)', cf[0] == cf[1], cf)

p = sp.expand((mu * I6 + K).det(method='berkowitz'))
P = sp.Poly(p, mu)
co = P.all_coeffs()
check('P1 p monic of degree 6', len(co) == 7 and co[0] == 1)
H = hurwitz_dets(co)
D5 = H[4]
a = sp.expand(C.adjugate() * sp.ones(n, 1))
avec = [sp.expand(a[i]) for i in range(n)]
check('P1 a = adj(C)1 does not depend on t', all(sp.Poly(x, b, t).degree(t) <= 0 for x in avec))
degt = sp.Poly(D5, t).degree()
L = t * (b - 2) - 1
Pb = sp.expand(sp.cancel((b - 2) ** degt * D5.subs(t, 1 / (b - 2))))
Res = sp.expand(sp.resultant(D5, L, t))
ratio = sp.cancel(Res / Pb)
check('P1 Res_t(Delta5, L) = constant * (beta-2)^deg_t * Delta5(beta, 1/(beta-2))', ratio.is_number, f'deg_t Delta5 = {degt}, ratio = {ratio}')
Pf = poly_frac(Pb, b)
lo, hi = Fr(17, 8), Fr(15, 7)
g = sp.gcd(sp.Poly(Pb, b), sp.Poly(sp.diff(Pb, b), b))
check('P1 P squarefree', g.degree() == 0)
nr = sturm_count(Pf, lo, hi)
check('P1 exactly one root of P (hence of Res_t(Delta5,L)) in (17/8, 15/7]', nr == 1 and peval(Pf, hi) != 0, f'Sturm count {nr}')
check('P1 cross-check: sympy count_roots on [17/8,15/7] = 1', sp.Poly(Pb, b).count_roots(sp.Rational(17, 8), sp.Rational(15, 7)) == 1)
facs = sp.factor_list(Pb, b)[1]
f14 = None
for fac, mult in facs:
    ff = poly_frac(fac, b)
    if sturm_count(ff, lo, hi) == 1:
        f14 = fac
cert_line = [json.loads(l) for l in open(os.path.join(ROOT, 'certificates/region2d_n6_shard2of3_HITS.jsonl'))]
G1c = [r for r in cert_line if r['idx'] == 1000255][0]
crit = [x for x in G1c['beta_critical'] if x['interval'] == ['17/8', '15/7']][0]
fc = sp.Poly(f14, b).all_coeffs()
check('P1 root lies on an irreducible factor equal to the certificate beta_critical minpoly on [17/8,15/7]',
      [int(x) for x in fc] == crit['minpoly_hi_to_lo'] or [-int(x) for x in fc] == crit['minpoly_hi_to_lo'],
      f'degree {sp.Poly(f14, b).degree()}; factor degrees {[sp.Poly(f_, b).degree() for f_, _ in facs]}')
check('P1 certificate: no recorded beta-cell meets [17/8,15/7]',
      all(not (Fr(cc['beta_cell'][0]) < hi and (cc['beta_cell'][1] is None or Fr(cc['beta_cell'][1]) > lo)) for cc in G1c['cert']['cells']),
      [cc['beta_cell'] for cc in G1c['cert']['cells']])
check('P1 certificate: res_Delta5_f1_L is a projection polynomial name', 'res_Delta5_f1_L' in G1c['projection'])
# refine beta0 by bisection
f14f = poly_frac(f14, b)
L0, H0 = lo, hi
for _ in range(40):
    mid = (L0 + H0) / 2
    if peval(f14f, L0) * peval(f14f, mid) <= 0: H0 = mid
    else: L0 = mid
OUT['beta0_isolating_interval'] = [str(L0), str(H0)]
OUT['beta0_display'] = float(L0)
# (c) positivity on the strip at t = 100 (Geneson's exact check)
okA = [mobius_positive(x, b, sp.Rational(17, 8), sp.Rational(15, 7)) for x in avec]
okH = [mobius_positive(h.subs(t, 100), b, sp.Rational(17, 8), sp.Rational(15, 7)) for h in H]
check('P1 at t=100 all six a_i > 0 on [17/8,15/7] (Mobius coefficient test)', all(o for o, _ in okA), [dg for _, dg in okA])
check('P1 at t=100 all six Hurwitz determinants > 0 on [17/8,15/7]', all(o for o, _ in okH), [dg for _, dg in okH])
check('P1 legality 100(beta0-2) > 1', 100 * (L0 - 2) > 1)

# (d) the whole legal vertical section beta = beta0: substitute t = 1/(beta-2) + s
def sign_at_beta0(gpoly):
    """exact sign of gpoly(beta0): 0 if f14 | gpoly, else refine the isolating interval until gpoly has no root in it."""
    gp = sp.Poly(gpoly, b)
    if gp.is_zero: return 0
    if sp.rem(gp, sp.Poly(f14, b)).is_zero: return 0
    gf = poly_frac(gp.as_expr(), b) if gp.degree() > 0 else None
    if gf is None: return 1 if gp.as_expr() > 0 else -1
    l, h = lo, hi
    while True:
        if peval(gf, l) != 0 and peval(gf, h) != 0 and sturm_count(gf, l, h) == 0:
            v = peval(gf, (l + h) / 2)
            return 1 if v > 0 else -1
        mid = (l + h) / 2
        if peval(f14f, l) * peval(f14f, mid) <= 0: h = mid
        else: l = mid

section = {}
allpos = True
for k, h in enumerate(H, start=1):
    dk = sp.Poly(h, t).degree()
    G = sp.expand(sp.cancel((b - 2) ** dk * h.subs(t, 1 / (b - 2) + s)))
    GP = sp.Poly(G, s)
    signs = [sign_at_beta0(cc) for cc in GP.all_coeffs()]   # high -> low in s
    nz = [x for x in signs if x != 0]
    ok = all(x > 0 for x in nz) and len(nz) > 0
    const = signs[-1]
    section[f'Delta{k}'] = {'deg_t': dk, 'coeff_signs_high_to_low_in_s': signs, 'constant_term_sign': const}
    allpos &= ok
check('P1 on beta=beta0, every Hurwitz det is a polynomial in s=t-1/(beta0-2) with nonnegative coefficients, not all 0 (=> >0 for s>0)', allpos, {k: v['coeff_signs_high_to_low_in_s'] for k, v in section.items()})
check('P1 Delta5 vanishes exactly at the corner (beta0, 1/(beta0-2)) and nowhere else on the closed section', section['Delta5']['constant_term_sign'] == 0 and allpos)
OUT['section_signs'] = section
check('P1 a_i(beta0) > 0', all(sign_at_beta0(x) > 0 for x in avec))
# D = det K > 0 on the section follows from Delta6 = c6*Delta5 > 0; also check directly
Dk = co[6]
GD = sp.Poly(sp.expand(sp.cancel((b - 2) * Dk.subs(t, 1 / (b - 2) + s))), s)
check('P1 D = det K > 0 on the section', all(x > 0 for x in [sign_at_beta0(cc) for cc in GD.all_coeffs()] if x != 0) and sign_at_beta0(GD.all_coeffs()[-1]) > 0)

# ================================================================ PART 2: E_0, Lemma 7.4 hypotheses, R = 1/4
rows7 = ['0011110', '0011100', '1101111', '1110111', '0011011', '0100101', '1100000']
A7 = sp.Matrix([[int(ch) for ch in r] for r in rows7])
check('P2 rows decode to digraph6 &FNFL~[udo?', to_d6(A7) == '&FNFL~[udo?', to_d6(A7))
I7, J7 = sp.eye(7), sp.ones(7, 7)
chi = sp.expand((z * I7 + A7).det(method='berkowitz'))
psi = sp.expand((sp.ones(1, 7) * (z * I7 + A7).adjugate() * sp.ones(7, 1))[0])
mdl = sp.expand((z * I7 + A7 - c_ * J7).det(method='berkowitz'))
check('P2 det(zI-(cJ-A)) = chi - c psi (matrix determinant lemma, exact)', sp.expand(mdl - (chi - c_ * psi)) == 0)
q = sp.Poly(z ** 2 - z + 1, z)
Pp = sp.Poly(psi, z); Pc = sp.Poly(chi, z)
check('P2 z^2-z+1 divides psi', sp.rem(Pp, q).is_zero)
psi2 = sp.quo(Pp, q)
check('P2 z^2-z+1 divides psi exactly once (z0 simple)', not sp.rem(psi2, q).is_zero)
check('P2 gcd(chi, z^2-z+1) = 1 (chi(z0) != 0)', sp.gcd(Pc, q).degree() == 0)
x = sp.symbols('x')
h = sp.Poly(sp.expand(psi2.as_expr().subs(z, sp.Rational(1, 2) - x)), x)
hc = h.all_coeffs()
if hc[0] < 0: hc = [-cc for cc in hc]
Hh = hurwitz_dets(hc)
check('P2 every other root of psi has Re z > 1/2 (Hurwitz dets of psi2(1/2-x) all > 0, leading coeff > 0)', hc[0] > 0 and all(v > 0 for v in Hh), [str(v) for v in Hh])
z0 = sp.Rational(1, 2) - sp.sqrt(3) / 2 * sp.I
val = sp.nsimplify(sp.radsimp(sp.expand(chi.subs(z, z0)) / sp.expand(sp.diff(psi, z).subs(z, z0))))
val = sp.expand(sp.radsimp(val))
re, im = sp.re(val), sp.im(val)
check('P2 chi(z0)/psi\'(z0) = 1/4 - (sqrt3/12) i at z0 = 1/2 - (sqrt3/2) i', sp.simplify(re - sp.Rational(1, 4)) == 0 and sp.simplify(im + sp.sqrt(3) / 12) == 0, str(val))
z0b = sp.conjugate(z0)
valb = sp.expand(sp.radsimp(sp.expand(chi.subs(z, z0b)) / sp.expand(sp.diff(psi, z).subs(z, z0b))))
check('P2 conjugate root gives the same real part', sp.simplify(sp.re(valb) - re) == 0)
shat = psi / chi
v2 = sp.expand(sp.radsimp(1 / sp.expand(sp.radsimp(sp.diff(shat, z).subs(z, z0)))))
check('P2 1/s^\'(z0) equals chi(z0)/psi\'(z0)', sp.simplify(v2 - val) == 0)
R = Fr(1, 4)
check('P2 R = 1/4 > 1/8 (Fraction)', R > Fr(1, 8))
w2 = (A7.T + sp.Rational(1, 2) * I7).inv() * sp.ones(7, 1)
check('P2 (iv) w(2) = (A^T + I/2)^{-1} 1 > 0', all(v > 0 for v in w2), [str(v) for v in w2])
# legality expansion: c - sqrt(c^2-c) = 1/(1+sqrt(1-e)) = (1 - sqrt(1-e))/e, series coefficients all positive
ser = sp.series((1 - sp.sqrt(1 - e_)) / e_, e_, 0, 8).removeO()
sc = [ser.coeff(e_, k) for k in range(8)]
check('P2 c - sqrt(c^2-c) = 1/2 + 1/(8c) + 1/(16c^2) + ... (all coefficients > 0)', sc[0] == sp.Rational(1, 2) and sc[1] == sp.Rational(1, 8) and all(v > 0 for v in sc), [str(v) for v in sc])
# IFT first-order check at c = 1000: gamma between 1/2+1/(8c) and 1/2+R/c stable and legal; above unstable
def hurwitz_frac(M):
    m = M.shape[0]
    cp = sp.Poly(sp.expand((mu * sp.eye(m) + M).det(method='berkowitz')), mu).all_coeffs()
    return [v for v in hurwitz_dets(cp)]
for cval, num, expect in [(1000, 3, True), (1000, 5, False)]:
    gam = sp.Rational(1, 2) + sp.Rational(num, 16 * cval)
    legal = gam ** 2 - 2 * cval * gam + cval < 0
    Hd = hurwitz_frac(cval * J7 - A7.T - gam * I7)
    st = all(v > 0 for v in Hd)
    # map back to (beta, t), eps, delta and check with W built from CLAUDE.md conventions directly
    beta_ = 1 / (1 - gam); t_ = beta_ * (cval - 1)
    eps = 1 / (1 + t_); dl = (beta_ - 1) / (1 + t_)
    W7 = sp.Matrix(7, 7, lambda i, j: 0 if i == j else (-1 + eps if A7[j, i] == 1 else -1 - dl))
    Hw = hurwitz_frac(-(-I7 + W7))  # p(mu) = det(mu I - (-I+W)) Hurwitz <=> -I+W stable
    stw = all(v > 0 for v in Hw)
    xfp = (I7 - W7).LUsolve(sp.ones(7, 1))
    pos = all(v > 0 for v in xfp)
    qv = dl * (1 - eps) / eps
    check(f'P2 c={cval}, gamma=1/2+{num}/(16c): legal={legal}, K stable={st}, W-path stable={stw}, x>0={pos} (expected stable={expect})',
          legal and st == expect and stw == expect and (pos if expect else True), f'eps={eps}, delta={dl}, q={qv}')

# ================================================================ PART 3: primitive proportions, recount with nauty labelg
def from_d6(s_):
    m = ord(s_[1]) - 63; pp = s_[2:]; M = [[0] * m for _ in range(m)]; bit = 0
    for i in range(m):
        for j in range(m):
            M[i][j] = ((ord(pp[bit // 6]) - 63) >> (5 - bit % 6)) & 1; bit += 1
    return M
def d6_of(M):
    return to_d6(sp.Matrix(M))
n6s = json.load(open(os.path.join(ROOT, 'certificates/n6_summary.json')))
h6 = [hh['A'] for hh in n6s['hits'] if sum(map(sum, hh['A'])) < 30]
h7 = [json.loads(l)['d6'] for l in open(os.path.join(ROOT, 'certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl'))]
cls = json.load(open(os.path.join(ROOT, 'certificates/n8_counterexamples_classes.json')))
h8 = [d for row in cls['rows'] for d in row['members_d6']]
check('P3 input sizes 3 / 123 / 9133', (len(h6), len(h7), len(h8)) == (3, 123, 9133), (len(h6), len(h7), len(h8)))
C6 = set(labelg([d6_of(M) for M in h6]))
C7 = set(labelg(h7 + ['&FNFL~[udo?']))
check('P3 canonical classes distinct', len(C6) == 3 and len(C7) == 124)
def subs_d6(d6s, k):
    out, owner = [], []
    for idx, d in enumerate(d6s):
        M = from_d6(d)
        for S in itertools.combinations(range(len(M)), k):
            out.append(d6_of([[M[i][j] for j in S] for i in S])); owner.append(idx)
    return out, owner
s76, o76 = subs_d6(h7, 6)
c76 = labelg(s76)
cont7 = {o for cf_, o in zip(c76, o76) if cf_ in C6}
prim7 = 123 - len(cont7)
s86, o86 = subs_d6(h8, 6)
c86 = labelg(s86)
s87, o87 = subs_d6(h8, 7)
c87 = labelg(s87)
cont8 = {o for cf_, o in zip(c86, o86) if cf_ in C6} | {o for cf_, o in zip(c87, o87) if cf_ in C7}
prim8 = 9133 - len(cont8)
OUT['primitive'] = {'n6': 3, 'n7': prim7, 'n8': prim8, 'n7_containing_n6': len(cont7), 'n8_containing_n6_or_n7(incl FNFL)': len(cont8)}
sh = [Fr(3, 3), Fr(prim7, 123), Fr(prim8, 9133)]
OUT['shares_exact'] = [str(x) for x in sh]
OUT['shares_display'] = [round(float(x), 4) for x in sh]
check('P3 primitive counts 3, 102, 6720 (independent labelg recount)', (prim7, prim8) == (102, 6720), (prim7, prim8))
check('P3 102/123 = 34/41 and 6720/9133 irreducible', Fr(102, 123) == Fr(34, 41) and Fr(6720, 9133).denominator == 9133)
check('P3 shares strictly decrease 1 > 34/41 > 6720/9133', sh[0] > sh[1] > sh[2])
check('P3 display roundings 0.83 and 0.74 correct to 2 d.p.', Fr(825, 1000) <= sh[1] < Fr(835, 1000) and Fr(735, 1000) <= sh[2] < Fr(745, 1000))

OUT['all_pass'] = all(cc['pass'] for cc in OUT['checks'])
json.dump(OUT, open(os.path.join(HERE, 'referee_check.json'), 'w'), indent=1, default=str)
print('ALL PASS' if OUT['all_pass'] else 'SOME FAILED')
