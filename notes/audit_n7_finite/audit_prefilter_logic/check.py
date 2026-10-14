"""Independent exact audit of filter7p (via instrumented harness) and resolve_U semantics.
Input: lines '<d6> <domfree> | coefs_i ... | E evals ... | FLAG' from ./harness with p-list on argv.
For each graph:
  * domination (own reimplementation of Def 1.1)            vs harness domfree
  * a_i(beta) by sympy Berkowitz det of symbolic C[i<-1]     vs harness m!*a_i(2+u) coefficients (exact equality)
  * 64^m a_i(p/64)                                           vs harness a_eval (exact equality)
  * N convention end-to-end: x solving (I-W)x=1 built from CLAUDE.md W at a random (eps,delta) with beta=p/64
    must be proportional to a(p/64)
  * ground truth G := exists beta>2 with all a_i(beta)>0, by own Sturm-sequence cell decomposition (no sympy isolation)
  * verdict consistency: D => dominated; N => G False and some i has all coefs <= 0; P => G True; U => record G,
    and fastregion.region(..., INF, positivity_only=True) must equal G.
  * Python bigint simulation of det128/descartes recording max |intermediate| (overflow headroom).
"""
import sys, random, math, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import sympy as sp

b, u = sp.symbols('beta u')

def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0]*n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A

def dominated(A):
    n = len(A)
    for j in range(n):
        for k in range(n):
            if j == k: continue
            # Def 1.1 (ii): j -> k and k -/-> j ; (i): every i not in {j,k} with i -> j has i -> k
            if A[j][k] == 1 and A[k][j] == 0 and all(not (A[i][j] == 1 and A[i][k] == 0) for i in range(n) if i not in (j, k)):
                return (k, j)
    return None

def a_polys(A):
    n = len(A)
    C = sp.Matrix(n, n, lambda i, j: 1 if i == j else (b if A[j][i] == 0 else 0))
    out = []
    for i in range(n):
        M = C.copy()
        for r in range(n): M[r, i] = 1
        out.append(sp.Poly(sp.expand(M.det(method='berkowitz')), b, domain='ZZ'))
    return out

# ---------------- own Sturm machinery (Fractions, low->high lists)
def trim(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p

def ev(p, x):
    r = F(0)
    for c in reversed(p): r = r * x + c
    return r

def prem_neg(a, bb):
    a = [F(x) for x in a]; bb = trim([F(x) for x in bb])
    while len(a) >= len(bb) and any(a):
        a = trim(a)
        if len(a) < len(bb): break
        c = a[-1] / bb[-1]; sh = len(a) - len(bb)
        for i, v in enumerate(bb): a[sh + i] -= c * v
        a = trim(a)
    r = [-x for x in trim(a)]
    if not r: return r
    den = 1
    for x in r: den = den * x.denominator // math.gcd(den, x.denominator)
    ints = [int(x * den) for x in r]
    g = 0
    for x in ints: g = math.gcd(g, abs(x))
    return [x // g for x in ints]          # positive rescaling only

def sturm_seq(p):
    seq = [p, trim([i * p[i] for i in range(1, len(p))])]
    while True:
        r = prem_neg(seq[-2], seq[-1])
        if not r: break
        seq.append(r)
    return seq

def V(seq, x):
    s = [v for v in (ev(q, x) for q in seq) if v != 0]
    return sum(1 for i in range(len(s) - 1) if (s[i] > 0) != (s[i + 1] > 0))

def ground_truth(apolys):
    if any(p.is_zero for p in apolys): return False, 0
    P = sp.Poly(1, b, domain='ZZ')
    for p in apolys:
        if p.degree() > 0:
            P = sp.lcm(P, sp.sqf_part(p))
    while P.degree() > 0 and P.eval(2) == 0:
        P = sp.Poly(sp.quo(P.as_expr(), b - 2), b, domain='ZZ')
    al = [[F(int(c)) for c in reversed(p.all_coeffs())] for p in apolys]
    S = []
    if P.degree() <= 0:
        S = [F(3)]
    else:
        pl = [int(c) for c in reversed(P.all_coeffs())]
        R = 1 + max(abs(F(c, pl[-1])) for c in pl[:-1])
        R = F(math.ceil(R) + 1)
        seq = sturm_seq(pl)
        cnt = lambda lo, hi: V(seq, lo) - V(seq, hi)
        assert ev(pl, F(2)) != 0 and ev(pl, R) != 0
        final = []; stack = [(F(2), R)]
        while stack:
            lo, hi = stack.pop()
            c0 = cnt(lo, hi)
            if c0 <= 1:
                final.append((lo, hi, c0)); continue
            c = (lo + hi) / 2; k = 2
            while ev(pl, c) == 0:
                k += 1; c = lo + (hi - lo) / k
            stack.append((lo, c)); stack.append((c, hi))
        final.sort()
        S = [hi for (lo, hi, c0) in final]
        lo, hi, c0 = final[0]
        if c0 == 0:
            S.append((lo + hi) / 2)
        else:
            c = hi
            while True:
                c = (F(2) + c) / 2
                if ev(pl, c) != 0 and cnt(F(2), c) == 0: break
            S.append(c)
        assert all(s > 2 for s in S)
    for s in S:
        if all(ev(p, s) > 0 for p in al):
            return True, len(S)
    return False, len(S)

# ---------------- bigint simulation of det128/descartes to record max |intermediate|
MAXABS = [0]
def note(*vs):
    for v in vs:
        if abs(v) > MAXABS[0]: MAXABS[0] = abs(v)

def det_sim(M):
    M = [r[:] for r in M]; k = len(M); prev = 1; sign = 1
    for c in range(k - 1):
        if M[c][c] == 0:
            sw = next((r for r in range(c + 1, k) if M[r][c] != 0), None)
            if sw is None: return 0
            M[c], M[sw] = M[sw], M[c]; sign = -sign
        for r in range(c + 1, k):
            for j in range(c + 1, k):
                x1 = M[r][j] * M[c][c]; x2 = M[r][c] * M[c][j]; note(x1, x2, x1 - x2)
                q, rem = divmod(x1 - x2, prev); assert rem == 0
                M[r][j] = q
        prev = M[c][c]
    return sign * M[k - 1][k - 1]

def aeval_sim(Nm, i, p, q):
    n = len(Nm)
    M = [[(q if r == c else 0) + p * Nm[r][c] for c in range(n)] for r in range(n)]
    for r in range(n): M[r][i] = q
    return det_sim(M)

def descartes_sim(Nm, i):
    m = len(Nm)
    f = [aeval_sim(Nm, i, 2 + k, 1) for k in range(m + 2)]
    w = f[:m + 1]; D = []
    for k in range(m + 1):
        D.append(w[0])
        for j in range(m - k): w[j] = w[j + 1] - w[j]; note(w[j])
    fact = math.factorial(m); coef = [0] * (m + 1); FF = [0] * (m + 2); FF[0] = 1
    for k in range(m + 1):
        if k > 0:
            for j in range(k, 0, -1): FF[j] = FF[j - 1] - (k - 1) * FF[j]
            FF[0] = -(k - 1) * FF[0]
        sc = D[k] * (fact // math.factorial(k)); note(sc)
        for j in range(k + 1):
            note(sc * FF[j]); coef[j] += sc * FF[j]; note(coef[j])
    val = 0
    for j in range(m, -1, -1):
        val = val * (m + 1); note(val); val += coef[j]; note(val)
    assert val == fact * f[m + 1]
    note(fact * f[m + 1])
    return coef

def main():
    ps = [int(x) for x in sys.argv[1].split(',')]
    regcheck = '--region' in sys.argv
    if regcheck:
        from fastregion import region, INF
    rng = random.Random(12345)
    tally = {}; bad = 0; nlines = 0
    for line in sys.stdin:
        nlines += 1
        parts = line.strip().split(' | ')
        d6, df = parts[0].split(); df = int(df)
        flag = parts[-1].strip()
        coef_c = [[int(x) for x in seg.split()] for seg in parts[1:-2]]
        evals = [int(x) for x in parts[-2].split()[1:]]
        A = from_d6(d6); n = len(A)
        errs = []
        dom = dominated(A)
        if (dom is None) != bool(df): errs.append('domfree mismatch')
        ap = a_polys(A)
        # exact coefficients of m! a_i(2+u)
        for i in range(n):
            e = sp.Poly(sp.expand(math.factorial(n) * ap[i].as_expr().subs(b, 2 + u)), u)
            cl = [int(e.coeff_monomial(u**j)) for j in range(n + 1)]
            if cl != coef_c[i]: errs.append(f'coef mismatch i={i}')
        # exact a_eval
        k = 0
        for p in ps:
            for i in range(n):
                exp = F(64) ** n * F(int(ap[i].eval(sp.Rational(p, 64))))if False else int(sp.Rational(64) ** n * ap[i].eval(sp.Rational(p, 64)))
                if exp != evals[k]: errs.append(f'a_eval mismatch i={i} p={p}')
                k += 1
        # N convention end-to-end through CLAUDE.md W
        p = rng.choice(ps); beta = F(p, 64); eps = F(rng.randint(1, 63), 64); delta = (beta - 1) * eps
        IW = [[F(1) if i == j else (1 - eps if A[j][i] == 1 else 1 + delta) for j in range(n)] for i in range(n)]
        Mx = sp.Matrix(n, n, lambda i, j: sp.Rational(IW[i][j].numerator, IW[i][j].denominator))
        if Mx.det() != 0:
            x = Mx.LUsolve(sp.ones(n, 1))
            av = [sp.Rational(64) ** n * ap[i].eval(sp.Rational(p, 64)) for i in range(n)]
            if any(x[i] * av[j] != x[j] * av[i] for i in range(n) for j in range(n)):
                errs.append('x not proportional to a: N convention')
        # bigint simulation
        Nm = [[0 if i == j else 1 - A[j][i] for j in range(n)] for i in range(n)]
        for i in range(n):
            if descartes_sim(Nm, i) != coef_c[i]: errs.append('sim coef mismatch')
        for (pp, qq) in [(1,1),(1,2),(1,4),(1,8),(2,1),(1,16),(4,1),(1,32),(8,1),(16,1),(48,1),(1,64),(3,8),(3,4)]:
            for i in range(n): aeval_sim(Nm, i, 2 * qq + pp, qq)
        # verdict
        G, nS = ground_truth(ap)
        if flag == 'D':
            if dom is None: errs.append('D without domination')
            if G: errs.append('dominated but positivity exists: contradicts Lemma D')
        elif flag == 'N':
            if G: errs.append('N-DROP BUT POSITIVITY EXISTS')
            if not any(all(c <= 0 for c in coef_c[i]) for i in range(n)): errs.append('N without nonpositive row')
        elif flag == 'P':
            if not G: errs.append('P but ground truth False')
        elif flag == 'U':
            if regcheck:
                c = region(A, list(range(n)), INF, positivity_only=True, replay=True)
                rv = c['verdict'] == 'positivity_possible'
                if rv != G: errs.append(f'region {c["verdict"]} vs ground truth {G}')
        # domination semantic: dominated graph never has positive a at legal beta? (Lemma D, record only)
        key = flag + ('/G' if G else '/-')
        tally[key] = tally.get(key, 0) + 1
        if errs:
            bad += 1; print('ERR', d6, errs, flush=True)
    print(json.dumps({'lines': nlines, 'bad': bad, 'tally': tally, 'max_abs_intermediate': str(MAXABS[0]),
                      'log2_max': MAXABS[0].bit_length()}))

if __name__ == '__main__':
    main()
