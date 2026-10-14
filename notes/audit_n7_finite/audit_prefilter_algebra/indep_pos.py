"""Independent exact decision of: exists beta > 2 with a_i(beta) > 0 for all i.
Input: harness output lines (m!*a_i(2+u) coefficients, low->high in u).  Own coprime-base +
Descartes-bisection root isolation (no sympy, no floats).  Output: '<d6> POS|NEVER <reason>'."""
import sys
from fractions import Fraction as F
from math import gcd, factorial

def trim(p):
    while p and p[-1] == 0: p.pop()
    return p
def ev(p, x):
    r = 0
    for c in reversed(p): r = r * x + c
    return r
def prim(p):
    p = trim([F(c) for c in p])
    if not p: return []
    den = 1
    for c in p: den = den * c.denominator // gcd(den, c.denominator)
    q = [int(c * den) for c in p]
    g = 0
    for c in q: g = gcd(g, abs(c))
    q = [c // g for c in q]
    if q[-1] < 0: q = [-c for c in q]
    return q
def pdivmod(a, b):
    a = [F(c) for c in a]; b = [F(c) for c in b]; trim(b)
    qd = [F(0)] * max(len(a) - len(b) + 1, 1)
    while len(trim(a)) >= len(b):
        c = a[-1] / b[-1]; s = len(a) - len(b); qd[s] = c
        for i, bv in enumerate(b): a[s + i] -= c * bv
        a.pop()
    return trim(qd), trim(a)
def pgcd(a, b):
    a, b = prim(a), prim(b)
    while b:
        _, r = pdivmod(a, b); a, b = b, prim(r)
    return prim(a)
def pexdiv(a, b):
    q, r = pdivmod(a, b); assert not r; return prim(q)
def deriv(p): return [i * p[i] for i in range(1, len(p))]
def sqf(p):
    g = pgcd(p, deriv(p))
    return prim(p) if len(g) <= 1 else pexdiv(p, g)

def coprime_base(polys):
    base = []
    for h in polys:
        h = prim(h)
        if len(h) <= 1: continue
        newbase = []
        for f in base:
            d = pgcd(f, h)
            if len(d) > 1:
                newbase.append(d)
                f2 = pexdiv(f, d)
                if len(f2) > 1: newbase.append(f2)
                h = pexdiv(h, d)
                # h may still share factors with d (if h had them squared) - h is squarefree, so not.
            else:
                newbase.append(f)
        if len(h) > 1: newbase.append(h)
        base = newbase
    # verify pairwise coprime
    for i in range(len(base)):
        for j in range(i + 1, len(base)):
            assert len(pgcd(base[i], base[j])) <= 1, "coprime base failed"
    return base

def var(coeffs):
    s = [c for c in coeffs if c != 0]
    return sum(1 for i in range(len(s) - 1) if (s[i] > 0) != (s[i + 1] > 0))
def taylor_var(f, l, r):
    """sign variations of (x+1)^d f((l + r x)/(1+x)) -- Descartes bound for roots in (l, r)."""
    d = len(f) - 1
    # g(x) = sum_k f_k (l + r x)^k (1+x)^(d-k)
    out = [F(0)] * (d + 1)
    for k, fk in enumerate(f):
        if fk == 0: continue
        p = [F(1)]
        for _ in range(k):  # * (l + r x)
            p = p + [F(0)]
            p = [p[i] * l + (p[i - 1] * r if i >= 1 else 0) for i in range(len(p))]
        for _ in range(d - k):  # * (1 + x)
            p = p + [F(0)]
            p = [p[i] + (p[i - 1] if i >= 1 else 0) for i in range(len(p))]
        for i, c in enumerate(p): out[i] += fk * c
    return var(out)

def isolate_pos(f):
    """disjoint isolating intervals (l, r) (l<r, f(l)f(r)<0, exactly one root) or exact roots (x, x) in (0, inf)."""
    f = list(f)
    while f[0] == 0: f = f[1:]          # deflate roots at u = 0 (beta = 2, excluded)
    if len(f) <= 1: return []
    B = 1 + max(abs(F(c, f[-1])) for c in f[:-1])
    B = F(int(B) + 1)
    assert ev(f, B) != 0
    out, stack = [], [(F(0), B)]
    while stack:
        l, r = stack.pop()
        v = taylor_var(f, l, r)
        if v == 0: continue
        if v == 1:
            assert ev(f, l) * ev(f, r) < 0 or ev(f, l) == 0
            out.append((l, r)); continue
        m = (l + r) / 2
        if ev(f, m) == 0: out.append((m, m))
        stack.append((l, m)); stack.append((m, r))
    # (l,r) from v==1 with l==0: f(0)!=0 after deflation; endpoints that are roots were recorded exactly and
    # an interval (l, r) with f(l)==0 cannot have v==1 for an open-interval root unless another root lies inside;
    # enforce strict sign change for every open interval:
    res = []
    for l, r in out:
        if l == r: res.append((l, r)); continue
        fl, fr = ev(f, l), ev(f, r)
        assert fl != 0 and fr != 0 and (fl > 0) != (fr > 0), ("bad isolating interval", f, l, r)
        res.append((l, r))
    return res

def refine(f, l, r):
    m = (l + r) / 2; fm = ev(f, m)
    if fm == 0: return (m, m)
    return (l, m) if (ev(f, l) > 0) != (fm > 0) else (m, r)

def decide(P):
    if any(len(trim(list(p))) == 0 for p in P): return "NEVER", "identically zero"
    base = coprime_base([sqf(p) for p in P])
    ivs = []
    for f in base:
        for l, r in isolate_pos(f): ivs.append([l, r, f])
    # refine until pairwise strictly disjoint
    while True:
        ivs.sort(key=lambda t: (t[0], t[1]))
        clash = False
        for k in range(len(ivs) - 1):
            a, b = ivs[k], ivs[k + 1]
            if a[1] >= b[0]:
                clash = True
                for t in (a, b):
                    if t[0] != t[1]: t[0], t[1] = refine(t[2], t[0], t[1])
                if a[0] == a[1] and b[0] == b[1] and a[0] == b[0]:
                    raise AssertionError("two coprime factors with a common root")
                break
        if not clash: break
    pts = []
    prev = F(0)
    for l, r, f in ivs:
        pts.append((prev + l) / 2); prev = r
    pts.append(prev + 1)
    for x in pts:
        if all(ev(p, x) > 0 for p in P): return "POS", str(x + 2)
    return "NEVER", f"{len(ivs)} roots, {len(pts)} cells"

def main():
    for line in sys.stdin:
        parts = line.split('#')
        d6 = parts[0].split()[0]
        n = ord(d6[1]) - 63
        P = []
        mf = factorial(n)
        for i in range(n):
            seg = parts[1 + i].split(':')[1].split()[:-1]
            c = [int(x) for x in seg]
            assert all(x % mf == 0 for x in c)
            P.append(trim([x // mf for x in c]))
        v, why = decide(P)
        print(d6, v, why, parts[0].split()[1], parts[1 + n].strip())

if __name__ == "__main__":
    main()
