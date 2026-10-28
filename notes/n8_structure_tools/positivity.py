"""Exact limit positivity set P = {beta > 2 : w(beta) > 0} of every candidate, w = b adj(M_b) 1 / det M_b,
M_b = b A^T + (b - 1) I.  det M_b (degree <= 8) and adj(M_b) 1 (degree <= 7) are obtained EXACTLY by integer
Lagrange interpolation at b = 0..8 (exact rational evaluation, then interpolation over Q; the result is checked
to be integral and re-checked at b = 9, 10).  Real roots in (2, oo) by sympy's exact real-root isolation;
signs on each sub-interval decided at an exact rational sample point.
Output: positivity.jsonl (per candidate: P as a list of components with exact endpoint data) and positivity.json."""
import json, os, sys
from fractions import Fraction as F
from multiprocessing import Pool
from collections import Counter
import sympy as sp
import lib

HERE = os.path.dirname(os.path.abspath(__file__))
X = sp.Symbol('b')


def interp(pts):
    """Exact Lagrange interpolation; pts = [(x, y)] -> coefficient list hi->lo (Fractions)."""
    n = len(pts)
    coeffs = [F(0)] * n
    for i, (xi, yi) in enumerate(pts):
        num = [F(1)]
        den = F(1)
        for j, (xj, _) in enumerate(pts):
            if j == i: continue
            num = [a - xj * b for a, b in zip(num + [F(0)], [F(0)] + num)]
            den *= (xi - xj)
        for k in range(n):
            coeffs[k] += yi * num[k] / den
    return coeffs


def det_adj_at(A, bv):
    n = len(A)
    M = [[bv * A[j][i] + (bv - 1) * (1 if i == j else 0) for j in range(n)] for i in range(n)]
    # det by fraction Gaussian elimination, adj(M) 1 = det * M^{-1} 1
    a = [[F(x) for x in row] for row in M]
    det = F(1)
    m = [row[:] + [F(1)] for row in a]
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c] != 0), None)
        if p is None:
            return F(0), None
        if p != c:
            m[c], m[p] = m[p], m[c]; det = -det
        det *= m[c][c]
        for r in range(c + 1, n):
            if m[r][c] != 0:
                f = m[r][c] / m[c][c]
                m[r] = [x - f * y for x, y in zip(m[r], m[c])]
    # back substitution
    x = [F(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (m[i][n] - sum(m[i][j] * x[j] for j in range(i + 1, n))) / m[i][i]
    return det, [det * xi for xi in x]


def polys(A):
    pts, vals = [], []
    cand = [F(k) for k in range(0, 40)] + [F(2 * k + 1, 2) for k in range(40)]
    extra = []
    for bv in cand:
        d, adj = det_adj_at(A, bv)
        if adj is None:          # det M_b = 0 at this sample: skip it (finitely many such b)
            continue
        if len(pts) < 9:
            pts.append(bv); vals.append((d, adj))
        else:
            extra.append(bv)
            if len(extra) == 2:
                break
    det = interp([(p, v[0]) for p, v in zip(pts, vals)])
    adj = [interp([(p, v[1][i]) for p, v in zip(pts, vals)]) for i in range(len(A))]
    for c in [det] + adj:
        assert all(x.denominator == 1 for x in c)
    def strip(c):
        c = [int(x) for x in c]
        while len(c) > 1 and c[0] == 0:
            c = c[1:]
        return c
    det = strip(det); adj = [strip(c) for c in adj]
    for bv in extra:
        d, a = det_adj_at(A, bv)
        ev = lambda c: sum(F(x) * bv ** (len(c) - 1 - k) for k, x in enumerate(c))
        assert a is not None
        assert ev(det) == d and all(ev(c) == ai for c, ai in zip(adj, a))
    return det, adj


def possets(d6):
    A = lib.from_d6(d6)
    det, adj = polys(A)
    ps = [sp.Poly(det, X)] + [sp.Poly(c, X) for c in adj]
    crit = []
    for p in ps:
        if p.is_zero:
            continue
        for r in p.real_roots():
            if r > 2:
                crit.append(r)
    crit = sorted(set(crit), key=lambda r: sp.Rational(str(sp.N(r, 60))))
    # strictly increasing check with exact comparison
    for u, v in zip(crit, crit[1:]):
        assert u < v
    edges = [sp.Integer(2)] + crit
    comps = []
    for k in range(len(edges)):
        lo = edges[k]; hi = edges[k + 1] if k + 1 < len(edges) else None
        if hi is None:
            s = sp.floor(sp.N(lo, 30)) + 1
            if not s > lo: s += 1
        else:
            s = (sp.Rational(str(sp.N(lo, 60))) + sp.Rational(str(sp.N(hi, 60)))) / 2
        assert lo < s and (hi is None or s < hi)
        dv = ps[0].eval(s)
        pos = dv != 0 and all(s * p.eval(s) / dv > 0 for p in ps[1:])
        comps.append((str(lo), str(hi) if hi is not None else 'oo', pos))
    # merge
    P = []
    for lo, hi, pos in comps:
        if pos:
            if P and P[-1][1] == lo:
                P[-1][1] = hi
            else:
                P.append([lo, hi])
    return dict(d6=d6, P=P, unbounded=bool(P) and P[-1][1] == 'oo', components=len(P),
                deg_det=len(det) - 1)


if __name__ == '__main__':
    hits = [json.loads(l)['d6'] for l in open(lib.REPO + '/n8/n8_limit_HITS.jsonl')]
    with Pool(2) as p:
        res = p.map(possets, hits, chunksize=20)
    with open(os.path.join(HERE, 'positivity.jsonl'), 'w') as f:
        for r in res:
            f.write(json.dumps(r) + '\n')
    s = dict(total=len(res), errors=sum(1 for r in res if 'error' in r),
             components=dict(Counter(r.get('components') for r in res)),
             unbounded=sum(1 for r in res if r.get('unbounded')))
    json.dump(s, open(os.path.join(HERE, 'positivity.json'), 'w'), indent=1)
    print(s)
