"""SANITY for notes/lemma_Bprime_2d.md: independent re-verification of the algebra region2d.projection
relies on (the trust base of §4), on a sample of (graph, support) cases. Exact integer arithmetic only.

For each case:
  N  : a_i(0) = 1, o_k(0) = 1, D(0,0) = 1, Delta_{m-1}(0,0) > 0                        (Lemma N)
  H  : region2d's Delta_{m-1} (sympy Berkowitz) == Bareiss det of the Hurwitz matrix built from
       c_l = phi_l + t omega_l, on a grid large enough to force polynomial identity
  RES: every resultant res_t(f,g) in the projection == Bareiss det of the Sylvester matrix (formal
       degrees) at deg-bound+1 integer beta points  (=> identity of polynomials)
  DSC: res_t(f, f_t) == (-1)^{d(d-1)/2} lc_t(f) disc_t(f) at enough points (sympy's disc convention)
  COV: each of D, Delta_{m-1}, L becomes a nonzero constant after dividing out (exactly) by all curves
       and all projection polynomials -> Z(h) ⊂ U Z(f) ∪ Z(P) x R  (covering hypothesis of B''')
  LC : no lc / disc / res in the projection is the zero polynomial
"""
import sys, json, random, itertools
sys.path.insert(0, "/home/user/CLTN")
from fractions import Fraction as F
import sympy as sp
import region2d, fastregion
from fastregion import beta_data, idet, peval
from graphs import from_key
from enum_cache import key_from_int

B, T = region2d.B, region2d.T
SIGNS = []

def coeffs_int(poly_in_T_at_b, d):
    """integer coefficient list (hi->lo) of length d+1 (formal degree d)"""
    c = [poly_in_T_at_b.coeff(T, k) for k in range(d, -1, -1)]
    return [int(v) for v in c]

def sylvester_det(fc, gc):
    p, q = len(fc) - 1, len(gc) - 1
    n = p + q
    M = []
    for i in range(q):
        M.append([0] * i + fc + [0] * (n - p - 1 - i))
    for i in range(p):
        M.append([0] * i + gc + [0] * (n - q - 1 - i))
    return idet(M)

def check_case(A, sigma, tag):
    out = {"tag": tag, "sigma": list(sigma)}
    bd = beta_data(A, sigma)
    m = bd["m"]
    # N
    assert all(peval(a, 0) == 1 for a in bd["a"]), "N: a_i(0) != 1"
    assert all(peval(o, 0) == 1 for o in bd["o"].values()), "N: o_k(0) != 1"
    curves0 = region2d.curves_lc(bd)
    Dp, Dl, Lp = [f for _, f in curves0]
    assert Dp.as_expr().subs({B: 0, T: 0}) == 1, "N: D(0,0) != 1"
    assert Dl.as_expr().subs({B: 0, T: 0}) > 0, "N: Delta(0,0) <= 0"
    # H : independent Delta_{m-1}
    k = m - 1
    if k >= 1:
        degb = max(Dl.degree(B), 0); degt = max(Dl.degree(T), 0)
        nb = max(degb, k * m) + 1; nt = max(degt, k) + 1
        for bi in range(nb):
            for ti in range(nt):
                c = [peval(bd["phi"][l], bi) + ti * peval(bd["omega"][l], bi) for l in range(m + 1)]
                H = [[(c[2 * i - j] if 0 <= 2 * i - j <= m else 0) for j in range(1, k + 1)] for i in range(1, k + 1)]
                assert idet([[int(x) for x in r] for r in H]) == Dl.as_expr().subs({B: bi, T: ti}), "H mismatch"
    out["H"] = "ok"
    curves, polys, names = region2d.projection(bd, curves0)
    # LC
    for p, nm in zip(polys, names):
        assert not p.is_zero, f"zero projection poly {nm}"
    # RES and DSC
    nres = ndsc = 0
    cur = [(nm, sp.Poly(f.as_expr(), B, T)) for nm, f in curves]
    for (n1, f), (n2, g) in itertools.combinations(cur, 2):
        p_, q_ = f.degree(T), g.degree(T)
        if p_ < 1 or q_ < 1:
            continue
        r = sp.Poly(sp.resultant(f.as_expr(), g.as_expr(), T), B)
        bound = q_ * f.degree(B) + p_ * g.degree(B)
        npts = max(bound, r.degree()) + 1
        sg = None
        for b in range(3, 3 + npts):
            fc = coeffs_int(sp.expand(f.as_expr().subs(B, b)), p_)
            gc = coeffs_int(sp.expand(g.as_expr().subs(B, b)), q_)
            sd, rv = sylvester_det(fc, gc), r.eval(b)
            # sympy 1.14 returns -det(Sylvester) when deg f < deg g and deg f * deg g is odd;
            # only the zero set matters, so compare up to ONE global sign, recorded
            if sg is None and rv != 0:
                sg = 1 if sd == rv else -1
            assert sd == (sg or 1) * rv, f"RES mismatch {n1},{n2} at {b}"
        SIGNS.append((p_, q_, sg))
        nres += 1
    for nm, f in cur:
        d = f.degree(T)
        if d < 2:
            continue
        disc = sp.Poly(sp.discriminant(f.as_expr(), T), B)
        lc = sp.Poly(f.as_expr().coeff(T, d), B)
        ft = sp.Poly(sp.diff(f.as_expr(), T), B, T)
        bound = (d - 1) * f.degree(B) + d * f.degree(B)
        npts = max(bound, (lc * disc).degree()) + 1
        sgn = (-1) ** (d * (d - 1) // 2)
        for b in range(3, 3 + npts):
            fc = coeffs_int(sp.expand(f.as_expr().subs(B, b)), d)
            gc = coeffs_int(sp.expand(ft.as_expr().subs(B, b)), d - 1)
            assert sylvester_det(fc, gc) == sgn * lc.eval(b) * disc.eval(b), f"DSC mismatch {nm} at {b}"
        ndsc += 1
    # COV
    pieces = [sp.Poly(f.as_expr(), B, T) for _, f in curves] + [sp.Poly(p.as_expr(), B, T) for p in polys if p.degree() >= 1]
    for nm, h in (("D", Dp), ("Delta", Dl), ("L", Lp)):
        h = sp.Poly(h.as_expr(), B, T)
        changed = True
        while changed and h.total_degree() >= 1:
            changed = False
            for pc in pieces:
                if pc.total_degree() < 1:
                    continue
                qq, rr = sp.div(h, pc)
                if rr.is_zero:
                    h = qq; changed = True
        assert h.total_degree() <= 0 and not h.is_zero, f"COV failed for {nm}: leftover {h}"
    out.update({"m": m, "n_curves": len(curves), "curve_tdeg": [f.degree(T) for _, f in curves],
                "n_proj": len(polys), "res_checked": nres, "disc_checked": ndsc, "names": names})
    return out

def cases():
    G = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
    yield G, list(range(6)), "geneson_full"
    recs = []
    for sh in range(3):
        for line in open(f"/home/user/CLTN/certificates/region2d_n6_shard{sh}of3.jsonl"):
            r = json.loads(line)
            if r.get("n_beta_cells", 0) > 0:
                recs.append((r["idx"], r["key"], r["verdict"], r["n_beta_cells"]))
    hits = [r for r in recs if r[2] != "never"]
    rnd = random.Random(20261013)
    big = sorted([r for r in recs if r[3] >= 8], key=lambda r: -r[3])[:6]
    pick = hits + big + rnd.sample(recs, 300)
    for idx, key, v, nb in pick:
        yield from_key(key_from_int(key, 6), 6), list(range(6)), f"n6_idx{idx}_{v}_nb{nb}"
    for j in range(10):
        idx, key, v, nb = rnd.choice(recs)
        A = from_key(key_from_int(key, 6), 6)
        size = rnd.choice([2, 3, 4, 5])
        yield A, sorted(rnd.sample(range(6), size)), f"n6_idx{idx}_sub{size}"

if __name__ == "__main__":
    res = []
    for A, sigma, tag in cases():
        r = check_case(A, sigma, tag)
        res.append(r)
        print(tag, "m=%d curves=%d tdeg=%s proj=%d res=%d disc=%d" % (r["m"], r["n_curves"], r["curve_tdeg"], r["n_proj"], r["res_checked"], r["disc_checked"]), flush=True)
    json.dump(res, open(__file__.replace(".py", ".json"), "w"), indent=1)
    import collections
    print("resultant sign vs Sylvester, by (deg_t f, deg_t g):", collections.Counter(SIGNS))
    print("ALL PROJECTION CHECKS PASSED on", len(res), "cases")
