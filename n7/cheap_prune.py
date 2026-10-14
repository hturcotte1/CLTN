"""
cheap_prune.py -- exact, sound NECESSARY-condition pass before the two-parameter CAD.

A full support sigma = V carries a stable positive fixed point at legal (beta, t) only if
  (i)   all a_i(beta) > 0                     (Lemma A.2; a = adj(I + beta N) 1 is t-free, Lemma A)
  (ii)  t > 1/(beta - 2)                       (legality, q > 1)
  (iii) every coefficient of det(mu I + K) = phi(mu) + t*omega(mu) is > 0
        (Stodola: a Hurwitz polynomial has all coefficients positive; K = I + beta N + tJ,
         phi = det(mu I + C), omega = 1^T adj(mu I + C) 1, matrix determinant lemma; the
         coefficient of mu^(m-k) is phi_k(beta) + t*omega_k(beta), exactly as
         fastregion.stability_coeffs).
For fixed beta the constraints (ii)-(iii) are linear in t, so the feasible t-set is an open
interval (max of lower bounds, min of upper bounds), empty if some omega_k(beta) = 0 with
phi_k(beta) <= 0. The feasible set in the (beta, t) plane is OPEN (strict inequalities of
continuous functions), and its projection to beta is a finite union of open intervals whose
endpoints are roots of: a_i, phi_k, omega_k, the pairwise numerators phi_k*omega_l -
phi_l*omega_k, and the legality numerators (beta-2)*phi_k + omega_k. Between consecutive
roots the ordering of all bounds is constant, so one exact rational sample per root-free
open cell (fastregion._cells on the exact isolating intervals) decides the cell, and the
whole set is nonempty iff some cell sample is feasible.

Verdict 'cheap_never' is a certificate that NO legal (beta, t) satisfies (i)-(iii), hence
no stable positive full-support fixed point exists (sound prune). 'survives' means the CAD
must decide. --selfcheck additionally tests every cheap_never graph at random exact legal
points (must all be infeasible) and every survivor at its feasible sample (must be feasible).
Lines in: "<digraph6> ..." ; lines out: "<digraph6> <cheap_never|survives> <n_critical>".
No floating point anywhere.
"""
import sys, os, json, time, argparse, random
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fastregion import beta_data, peval, padd, pscale, _cells
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy import ZZ


def hi_to_lo(p):
    q = list(p)
    while q and q[-1] == 0:
        q.pop()
    return [int(c) for c in reversed(q)] if q else []


def pmul(p, q):
    if not p or not q:
        return []
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                r[i + j] += a * b
    return r


def psub(p, q):
    return padd(p, pscale(q, -1))


def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A


def t_interval(bd, b):
    """Exact feasible t-interval (lo, hi) at rational beta b for (ii)-(iii); hi None = +inf; returns None if empty."""
    m = bd["m"]
    lo = 1 / (b - 2); hi = None
    for k in range(1, m + 1):
        ph = F(peval(bd["phi"][k], b)); om = F(peval(bd["omega"][k], b))
        if om > 0:
            lo = max(lo, -ph / om)
        elif om < 0:
            ub = -ph / om; hi = ub if hi is None else min(hi, ub)
        elif ph <= 0:
            return None
    if hi is not None and not (lo < hi):
        return None
    return (lo, hi)


def feasible_at(bd, b, t):
    """Exact pointwise test of (i)-(iii) at legal (b, t)."""
    if any(peval(ai, b) <= 0 for ai in bd["a"]):
        return False
    if not (t * (b - 2) > 1):
        return False
    for k in range(1, bd["m"] + 1):
        if F(peval(bd["phi"][k], b)) + t * F(peval(bd["omega"][k], b)) <= 0:
            return False
    return True


def cheap_pass(A, bd=None):
    m = len(A)
    bd = bd or beta_data(A, list(range(m)))
    a, phi, omega = bd["a"], bd["phi"], bd["omega"]
    ks = list(range(1, m + 1))
    polys = list(a) + [phi[k] for k in ks] + [omega[k] for k in ks]
    for i in ks:
        for j in ks:
            if i < j:
                polys.append(psub(pmul(phi[i], omega[j]), pmul(phi[j], omega[i])))
    bm2 = [-2, 1]
    for k in ks:
        polys.append(padd(pmul(bm2, phi[k]), omega[k]))
    polys = [hi_to_lo(p) for p in polys if any(p)]
    b_min = F(2)
    iso = dup_isolate_real_roots_list(polys, ZZ, inf=2, strict=True, basis=True)
    iso = [iv for iv in iso if not (iv[0][0] == iv[0][1] == b_min)]
    for cell in _cells(iso, b_min):
        b = cell[4]
        if any(peval(ai, b) <= 0 for ai in a):
            continue
        iv = t_interval(bd, b)
        if iv is not None:
            return "survives", len(iso), (b, iv)
    return "cheap_never", len(iso), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--selfcheck", type=int, default=0, help="random exact legal points per cheap_never graph")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    k, K = map(int, args.shard.split("/"))
    rng = random.Random(args.seed)
    t0 = time.time(); st = {"graphs": 0, "survives": 0, "cheap_never": 0, "selfcheck_points": 0, "selfcheck_failures": 0}
    done = set(); mode = "w"
    if args.resume and os.path.exists(args.out):
        good = []
        for line in open(args.out):
            parts = line.split()
            if len(parts) >= 3 and parts[1] in ("survives", "cheap_never"):
                good.append(line if line.endswith("\n") else line + "\n"); done.add(parts[0])
                st["graphs"] += 1; st[parts[1]] += 1
        open(args.out, "w").writelines(good); mode = "a"; st["resumed"] = len(good)
    with open(args.out, mode) as fh:
        for li, line in enumerate(open(args.src)):
            if li % K != k:
                continue
            s = line.split()[0]
            if s in done:
                continue
            A = from_d6(s); m = len(A)
            bd = beta_data(A, list(range(m)))
            v, ni, wit = cheap_pass(A, bd)
            st["graphs"] += 1; st[v] += 1
            if args.selfcheck:
                if v == "survives":
                    b, (lo, hi) = wit
                    t = lo + 1 if hi is None else (lo + hi) / 2
                    st["selfcheck_points"] += 1
                    if not feasible_at(bd, b, t):
                        st["selfcheck_failures"] += 1; print("SELFCHECK FAIL (survivor sample infeasible)", s, flush=True)
                else:
                    for _ in range(args.selfcheck):
                        b = F(2) + F(rng.randint(1, 4000), rng.choice([1000, 999, 512, 97]))
                        t = 1 / (b - 2) + F(rng.randint(1, 100000), rng.choice([1, 7, 100, 1000]))
                        st["selfcheck_points"] += 1
                        if feasible_at(bd, b, t):
                            st["selfcheck_failures"] += 1; print("SELFCHECK FAIL (cheap_never but feasible point)", s, b, t, flush=True)
            fh.write(f"{s} {v} {ni}\n")
            if st["graphs"] % 2000 == 0:
                print(st["graphs"], st["survives"], round(time.time() - t0, 1), flush=True); fh.flush()
    st["seconds"] = round(time.time() - t0, 1)
    json.dump(st, open(args.out + ".stats.json", "w"), indent=1)
    print(json.dumps(st))


if __name__ == "__main__":
    main()
