#!/usr/bin/env python3
"""
aggregate_n8.py -- aggregate certificate of the n = 8 eps -> 0 LIMIT map of the target-free clique
conjecture (full supports).  Modelled on n7/aggregate_n7.py and on the writers of
certificates/n7_finite_summary.json and certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit_summary.json.

Inputs (every number in the output comes from one of these; see the 'provenance' map in the output):
  n8/n8_limit_HITS.jsonl                 one record per non-complete limit-stable full support (the 9,133)
  n8/n8_limit_results.jsonl              every strict spectral survivor with its exact limit verdict
  n8/n8_limit_tally.json                 the producer's tally (recounted here from the results file)
  n8/n8_prefilter_TOTAL.txt, n8/n8_accounting.json, n8/n8_restricted_counts.json   the enumeration chain
  n8/check_k8decomp.log, n8/check_orbit.log                                        the two per-B ledgers
  n8/limit_sharding/shard{0..3}_tally.json                                         the sharded limit stage
  certificates/n8_limit_summary.json     the limit-stage certificate (verification blocks are copied from it)
  notes/verify_n8/SYNTHESIS.json, notes/verify_n8_2/SYNTHESIS.json                the two machine campaigns
  notes/verify_n8/v_structure/q_upper_ends_refined.json                            independent q-end refinement
  notes/verify_n8_2/completeness_twins/part2_blowups.jsonl, part34_summary.json    lifts and twin closure
  n7/n7_limit_HITS.jsonl                 the 123 n = 7 limit counterexamples (parents of the lifts)
  certificates/qstar/qstar_n7limit.json  q_7* (plastic number) for the comparison
  n7/labelg                              nauty canonical forms
Outputs:
  certificates/n8_summary.json                 the aggregate certificate
  certificates/n8_counterexamples_classes.json the D-equivalence classes with their (verified identical) interval data

Exact arithmetic only: fractions.Fraction, integer polynomials, Sturm sequences (fastregion.sturm_count),
exact bisection.  Decimals appear only in fields whose name says 'approx' and are display approximations.

D-equivalence is the Lemma J key EXACTLY as n7/aggregate_n7.py defines it (function d_key copied verbatim):
the integer characteristic polynomial of the adjacency matrix A (fastregion.icharpoly) together with the sorted
multiset over vertices of the in-walk count sequences ((A^T)^k 1)_i for k = 1..n (n = 8 here; Theorem J.2 of
notes/lemma_J_equivalence.md shows k <= n-1 already suffices, n7/aggregate_n7.py uses k = 1..n and so do we).
By Lemma J.0/J.2 two D-equivalent graphs have the same full-support data (d, phi, omega, multiset of a_i) and hence
the same limit-stable q-interval; this is verified on every class below.
"""
import sys, os, json, re, subprocess, collections, time
from fractions import Fraction as F
from multiprocessing import Pool

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE)
from n7.curve_test import from_d6                                   # digraph6 decoder (CLAUDE.md convention A[i][j] = 1 iff i -> j)
from fastregion import icharpoly, peval, pcompose_affine, primitive_int, sturm_count
from graphs import is_domination_free, is_complete_bidirected
from ctln import is_clique

LABELG = os.path.join(BASE, "n7", "labelg")
N = 8
WIDTH = F(1, 10 ** 9)          # target width of every refined isolating interval
Q_THRESHOLDS = [F(1), F(10001, 10000), F(1001, 1000), F(101, 100), F(21, 20), F(11, 10), F(6, 5), F(13, 10),
                F(7, 5), F(3, 2), F(8, 5), F(17, 10), F(7, 4), F(2)]


# ----------------------------------------------------------------------------------------------- graphs

def canon(d6s):
    """nauty canonical digraph6 strings (same labelg binary as n7/aggregate_n7.py)."""
    out = []
    for i in range(0, len(d6s), 2000):
        chunk = d6s[i:i + 2000]
        res = subprocess.run([LABELG], input="\n".join(chunk) + "\n", capture_output=True, text=True).stdout
        lines = [l for l in res.splitlines() if l and not l.startswith(">")]
        assert len(lines) == len(chunk), (len(lines), len(chunk))
        out += lines
    return out


def to_d6(A):
    """digraph6 encoder, inverse of n7.curve_test.from_d6 (row-major bits, 6 per character)."""
    n = len(A)
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6:
        bits.append(0)
    s = "&" + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        s += chr(v + 63)
    return s


def d_key(A):
    """VERBATIM from n7/aggregate_n7.py (Lemma J key)."""
    n = len(A)
    cp = tuple(icharpoly(A))
    v = [1] * n; walks = []
    for _ in range(n):
        v = [sum(A[j][i] * v[j] for j in range(n)) for i in range(n)]   # (A^T v)_i = sum_j A[j][i] v_j : in-walks ending at i
        walks.append(tuple(v))
    per_node = tuple(sorted(tuple(w[i] for w in walks) for i in range(n)))
    return (cp, per_node)


def bidirected_in_twins(A):
    """VERBATIM from n7/aggregate_n7.py."""
    n = len(A)
    return [(u, v) for u in range(n) for v in range(u + 1, n)
            if A[u][v] and A[v][u] and all(A[k][u] == A[k][v] for k in range(n) if k not in (u, v))]


def blow_up(A, v):
    """Closed-twin blow-up: new vertex n = v' with the in- and out-neighbourhood of v, plus v <-> v'
    (the construction of notes/verify_n8_2/completeness_twins/lib.py and of Lemma W)."""
    n = len(A)
    B = [row[:] + [0] for row in A] + [[0] * (n + 1)]
    for j in range(n):
        if j != v:
            B[n][j] = A[v][j]
            B[j][n] = A[j][v]
    B[v][n] = B[n][v] = 1
    return B


# ------------------------------------------------------------------------------- exact algebraic numbers

def pval(p_hi_to_lo, x):
    r = F(0)
    for c in p_hi_to_lo:
        r = r * x + c
    return r


class Alg:
    """A real algebraic number given by an integer polynomial (hi->lo) and an isolating interval.
    lo == hi means an exact rational.  All comparisons are exact."""
    __slots__ = ("p", "lo", "hi")

    def __init__(self, p_hi_to_lo, lo, hi):
        self.p = [int(c) for c in p_hi_to_lo]
        self.lo, self.hi = F(lo), F(hi)
        if self.lo != self.hi:
            if len(self.p) == 2:                       # degree 1: exact rational root
                self.lo = self.hi = F(-self.p[1], self.p[0])
            else:
                plo, phi = pval(self.p, self.lo), pval(self.p, self.hi)
                assert plo != 0 and phi != 0 and (plo > 0) != (phi > 0), "isolating interval must have a strict sign change"

    def exact(self):
        return self.lo == self.hi

    def refine(self, width):
        while self.hi - self.lo >= width:
            m = (self.lo + self.hi) / 2
            pm = pval(self.p, m)
            if pm == 0:                                 # rational root: only possible for degree 1 (irreducible), kept exact anyway
                self.lo = self.hi = m
                return
            if (pm > 0) == (pval(self.p, self.lo) > 0):
                self.lo = m
            else:
                self.hi = m

    def cmp_rational(self, r):
        """-1, 0, 1 for root <, =, > r.  One exact sign evaluation."""
        r = F(r)
        if self.exact():
            return (self.lo > r) - (self.lo < r)
        if r <= self.lo:
            return 1
        if r >= self.hi:
            return -1
        pr = pval(self.p, r)
        if pr == 0:
            self.lo = self.hi = r
            return 0
        if (pr > 0) == (pval(self.p, self.lo) > 0):
            self.lo = r                                 # root in (r, hi)
            return 1
        self.hi = r                                     # root in (lo, r)
        return -1

    def one_root_in(self, lo, hi):
        """Exactly one distinct real root of p in (lo, hi] (Sturm, exact)."""
        return sturm_count(list(reversed(self.p)), F(lo), F(hi)) == 1

    def cmp(self, other):
        """-1, 0, 1 for self <, =, > other (exact).  Equality of two interval-given numbers is decided by an
        identical polynomial plus a single root (Sturm) in the hull of the two intervals; different irreducible
        polynomials share no root, so refinement separates them."""
        if self.exact():
            return -other.cmp_rational(self.lo)
        if other.exact():
            return self.cmp_rational(other.lo)
        for _ in range(400):
            if self.hi <= other.lo:
                return -1
            if other.hi <= self.lo:
                return 1
            if self.p == other.p:
                lo, hi = min(self.lo, other.lo), max(self.hi, other.hi)
                if self.one_root_in(lo, hi):
                    return 0
            w = min(self.hi - self.lo, other.hi - other.lo) / 4
            self.refine(w); other.refine(w)
            if self.exact() or other.exact():
                return self.cmp(other)
        raise ArithmeticError("comparison did not terminate")

    def as_json(self):
        return {"minpoly_hi_to_lo": self.p, "isolating_interval": [str(self.lo), str(self.hi)],
                "degree": len(self.p) - 1, "approx_decimal_NOT_A_CERTIFICATE": float((self.lo + self.hi) / 2)}


def q_sturm_check(args):
    """Worker: verify the record's q end-point data exactly and refine it.  Returns plain data."""
    idx, hi, lo = args
    out = {"idx": idx}
    # (1) the q minimal polynomial stored in the record is exactly the beta minimal polynomial composed with beta = 1 + q
    bp = [F(v) for v in reversed(hi["beta_minpoly_hi_to_lo"])]
    qp = list(reversed(primitive_int(pcompose_affine(bp, F(1), F(1)))))
    out["q_minpoly_matches_beta_minpoly"] = (qp == hi["q_minpoly_hi_to_lo"])
    blo, bhi = F(hi["beta_interval"][0]), F(hi["beta_interval"][1])
    out["q_interval_matches_beta_interval"] = ([str(blo - 1), str(bhi - 1)] == hi["q_interval"])
    # (2) the stored isolating interval contains exactly one root (Sturm, exact)
    qlo, qhi = F(hi["q_interval"][0]), F(hi["q_interval"][1])
    if qlo == qhi:
        out["hi_one_root"] = (pval(hi["q_minpoly_hi_to_lo"], qlo) == 0)
    else:
        out["hi_one_root"] = (pval(hi["q_minpoly_hi_to_lo"], qlo) != 0 and
                              sturm_count(list(reversed(hi["q_minpoly_hi_to_lo"])), qlo, qhi) == 1)
    a = Alg(hi["q_minpoly_hi_to_lo"], qlo, qhi)
    a.refine(WIDTH)
    out["hi_refined"] = [str(a.lo), str(a.hi)]
    out["hi_cmp"] = [a.cmp_rational(thr) for thr in Q_THRESHOLDS]
    out["hi_lt_2001_over_1000_beta"] = a.cmp_rational(F(1001, 1000)) < 0     # beta < 2 + 1/1000  <=>  q < 1 + 1/1000
    # (3) lower end
    if lo["kind"] == "beta_min":
        out["lo"] = {"kind": "beta_min", "q": lo["q"], "beta": lo["beta"]}
        out["lo_one_root"] = True
        out["width_lt_1_over_1000"] = out["hi_lt_2001_over_1000_beta"]
        out["width_decided"] = True
    else:
        llo, lhi = F(lo["q_interval"][0]), F(lo["q_interval"][1])
        bp2 = [F(v) for v in reversed(lo["beta_minpoly_hi_to_lo"])]
        qp2 = list(reversed(primitive_int(pcompose_affine(bp2, F(1), F(1)))))
        out["q_minpoly_matches_beta_minpoly"] &= (qp2 == lo["q_minpoly_hi_to_lo"])
        if llo == lhi:
            out["lo_one_root"] = (pval(lo["q_minpoly_hi_to_lo"], llo) == 0)
        else:
            out["lo_one_root"] = (pval(lo["q_minpoly_hi_to_lo"], llo) != 0 and
                                  sturm_count(list(reversed(lo["q_minpoly_hi_to_lo"])), llo, lhi) == 1)
        b = Alg(lo["q_minpoly_hi_to_lo"], llo, lhi)
        b.refine(WIDTH)
        out["lo"] = {"kind": "root", "q_minpoly_hi_to_lo": lo["q_minpoly_hi_to_lo"], "refined": [str(b.lo), str(b.hi)],
                     "vanishing": lo["vanishing"], "degree": len(lo["q_minpoly_hi_to_lo"]) - 1}
        out["lo_lt_hi"] = (b.cmp(a) < 0)
        # width = hi - lo < 1/1000 ?  bounds: (a.lo - b.hi, a.hi - b.lo)
        w_lo, w_hi = a.lo - b.hi, a.hi - b.lo
        if w_hi < F(1, 1000):
            out["width_lt_1_over_1000"], out["width_decided"] = True, True
        elif w_lo > F(1, 1000):
            out["width_lt_1_over_1000"], out["width_decided"] = False, True
        else:
            out["width_lt_1_over_1000"], out["width_decided"] = None, False
    return out


def norm_vanishing(names):
    """Critical-polynomial names up to the vertex index of the positivity polynomials a_i = (adj(C) 1)_i:
    D-equivalence (Lemma J.0) preserves the multiset of the a_i, not their vertex labels."""
    return sorted("a*" if re.fullmatch(r"a\d+", nm) else nm for nm in names)


def irreducible(p):
    from sympy import Poly, symbols, ZZ
    x = symbols("x")
    return bool(Poly(p, x, domain=ZZ).is_irreducible)


def irreducible_worker(p):
    return (p, irreducible(list(p)))


# ------------------------------------------------------------------------------------------- main

def main():
    os.chdir(BASE)
    t0 = time.time()
    prov = {}

    # ---- enumeration chain
    tot = dict(zip(*[iter(open("n8/n8_prefilter_TOTAL.txt").read().split())] * 2))
    tot = {k: int(v) for k, v in tot.items()}
    acc = json.load(open("n8/n8_accounting.json"))
    rc = json.load(open("n8/n8_restricted_counts.json"))
    lim_sum = json.load(open("certificates/n8_limit_summary.json"))
    k8_ledger = open("n8/check_k8decomp.log").read().strip().splitlines()[-1]
    orbit_ledger = open("n8/check_orbit.log").read().strip().splitlines()
    assert tot["total"] == acc["aggregate"]["total"] == acc["space_total"] == rc["digraphs_total"] == 801558826511
    assert acc["ok"] and acc["n_problems"] == 0
    assert "OK" in k8_ledger and orbit_ledger[-1].strip() == "orbitenum ledger ok"
    prov["chain.enumeration"] = "n8/n8_prefilter_TOTAL.txt, n8/n8_accounting.json, n8/n8_restricted_counts.json"
    prov["chain.ledgers"] = "n8/check_k8decomp.log, n8/check_orbit.log"
    prov["chain.schemes"] = "certificates/n8_limit_summary.json (enumeration.schemes), ANOMALIES.md 2026-10-27 (two entries), notes/k8_decomposition.md"

    # ---- limit stage: tally recounted from the results file; sharded stage
    tally = json.load(open("n8/n8_limit_tally.json"))
    rec_count = collections.Counter(); n_results = 0; results_d6 = set()
    for line in open("n8/n8_limit_results.jsonl"):
        r = json.loads(line); n_results += 1; results_d6.add(r["d6"])
        if r["verdict"] == "never":
            rec_count["never_positivity" if r.get("reason") else "never_full"] += 1
        elif r["verdict"] == "stable_somewhere":
            rec_count["complete" if r["complete"] else "stable_noncomplete"] += 1
        else:
            rec_count[r["verdict"]] += 1
        rec_count["flag_" + r["flag"]] += 1
    assert n_results == len(results_d6) == tally["S"] + tally["O"] == 347639
    assert rec_count["never_positivity"] == tally["never_positivity"] and rec_count["never_full"] == tally["never_full"]
    assert rec_count["stable_noncomplete"] == tally["stable_noncomplete"] and rec_count["complete"] == tally["complete"]
    assert rec_count.get("escalate", 0) == tally["escalate"] == 0
    assert tot["survivors"] - tot["boundary"] == tally["S"] == rec_count["flag_S"]
    shards = [json.load(open(f"n8/limit_sharding/shard{k}_tally.json")) for k in range(4)]
    shard_sum = {k: sum(s[k] for s in shards) for k in ("S", "O", "never_positivity", "never_full", "stable_somewhere", "escalate", "complete", "stable_noncomplete")}
    assert all(shard_sum[k] == tally[k] for k in shard_sum)
    prov["counts.limit_stage"] = "n8/n8_limit_results.jsonl (recounted), n8/n8_limit_tally.json, n8/limit_sharding/shard{0..3}_tally.json"

    # ---- the hits
    hits = [json.loads(l) for l in open("n8/n8_limit_HITS.jsonl")]
    assert len(hits) == 9133 == rec_count["stable_noncomplete"]
    assert len({h["d6"] for h in hits}) == 9133 and all(h["d6"] in results_d6 for h in hits)
    mats = [from_d6(h["d6"]) for h in hits]
    assert all(len(A) == N for A in mats)
    assert all(to_d6(A) == h["d6"] for A, h in zip(mats, hits)), "digraph6 round trip"
    arcs = [sum(map(sum, A)) for A in mats]
    assert all(a == h["edges"] for a, h in zip(arcs, hits))
    assert not any(is_complete_bidirected(A) for A in mats) and not any(h["complete"] for h in hits)
    assert not any(is_clique(A, list(range(N))) for A in mats)
    domfree = sum(is_domination_free(A) for A in mats)
    canons = canon([h["d6"] for h in hits])
    assert len(set(canons)) == 9133, "hits pairwise non-isomorphic"
    assert all(h["verdict"] == "stable_somewhere" and h["flag"] == "S" and len(h["cert"]["stable_intervals"]) == 1 for h in hits)
    assert all(h["cert"]["eps"] == "0 (limit)" and h["cert"]["t"] is None and h["cert"]["beta_min"] == "2" for h in hits)
    assert all(h["cert"]["stable_intervals"][0]["hi"]["kind"] == "root" for h in hits)
    prov["hits"] = "n8/n8_limit_HITS.jsonl"

    # ---- D-equivalence classes (Lemma J key, verbatim aggregate_n7.d_key)
    keys = {}
    for i, A in enumerate(mats):
        keys.setdefault(d_key(A), []).append(i)
    classes = sorted(keys.values(), key=lambda c: (-len(c), [canons[i] for i in c]))
    class_of = {}
    for ci, c in enumerate(classes):
        for i in c:
            class_of[i] = ci
    class_sizes = sorted(collections.Counter(len(c) for c in classes).items())
    twins = [bidirected_in_twins(A) for A in mats]
    prov["D_equivalence"] = "computed here with d_key (verbatim n7/aggregate_n7.py) from n8/n8_limit_HITS.jsonl; definition notes/lemma_J_equivalence.md Theorem J.2"

    # ---- q end-points: exact checks and refinement (2 workers, nice)
    args = [(i, h["cert"]["stable_intervals"][0]["hi"], h["cert"]["stable_intervals"][0]["lo"]) for i, h in enumerate(hits)]
    with Pool(2) as pool:
        qres = pool.map(q_sturm_check, args, chunksize=50)
    qres.sort(key=lambda r: r["idx"])
    assert all(r["q_minpoly_matches_beta_minpoly"] for r in qres)
    assert all(r["q_interval_matches_beta_interval"] for r in qres)
    assert all(r["hi_one_root"] and r["lo_one_root"] for r in qres)
    assert all(r.get("lo_lt_hi", True) for r in qres)
    lo_kinds = collections.Counter(r["lo"]["kind"] for r in qres)
    hi_polys = {}
    for r, h in zip(qres, hits):
        p = tuple(h["cert"]["stable_intervals"][0]["hi"]["q_minpoly_hi_to_lo"])
        hi_polys.setdefault(p, []).append(r["idx"])
    lo_polys = {tuple(r["lo"]["q_minpoly_hi_to_lo"]) for r in qres if r["lo"]["kind"] == "root"}
    with Pool(2) as pool:
        irr = dict(pool.map(irreducible_worker, sorted(set(hi_polys) | lo_polys, key=len), chunksize=20))
    assert all(irr.values()), "every stored minimal polynomial must be irreducible"
    prov["q_data"] = "n8/n8_limit_HITS.jsonl (cert.stable_intervals[0].hi / .lo), refined here by exact bisection; cross-checked against notes/verify_n8/v_structure/q_upper_ends_refined.json"

    # ---- q_8*(limit): exact supremum over the 9,133 upper ends
    algs = [Alg(h["cert"]["stable_intervals"][0]["hi"]["q_minpoly_hi_to_lo"], *r["hi_refined"]) for h, r in zip(hits, qres)]
    best = 0
    for i in range(1, len(algs)):
        if algs[i].cmp(algs[best]) > 0:
            best = i
    ties = [i for i in range(len(algs)) if i != best and algs[i].cmp(algs[best]) == 0]
    others_below = all(algs[i].cmp(algs[best]) < 0 for i in range(len(algs)) if i != best and i not in ties)
    assert others_below
    algs[best].refine(WIDTH)
    qstar = algs[best]
    # certified ordering: after refinement, every non-tied upper end's interval lies strictly below the champion's
    ordering_ok = all(algs[i].hi <= qstar.lo for i in range(len(algs)) if i != best and i not in ties)
    assert ordering_ok
    hi_best = hits[best]["cert"]["stable_intervals"][0]["hi"]
    # n = 7 comparison
    q7 = json.load(open("certificates/qstar/qstar_n7limit.json"))["n7_limit_CANDIDATES_lower_bound"]["max_beta_c_minus_1"]
    q7a = Alg(q7["minpoly_hi_to_lo"], *q7["isolating_interval"])
    q8_gt_q7 = qstar.cmp(q7a) > 0
    prov["q_star.n7_comparison"] = "certificates/qstar/qstar_n7limit.json (plastic number), ANOMALIES.md 2026-10-15 'q_7* is the plastic number'"

    # ---- cross-check with the first campaign's refinement (v_structure)
    ref = json.load(open("notes/verify_n8/v_structure/q_upper_ends_refined.json"))
    ref_all = {r["d6"]: r for r in ref["all"]}
    assert set(ref_all) == {h["d6"] for h in hits}
    overlap_ok = 0
    for h, a in zip(hits, algs):
        rlo, rhi = [F(v) for v in ref_all[h["d6"]]["q_hi_refined_interval"]]
        if max(rlo, a.lo) <= min(rhi, a.hi):
            overlap_ok += 1
    assert overlap_ok == 9133, "every independently refined interval must overlap ours"
    assert ref["largest"]["d6"] in {hits[best]["d6"]} | {hits[i]["d6"] for i in ties}

    # ---- distribution of the upper ends
    above = {str(thr): sum(1 for r in qres if r["hi_cmp"][k] > 0) for k, thr in enumerate(Q_THRESHOLDS)}
    equal = {str(thr): sum(1 for r in qres if r["hi_cmp"][k] == 0) for k, thr in enumerate(Q_THRESHOLDS)}
    bins = []
    for k in range(len(Q_THRESHOLDS) - 1):
        lo_t, hi_t = Q_THRESHOLDS[k], Q_THRESHOLDS[k + 1]
        cnt = sum(1 for r in qres if r["hi_cmp"][k] > 0 and r["hi_cmp"][k + 1] <= 0)
        bins.append({"q_in": f"({lo_t}, {hi_t}]", "count": cnt})
    assert sum(b["count"] for b in bins) + above[str(Q_THRESHOLDS[-1])] + equal[str(Q_THRESHOLDS[0])] + sum(1 for r in qres if r["hi_cmp"][0] < 0) == 9133
    deg_hist = sorted(collections.Counter(len(h["cert"]["stable_intervals"][0]["hi"]["q_minpoly_hi_to_lo"]) - 1 for h in hits).items())
    van_hist = sorted(collections.Counter(",".join(h["cert"]["stable_intervals"][0]["hi"]["vanishing"]) for h in hits).items())
    lo_van_hist = sorted(collections.Counter(",".join(h["cert"]["stable_intervals"][0]["lo"]["vanishing"]) for h in hits if h["cert"]["stable_intervals"][0]["lo"]["kind"] == "root").items())
    minpoly_counts = sorted(((len(v), list(p)) for p, v in hi_polys.items()), key=lambda x: (-x[0], x[1]))
    narrow_from_2 = sum(1 for r in qres if r["lo"]["kind"] == "beta_min" and r["width_lt_1_over_1000"])
    narrow_root = sum(1 for r in qres if r["lo"]["kind"] == "root" and r["width_lt_1_over_1000"] is True)
    narrow_undecided = sum(1 for r in qres if not r["width_decided"])
    # smallest upper end
    worst = 0
    for i in range(1, len(algs)):
        if algs[i].cmp(algs[worst]) < 0:
            worst = i
    algs[worst].refine(WIDTH)

    # ---- lifts: closed-twin blow-ups of the 123 n = 7 limit counterexamples (recomputed) vs the second campaign's file
    n7hits = [json.loads(l) for l in open("n7/n7_limit_HITS.jsonl")]
    assert len(n7hits) == 123
    bl_d6, bl_parent = [], []
    for pi, r in enumerate(n7hits):
        P = from_d6(r["d6"]); assert len(P) == 7
        for v in range(7):
            bl_d6.append(to_d6(blow_up(P, v))); bl_parent.append((pi, v))
    bl_canon = canon(bl_d6)
    assert len(bl_d6) == 861
    distinct_blowups = sorted(set(bl_canon))
    hit_by_canon = {c: i for i, c in enumerate(canons)}
    lifted = {}                                                     # hit index -> list of (parent d6, v)
    for c, (pi, v) in zip(bl_canon, bl_parent):
        if c in hit_by_canon:
            lifted.setdefault(hit_by_canon[c], []).append([n7hits[pi]["d6"], v])
    lift_parents = sorted({p for L in lifted.values() for p, _ in L})
    # the second campaign's file
    p2 = [json.loads(l) for l in open("notes/verify_n8_2/completeness_twins/part2_blowups.jsonl")]
    p2_in = {r["labelg"] for r in p2 if r["in_hits_labelg"]}
    p2_parents = {r["parent_d6"] for r in p2 if r["in_hits_labelg"]}
    p2_distinct = {r["labelg"] for r in p2}
    assert p2_in == {canons[i] for i in lifted}, "lifted set must equal the second campaign's"
    assert len(p2_distinct) == len(distinct_blowups) == 858 and len(p2_in) == 61
    assert set(p2_parents) == set(lift_parents)
    p34 = json.load(open("notes/verify_n8_2/completeness_twins/part34_summary.json"))
    prov["lifts"] = "recomputed here (blow_up + n7/labelg) from n7/n7_limit_HITS.jsonl; cross-checked against notes/verify_n8_2/completeness_twins/part2_blowups.jsonl"
    prov["with_bidirected_in_twins"] = "computed here (bidirected_in_twins, verbatim n7/aggregate_n7.py); cross-checked against notes/verify_n8_2/completeness_twins/part34_summary.json part3.records_with_twins"
    with_twins = sum(1 for t in twins if t)
    assert with_twins == p34["part3"]["records_with_twins"] == 4271
    assert sum(len(t) for t in twins) == p34["part3"]["twin_pairs"]

    # ---- class consistency: identical interval data within every D-class (Lemma J.0)
    exceptions = []
    class_rows = []
    label_only_diffs = 0        # members whose 'vanishing' names differ only in the vertex index of a_i (a_i are vertex-indexed; D-equivalence preserves their multiset only)
    for ci, c in enumerate(classes):
        i0 = c[0]; r0 = qres[i0]; h0 = hits[i0]["cert"]["stable_intervals"][0]
        rep_hi = Alg(h0["hi"]["q_minpoly_hi_to_lo"], *r0["hi_refined"])
        rep_lo = None if r0["lo"]["kind"] == "beta_min" else Alg(r0["lo"]["q_minpoly_hi_to_lo"], *r0["lo"]["refined"])
        same = True
        for i in c[1:]:
            r = qres[i]; h = hits[i]["cert"]["stable_intervals"][0]
            prob = []
            if arcs[i] != arcs[i0]:
                prob.append("arcs differ")
            if r["lo"]["kind"] != r0["lo"]["kind"]:
                prob.append("lo kind differs")
            elif r["lo"]["kind"] == "root":
                if r["lo"]["q_minpoly_hi_to_lo"] != r0["lo"]["q_minpoly_hi_to_lo"]:
                    prob.append("lo minpoly differs")
                elif Alg(r["lo"]["q_minpoly_hi_to_lo"], *r["lo"]["refined"]).cmp(rep_lo) != 0:
                    prob.append("lo root differs")
            if h["hi"]["q_minpoly_hi_to_lo"] != h0["hi"]["q_minpoly_hi_to_lo"]:
                prob.append("hi minpoly differs")
            elif Alg(h["hi"]["q_minpoly_hi_to_lo"], *r["hi_refined"]).cmp(rep_hi) != 0:
                prob.append("hi root differs")
            if norm_vanishing(h["hi"]["vanishing"]) != norm_vanishing(h0["hi"]["vanishing"]):
                prob.append("hi vanishing set differs (up to vertex index of a_i)")
            if h["hi"]["vanishing"] != h0["hi"]["vanishing"]:
                label_only_diffs += 1
            if prob:
                same = False
                exceptions.append({"class": ci, "member": hits[i]["d6"], "representative": hits[i0]["d6"], "problems": prob})
        rep_hi.refine(WIDTH)
        if rep_lo is not None:
            rep_lo.refine(WIDTH)
        class_rows.append({
            "class": ci, "size": len(c), "members_d6": [hits[i]["d6"] for i in c], "members_canonical": [canons[i] for i in c],
            "arcs": arcs[i0], "missing_arcs": N * (N - 1) - arcs[i0],
            "members_with_bidirected_in_twins": sum(1 for i in c if twins[i]),
            "members_lifted_from_n7": sum(1 for i in c if i in lifted),
            "interval_data_identical_within_class": same,
            "stable_q_interval": {
                "lo": {"kind": "beta_min", "q": "1", "beta": "2", "note": "open end at the legality limit beta -> 2+, q -> 1+"} if rep_lo is None else
                      {"kind": "root", **rep_lo.as_json(), "vanishing": r0["lo"]["vanishing"]},
                "hi": {"kind": "root", **rep_hi.as_json(), "vanishing": h0["hi"]["vanishing"],
                       "beta_minpoly_hi_to_lo": h0["hi"]["beta_minpoly_hi_to_lo"], "beta_interval_stored": h0["hi"]["beta_interval"]},
                "open": True,
            },
        })

    # ---- per-record rows (compact: minimal polynomials by table index)
    poly_table = [list(p) for _, p in minpoly_counts]
    poly_index = {tuple(p): k for k, p in enumerate(poly_table)}
    rows = []
    for i, (h, r) in enumerate(zip(hits, qres)):
        hi = h["cert"]["stable_intervals"][0]["hi"]
        rows.append({
            "d6": h["d6"], "canonical": canons[i], "arcs": arcs[i], "class": class_of[i],
            "bidirected_in_twin_pairs": twins[i],
            "lifted_from_n7": lifted.get(i),
            "lo": r["lo"],
            "hi_q_minpoly_index": poly_index[tuple(hi["q_minpoly_hi_to_lo"])],
            "hi_q_interval_stored": hi["q_interval"], "hi_q_interval_refined": r["hi_refined"],
            "hi_q_interval_v_structure": ref_all[h["d6"]]["q_hi_refined_interval"],
            "hi_vanishing": hi["vanishing"], "hi_degree": len(hi["q_minpoly_hi_to_lo"]) - 1,
            "beta_interval_narrower_than_1_over_1000": r["width_lt_1_over_1000"],
            "witnesses_all_stable_positive": all(w["stable"] and w["positive"] for w in h["witnesses"]),
            "has_legal_stable_positive_stored_witness": any(w["stable"] and w["positive"] and F(w["q"]) > 1 for w in h["witnesses"]),
        })

    syn1 = json.load(open("notes/verify_n8/SYNTHESIS.json"))
    syn2 = json.load(open("notes/verify_n8_2/SYNTHESIS.json"))

    summary = {
        "n": N,
        "date": "2026-10-28",
        "scope": ("eps -> 0 LIMIT map of the FULL support of every digraph on 8 vertices with >= 29 arcs (equivalently, after the e_2 edge "
                  "bound of Lemma A.2, notes/lean/CtlnCert/CtlnCert/EdgeBound.lean: underlying graph with 15..28 edges and >= 29 arcs; "
                  "every sparser digraph has no stable full support at any legal parameter). NOT done: the finite-eps map at n = 8; the "
                  "417,603,634 boundary (B) lines of the prefilter (alpha(A|perp) >= -1/2, never limit-stable by Lemma I.3(c2)) were not "
                  "stored, so finite-eps-only counterexamples among them (alpha = -1/2 exactly, the FNFL type of n = 7) are outside this run. "
                  "Proper supports are not covered (Theorem F reduces them to induced subgraphs on <= 7 nodes, mapped in notes/REPORT_n7.md)."),
        "status": lim_sum["status"],
        "chain": {
            "restricted_space": {"underlying_edges": rc["geng_edges"], "min_arcs": rc["min_arcs"], "underlying_graphs": rc["underlying_graphs"],
                                 "digraphs_exact_burnside": rc["digraphs_total"], "method": rc["method"],
                                 "A000273_8_for_reference": 1793359192848},
            "digraphs_read_exactly_once": tot["total"],
            "accounting": {"ok": acc["ok"], "problems": acc["problems"], "shards_checked": acc["shards_checked"], "jobs_checked": acc["jobs_checked"],
                           "split_graphs_checked": acc["split_graphs_checked"], "kept_lines": acc["kept_lines"]},
            "enumeration_schemes": lim_sum["enumeration"]["schemes"],
            "swaps": [
                {"what": "K8 (575,009,339 digraphs, |Aut| = 40,320): 348 directg parts replaced by 40 batches enumerated by bidirected subgraph B (n8/k8decomp)",
                 "lemma": "notes/k8_decomposition.md (isomorphism classes of digraphs with underlying graph G <-> pairs (B up to Aut(G), orientation of G - B up to the stabiliser of B))",
                 "ledger": k8_ledger, "anomaly_entry": "ANOMALIES.md 2026-10-27 'K8 enumerated by its bidirected subgraph'"},
                {"what": "K8-K3, K8-K_{1,5}, K8-(K2+K_{1,5}): 99 directg parts replaced by 14 orbitenum parts (n8/orbitenum)",
                 "ledger": orbit_ledger, "anomaly_entry": "ANOMALIES.md 2026-10-27 'Three more high-symmetry graphs enumerated by orbits'",
                 "post_hoc_content_check": "K8-e re-enumerated by orbitenum: 143,266 canonical strict-survivor lines identical to the 321 directg parts (ANOMALIES.md 2026-10-28, n8/orbitenum/k8e_check)"},
            ],
            "prefilter": {"domination_free_LemmaD": tot["domfree"], "spectral_survivors": tot["survivors"],
                          "boundary_B_alpha_ge_minus_half_never_limit_stable": tot["boundary"], "strict_S": tot["survivors"] - tot["boundary"],
                          "overflow_O": tot["overflow"], "B_lines_stored": False},
            "limit_stage": {"engine": "fastregion.region(A, range(8), INF, replay=True) on every S/O line, n8/run_n8_limit.py, exact (no floating point)",
                            "sharded": "interrupted single run at 59,610 lines; finished as 4 shards of the unmodified script over disjoint job subsets, merged by n8/limit_sharding/merge_finalize.py",
                            "shard_tallies": shards, "shard_sum_equals_merged_tally": True,
                            "interrupted_single_run_records_cross_checked": tally["interrupted_single_run_records_cross_checked"],
                            "seconds": tally["seconds"]},
        },
        "counts": {
            "digraphs_read": tot["total"], "domination_free": tot["domfree"], "spectral_survivors": tot["survivors"],
            "boundary": tot["boundary"], "strict": tot["survivors"] - tot["boundary"], "overflow": tot["overflow"],
            "limit_records": n_results, "never_positivity": rec_count["never_positivity"], "never_full": rec_count["never_full"],
            "stable_somewhere": rec_count["stable_noncomplete"] + rec_count["complete"], "escalate": rec_count.get("escalate", 0),
            "complete": rec_count["complete"], "candidates": rec_count["stable_noncomplete"],
            "accounting_identity": "strict == never_positivity + never_full + stable_somewhere + escalate",
            "accounting_identity_holds": tot["survivors"] - tot["boundary"] == rec_count["never_positivity"] + rec_count["never_full"] + rec_count["stable_noncomplete"] + rec_count["complete"] + rec_count.get("escalate", 0),
        },
        "hits": len(hits),
        "hits_pairwise_nonisomorphic_labelg": True,
        "hits_all_noncomplete_nonclique": True,
        "hits_domination_free": domfree,
        "D_equivalence_definition": ("Lemma J key exactly as n7/aggregate_n7.py d_key: (integer characteristic polynomial of A, sorted multiset over "
                                     "vertices of the in-walk sequences ((A^T)^k 1)_i, k = 1..8).  Theorem J.2 (notes/lemma_J_equivalence.md): "
                                     "equal keys <=> same full-support data D(A) = (d, phi, omega, multiset a_i) <=> identical full-support region "
                                     "for every (eps, delta)."),
        "D_equivalence_classes": len(classes),
        "class_sizes": class_sizes,
        "classes_singletons": sum(1 for c in classes if len(c) == 1),
        "largest_class": len(classes[0]),
        "class_interval_data_identical": not exceptions,
        "class_exceptions": exceptions,
        "class_members_whose_vanishing_labels_differ_only_by_vertex_index": label_only_diffs,
        "edge_counts": sorted(collections.Counter(arcs).items()),
        "missing_arc_counts": sorted(collections.Counter(N * (N - 1) - a for a in arcs).items()),
        "min_arcs_hits": [hits[i]["d6"] for i in range(len(hits)) if arcs[i] == min(arcs)],
        "with_bidirected_in_twins": with_twins,
        "bidirected_in_twin_pairs_total": sum(len(t) for t in twins),
        "lemma_J_switch_closure_second_campaign": p34["part3"],
        "lifts": {
            "definition": "closed-twin blow-up of an n = 7 limit counterexample at a vertex v: add v' with the in- and out-neighbourhood of v and v <-> v'",
            "parents": 123, "blowups_built": 861, "distinct_blowups_labelg": len(distinct_blowups),
            "distinct_blowups_in_hits": len(lifted), "parents_with_a_lift_in_hits": len(lift_parents),
            "hits_not_a_single_closed_twin_lift": len(hits) - len(lifted),
            "blowups_not_in_hits": len(distinct_blowups) - len(lifted),
            "blowups_not_in_hits_status": "all exactly decided NOT limit-stable at legal points inside the parent's interval (notes/verify_n8_2/completeness_twins/part2_summary.json: 797, 0 missed)",
            "agrees_with_second_campaign": True,
            "lifted_hits": [{"d6": hits[i]["d6"], "canonical": canons[i], "parents": lifted[i]} for i in sorted(lifted)],
            "lift_parents_d6": lift_parents,
        },
        "q_data": {
            "coordinates": "limit t -> infinity (eps -> 0): q = beta - 1, legal iff beta > 2 i.e. q > 1; every candidate has exactly ONE stable interval (open cell of beta between consecutive critical roots)",
            "what_an_endpoint_of_kind_root_certifies": ("the end-point is an isolated real root (Sturm count 1 on the stored interval, re-verified here) of the "
                                                        "irreducible integer polynomial stored as minpoly (irreducibility re-verified here for all distinct polynomials); "
                                                        "the 'vanishing' field names the critical polynomial(s) that vanish there (Delta = the Hurwitz determinant of the "
                                                        "characteristic polynomial of -I+W, a_k = coefficient k of that characteristic polynomial, D = det K (positivity), "
                                                        "o_k = off-support drive). The record certifies stability on the OPEN cell by an exact sample (constant verdict on a "
                                                        "cell with no critical root inside, replay-checked isolation); AT the end-point a Hurwitz minor vanishes, which is "
                                                        "three-valued undecided, so the end-point itself is never certified stable."),
            "lower_end_kinds": dict(lo_kinds),
            "intervals_starting_at_beta_2": lo_kinds["beta_min"],
            "intervals_starting_above_beta_2": lo_kinds["root"],
            "intervals_narrower_than_1_over_1000_in_beta": {"starting_at_beta_2": narrow_from_2, "starting_above_2": narrow_root, "undecided_at_width_1e-9": narrow_undecided},
            "first_campaign_count_for_comparison": lim_sum["machine_confirmation"]["findings"]["graphs_whose_stored_witnesses_all_lie_outside_a_stable_interval_narrower_than_1_over_1000"],
            "stored_minpolys_checked": {"q_minpoly_equals_beta_minpoly_at_beta_1_plus_q": True, "q_interval_equals_beta_interval_minus_1": True,
                                        "sturm_count_1_on_every_stored_isolating_interval": True, "all_distinct_minpolys_irreducible": True,
                                        "distinct_hi_minpolys": len(hi_polys), "distinct_lo_minpolys": len(lo_polys)},
            "hi_degree_histogram": deg_hist,
            "hi_vanishing_histogram": van_hist,
            "lo_vanishing_histogram_root_cases": lo_van_hist,
            "upper_ends_above": above, "upper_ends_equal_to": equal,
            "upper_end_histogram": bins,
            "upper_ends_above_3_over_2": above["3/2"], "upper_ends_above_11_over_10": above["11/10"], "upper_ends_above_6_over_5": above["6/5"],
            "all_upper_ends_below_2": above["2"] == 0 and equal["2"] == 0,
            "all_upper_ends_above_1": above["1"] == 9133,
            "smallest_upper_end": {"d6": hits[worst]["d6"], **algs[worst].as_json(), "vanishing": hits[worst]["cert"]["stable_intervals"][0]["hi"]["vanishing"]},
            "smallest_upper_end_v_structure_interval": ref_all[hits[worst]["d6"]]["q_hi_refined_interval"],
            "smallest_upper_end_discrepancy_note": ("the first campaign's synthesis line (notes/verify_n8/SYNTHESIS.json: 'smallest ~1.000276', copied into "
                                                    "certificates/n8_limit_summary.json machine_confirmation.findings.q_upper_end_min_approx) disagrees with that "
                                                    "campaign's own per-graph file notes/verify_n8/v_structure/q_upper_ends_refined.json, whose smallest entry is "
                                                    "&GMBQMHtMH~]W in [77079731/77070336, 1079116235/1078984704] (about 1.000122); the exact refinement here agrees "
                                                    "with the per-graph file. A display value only; no certificate depends on it. To be logged in ANOMALIES.md by the orchestrator."),
            "most_common_hi_minpolys": [{"count": c, "degree": len(p) - 1, "minpoly_hi_to_lo": p} for c, p in minpoly_counts[:12]],
            "q_minpoly_table": poly_table,
        },
        "q_star_limit": {
            "definition": "q_8*(limit) := sup over the 9,133 candidates of the upper end of the stable q-interval (limit map; the finite-eps q_8* is NOT determined by this run)",
            **qstar.as_json(),
            "isolating_interval_width_lt": str(WIDTH),
            "beta_minpoly_hi_to_lo": hi_best["beta_minpoly_hi_to_lo"],
            "beta_interval_stored": hi_best["beta_interval"],
            "vanishing_at_q_star": hi_best["vanishing"],
            "attained_by": sorted([hits[best]["d6"]] + [hits[i]["d6"] for i in ties]),
            "attained_by_canonical": sorted([canons[best]] + [canons[i] for i in ties]),
            "attained_by_same_D_class": len({class_of[best]} | {class_of[i] for i in ties}) == 1,
            "arcs_of_attaining_graphs": sorted({arcs[best]} | {arcs[i] for i in ties}),
            "ordering_certified": ordering_ok,
            "ordering_certificate": "every other upper end has an isolating interval (width < 1e-9, exact bisection) lying entirely below the champion's; the tie is the same root (identical irreducible polynomial, Sturm count 1 on the hull of the two intervals)",
            "attained": False,
            "attained_note": ("NOT attained in the limit map: the stable set of each candidate is the open cell below this root (the record's sample point and "
                              "replay-checked isolation certify the open cell only); at q = q_8* the Hurwitz determinant Delta of the attaining graphs vanishes "
                              "(three-valued undecided, not stable), and every other candidate's upper end is certified strictly smaller, so no candidate is "
                              "certified stable at or above q_8*."),
            "compare_with_n7": {"q_7_star_limit": q7, "q_8_star_gt_q_7_star": q8_gt_q7, "q_7_star_note": "plastic number, not attained (ANOMALIES.md 2026-10-15)"},
            "agrees_with_first_campaign_refinement": ref["largest"]["d6"] in {hits[best]["d6"]} | {hits[i]["d6"] for i in ties},
        },
        "witness_summary": {
            "hits_with_all_four_pointwise_witnesses_stable_and_positive": sum(1 for r in rows if r["witnesses_all_stable_positive"]),
            "hits_with_a_legal_stable_positive_stored_witness_q_gt_1": sum(1 for r in rows if r["has_legal_stable_positive_stored_witness"]),
            "legal_witness_certificate": lim_sum["machine_confirmation"]["legal_witness_certificate"],
        },
        "verification": {
            "status": "VERIFIED BY MACHINE (two independent campaigns); NO human hand check at n = 8",
            "first_campaign": {"workflow": lim_sum["machine_confirmation"]["workflow"], "date": lim_sum["machine_confirmation"]["date"],
                               "synthesis": "notes/verify_n8/SYNTHESIS.json", "overall_status": syn1["overall_status"], "status_line": syn1["status_line"],
                               "verifiers": list(syn1["verifiers"].keys()), "per_verifier": lim_sum["machine_confirmation"]["per_verifier"],
                               "findings": lim_sum["machine_confirmation"]["findings"]},
            "second_campaign": {"workflow": lim_sum["second_verification"]["workflow"], "date": lim_sum["second_verification"]["date"],
                                "synthesis": "notes/verify_n8_2/SYNTHESIS.json", "status": lim_sum["second_verification"]["status"],
                                "lenses": lim_sum["second_verification"]["lenses"], "defects": lim_sum["second_verification"]["defects"],
                                "agent_directories": sorted(d for d in os.listdir("notes/verify_n8_2") if os.path.isdir(os.path.join("notes/verify_n8_2", d))),
                                "hand_check_kit": "notes/verify_n8_2/handcheck_kit/HANDCHECK.md"},
            "human_decision": lim_sum["second_verification"]["human_decision"],
            "anomaly_entries": ["ANOMALIES.md 2026-10-27 (seventh restart; K8 by bidirected subgraph; three graphs by orbits)",
                                "ANOMALIES.md 2026-10-28 (limit stage complete; K8-e orbit check; first campaign; hand check declined; chain referee corrections; second campaign clean)"],
        },
        "n7_comparison": {"n7_limit_counterexamples": 123, "n7_finite_eps_only": 1, "n7_finite_eps_only_d6": "&FNFL~[udo?",
                          "n8_limit_counterexamples": len(hits), "n8_lifted_from_n7": len(lifted), "n8_new": len(hits) - len(lifted),
                          "source": "certificates/n7_finite_summary.json, certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit_summary.json"},
        "rows": rows,
        "files": {"hits": "n8/n8_limit_HITS.jsonl", "results": "n8/n8_limit_results.jsonl", "tally": "n8/n8_limit_tally.json",
                  "limit_summary": "certificates/n8_limit_summary.json", "legal_witnesses": "certificates/n8_candidates_legal_witnesses.json",
                  "classes": "certificates/n8_counterexamples_classes.json", "this": "certificates/n8_summary.json", "writer": "n8/aggregate_n8.py"},
        "provenance": prov,
        "seconds": round(time.time() - t0, 1),
    }
    prov["counts.chain"] = "n8/n8_prefilter_TOTAL.txt (total, domfree, survivors, boundary, overflow); strict = survivors - boundary"
    prov["status"] = "certificates/n8_limit_summary.json"
    prov["verification"] = "certificates/n8_limit_summary.json (machine_confirmation, second_verification), notes/verify_n8/SYNTHESIS.json, notes/verify_n8_2/SYNTHESIS.json"
    prov["witness_summary"] = "n8/n8_limit_HITS.jsonl (witnesses), certificates/n8_limit_summary.json"
    prov["edge_counts, missing_arc_counts"] = "computed here from the decoded adjacency matrices of n8/n8_limit_HITS.jsonl (equal to the record field 'edges')"
    prov["class_sizes, classes"] = "computed here; written to certificates/n8_counterexamples_classes.json"
    prov["q_data.first_campaign_count_for_comparison"] = "certificates/n8_limit_summary.json machine_confirmation.findings"

    json.dump(summary, open("certificates/n8_summary.json", "w"), indent=1)
    json.dump({"n": N, "source": "n8/n8_limit_HITS.jsonl", "writer": "n8/aggregate_n8.py",
               "definition": summary["D_equivalence_definition"], "classes": len(classes), "class_sizes": class_sizes,
               "all_classes_have_identical_interval_data": not exceptions, "exceptions": exceptions,
               "members_whose_vanishing_labels_differ_only_by_vertex_index": label_only_diffs,
               "note_on_vanishing_labels": "a_i names are vertex-indexed positivity polynomials; D-equivalence preserves their multiset only, so a2 vs a3 at the same root is the same datum",
               "rows": class_rows}, open("certificates/n8_counterexamples_classes.json", "w"), indent=1)
    key = {k: v for k, v in summary.items() if k not in ("rows", "q_data", "lifts", "chain", "verification", "class_exceptions")}
    key["q_data_short"] = {k: v for k, v in summary["q_data"].items() if k not in ("q_minpoly_table", "most_common_hi_minpolys", "what_an_endpoint_of_kind_root_certifies")}
    key["lifts_short"] = {k: v for k, v in summary["lifts"].items() if k not in ("lifted_hits", "lift_parents_d6")}
    print(json.dumps(key, indent=1, default=str))


if __name__ == "__main__":
    main()
