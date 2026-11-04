# Geneson report 2026-11-04, very minor point 1: off-support curves missing from the Section 3.1 list

Scope: the paragraph of `notes/paper/paper.tex` beginning "The two-parameter version" (lines 422-425) lists the
curves a_i = 0, D = 0, Delta_{m-1} = 0, L = 0 and says the SFP verdict is constant on connected components of
their complement. For a proper support the off-support curves o_k = 0 (k not in sigma) are missing. Lemma 3.2
(lem:B, line 409), the definition of Z (line 437), and Lemma B' / Z_0 in `notes/lemma_Bprime_2d.md` (line 70)
all include them. The paper was not edited here.

Inputs: the uploaded report is byte-identical to `notes/referee_external/geneson_report_2026-11-04.pdf`
(sha256 e5d427a7...0beca).

## (a) The example, reproduced exactly

Script: `notes/referee_external/fixes_2026-11-04/prover/offsupport_check.py`
(output `offsupport_check.json`, stdout `offsupport_check_stdout.txt`; about 6 s).
It reimplements N, C, K, a, d, D, p, o_k and L from Section 2 and checks each identity it uses: the closed form
(2.2) with A_sigma^T, I - W_sigma = eps K, (I - W_sigma) x = theta 1 for x = (theta/eps) a/D, and the outside
drive of Definition 2.2(2) equal to theta o_7 / D. That last identity is derived symbolically from W, so it is not
assumed from Lemma 2.1(v).

Graph: G_1 in the **Table 1 labelling** (rows 011000 000011 000111 110011 110101 101110). Vertex 7 (index 6)
receives arcs from 2, 3, 5, 6 (indices 1, 2, 4, 5) and sends none (A[i][6] = 1 for those i, row 6 is zero).
Then N_{7,b} = 1 exactly for b in {1, 4}, so o_7 = d - beta (a_1 + a_4).

- o_7 = -beta^5 + 3 beta^4 - 3 beta^3 + 3 beta^2 - 2 beta + 1, which is identical to the report's polynomial
  (ratio 1, no extra normalisation). It is irreducible over Q.
- o_7(52/25) = 417793/9765625 > 0 and o_7(209/100) = -1006161749/10^10 < 0, as the report states. o_7 has exactly
  one root in [52/25, 209/100], at beta* ~ 2.0830248 (printed only). An exact isolating interval is
  (13161411/6318413, 13052825/6266284), and o_7(25/12) = -1093/248832 < 0.
- At t = 100, the following have 0 roots in the closed segment [52/25, 209/100], by sympy's exact count and by an
  independent Sturm chain over Fraction, and are positive at both endpoints:
  a_1..a_6 (degrees 4,3,4,4,2,2), D (degree 5), Delta_1..Delta_6 (degrees 0,3,6,10,15,20) and L = 100 beta - 201.
- Verdict by Lemma 2.3 at (beta, t) = (52/25, 100): positivity yes, Hurwitz yes, legal yes, o_7 < 0 **no**, so
  sigma is not an SFP. At (209/100, 100) all four conditions hold, so sigma **is** an SFP.
- Cross-check with `ctln.analyze_support` (project W convention) on the 7-node digraph. At eps = 1/101,
  delta = 27/2525 (q = 108/101), the result is positive, stable, off_support_ok = False, outside drive
  417793/18223613457 > 0, so it is not an SFP. At eps = 1/101, delta = 109/10100 (q = 109/101) the result is
  positive, stable, off_support_ok = True, outside drive -1006161749/19005064281749 < 0, so it is an SFP. Both
  drives equal o_7/D exactly. Neither point is BOUNDARY.

So the report is confirmed: the segment lies in one component of the complement of the listed curves, yet the
verdict changes across o_7 = 0. The paragraph's invariance claim is false for proper supports as written.

The labelling matters. In the stage-0 labelling of Geneson's graph (`certificates/stage0_geneson_2607.21396.json`),
the same 1-based in-set {2,3,5,6} gives a different o_7 = -b^5+3b^4-3b^3+4b^2-2b+1. That polynomial has no root in
the segment and is positive there, so sigma is an SFP at neither endpoint. The only isomorphism from Table-1 G_1 to
the stage-0 labelling is phi = (0,1,2,3,4,5) -> (5,4,3,0,1,2). It maps the report's in-set to stage-0 indices
{1,2,3,4}.

This stable non-clique proper support is a counterexample instance, but it is not new. It is record
`base = geneson_idx1000255, pattern = 30` of the certified n = 7 proper-support list (below), and its base G_1 is
hand-verified. No STOP is needed.

## (b) Implementation: the off-support polynomials are in both engines

- `fastregion.beta_data` (lines ~265-275) builds o[k] = d - beta * sum_{j: A[sigma_j][k] = 0} a_j for every k not in sigma.
  `fastregion.region` puts every o[k] into the isolation list in stage 1 (names `o{k}`, positivity and off-support
  cells, line ~354) and again in stage 2 (`names2 = names + ["D", "Delta"]`, so the stage-2 critical set is
  {a_i, o_k, D, Delta_{m-1}}). `decide_point` tests o_k/D < 0 at every sample.
- `region2d.positivity_intervals` (line ~257) isolates the roots of all a_i **and all o_k**. The positivity set it
  returns is {a_i > 0, o_k < 0}. In `region2d.region2d` (line ~288) `pos_polys` = a_i together with o_k, and these
  are prepended to the projection polynomials in `intp`. So every real o_k root is a wall of the beta-decomposition
  (a vertical line in the (beta, t) plane). The docstring and the certificate `method` field both read
  "CAD on {a_i, o_k, D, Delta_{m-1}, L}". `symbolic.q_region` also includes `off_polys` (lines 109-115, 163, 172).
- The n = 7 proper-support classification was run with **region2d**, not fastregion.region:
  `n7/run_n7_proper_supports.py` lines 17 and 73, and `certificates/n7_proper_supports_summary.json` field `method`.
  The engine includes o_6 either way.
- Certificate check: `notes/referee_external/fixes_2026-11-04/prover/certificate_check.py`
  (output `certificate_check.json`, stdout `certificate_check_stdout.txt`; about 72 s).
  - Report example = record pattern 30: verdict stable_somewhere. Its `positivity_critical` contains the irreducible
    factor [1,-3,3,-3,2,-1] (= -o_7) with isolating interval [2, 19/9]. Its positivity cell is (25/12, 9/4), which
    starts above beta*. 209/100 is in that cell and 52/25 is in no positivity cell. The stable cell
    beta in (25/12, 17/8), t > 48/5 contains (209/100, 100), and no stable beta-cell contains 52/25.
    This agrees with (a).
  - All 7,744 stable_somewhere records: o_6 was recomputed from the stored adjacency with `fastregion.beta_data`.
    Each of its 2,496 real roots above 2, irreducible factor by factor, was matched to an entry of the record's
    `positivity_critical` (same factor up to sign, overlapping isolating interval containing the root). **0 missing.**
  - The 4,544 never records store no cell decomposition (`n7/run_n7_proper_supports.py` line 79), so for them the
    claim rests on reading the code: the same `region2d` path with o_6 in `positivity_intervals` and `intp`.
- Conclusion: the omission is in the paper's prose only. The n = 7 proper-support results are unaffected
  (12,288 pairs: 7,744 stable somewhere, 4,544 never, 0 escalations). Full supports have no o_k, so the full-support
  counts and thresholds are unaffected too.

## (c) Text edit (proposed, not applied)

In `notes/paper/paper.tex` lines 423-425, change

    the complement of the curves $\{a_i=0\}$, $\{D=0\}$, $\{\Delta_{m-1}=0\}$ and the legality boundary
    $L(\beta,t)=t(\beta-2)-1=0$ in the open legal region, and that

to

    the complement of the curves $\{a_i=0\}$, $\{o_k=0\}$ ($k\notin\sigma$), $\{D=0\}$, $\{\Delta_{m-1}=0\}$ and the
    legality boundary $L(\beta,t)=t(\beta-2)-1=0$ in the open legal region (for $\sigma=[n]$ there are no $o_k$ and the
    list reduces to the $a_i$, $D$, $\Delta_{m-1}$ and $L$), and that

With this change, "these polynomials" in the projection-set clause includes the t-free o_k. That matches region2d,
which adds the a_i and o_k to the beta-isolation list. The list then agrees with Lemma 3.2's P_t and with Z (line 437).

Optional consistency edit (line 429-430): "the vertical line above a positivity root" could read "the vertical line
above a root of some $a_i$ or $o_k$". In region2d a "positivity root" already includes the o_k roots, but the paper
never defines the term.

Other occurrences: a recursive search of notes/ (*.md, *.tex, *.txt) and README.md for "connected components of the
complement", "complement of the curves" and the list "a_i=0" finds only this paragraph. The one other hit is
`notes/paper/paper_draft_before_check.tex` lines 384-385, an archived pre-check draft with the same omission. It is
not the built paper and should be left as is. `notes/lemma_Bprime_2d.md` already includes the o_k.

Credit: the relayed message says the credit as written is fine, so no new acknowledgement text is proposed.

Run commands:

    cd /home/user/CLTN && nice -n 10 python3 notes/referee_external/fixes_2026-11-04/prover/offsupport_check.py
    cd /home/user/CLTN && nice -n 10 python3 notes/referee_external/fixes_2026-11-04/prover/certificate_check.py
