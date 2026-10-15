# The target-free clique conjecture at n = 7: a complete exact map, and what it reveals for all n

Date: 2026-10-15. Companion to `notes/REPORT_n6.md` (definitions, conventions, the n = 6 theorem and
its audit). Conventions are CLAUDE.md section 2 throughout; everything below is exact (integers,
rationals, exact root isolation, Hurwitz determinants); decimals appear only as labelled approximations.

> **Status of the n >= 7 graphs.** Every 7-node (and larger) graph named below as a counterexample is a
> **candidate pending human hand-verification** (CLAUDE.md section 3). The n = 6 bases (Geneson's graph,
> idx 1387528, idx 1424629) were hand-verified and signed off on 2026-10-14. The theorems about infinite
> families are exact algebra on those verified bases and were refereed by independent agents; their
> n >= 7 instances inherit the candidate label.

## 0. Headline results

| # | result | status |
|---|---|---|
| 1 | **Complete finite-eps classification at n = 7.** Among all 882,033,440 digraphs on 7 nodes, exactly **124 non-clique graphs** have a full support that is a stable positive fixed point at some legal (eps, delta), each with an exact two-parameter region certificate; K7 is the only clique. 0 escalations. | CERTIFIED (candidates) |
| 2 | The 124 are **the 123 limit-stable graphs plus one graph, `&FNFL~[udo?`, invisible in the eps -> 0 limit** (spectral abscissa exactly -1/2; stable only for a bounded range of t at each beta). | CERTIFIED (candidates) |
| 3 | **The limit does not determine the finite-eps map**: "stable at finite eps => limit-stable" (Lemma M) is false, with exact certificates at m = 3 (pure stability) and n = 7 (with positivity). 40 of the 124 have stable regions that bulge past beta_c and are bounded above in t (Lemma AS: first-order coefficient R* > 0). | PROVED / CERTIFIED |
| 4 | **Lemma R (ray reduction).** The full support is stable at some legal parameter iff it is stable on the legality curve t(beta - 2) = 1; an exact one-dimensional test that makes the finite-eps map affordable (16 CPU-hours at n = 7). | PROVED, audited |
| 5 | **Theorem F.** Geneson's graph plus a bidirected clique class of k nodes is limit-stable for every k >= 0: **the conjecture fails for every n >= 6**, with sup q = 1.1664 (n = 6), 1.2747 (n = 7), the plastic number 1.3247 (n = 8..14), then decreasing to 1. | PROVED, refereed |
| 6 | **Theorem E.** An infinite family E_a (n = 7 + a) of finite-eps-only counterexamples: spectral abscissa exactly -1/2 for every a, stable positive at legal points with eps -> 0, beta -> 2+. `&FNFL~[udo?` is E_0 = the n = 6 hit idx 1387528 with one node doubled. | PROVED, refereed |
| 7 | **Growing thresholds.** Families Q (n = 7 + k, q -> 1.8442), Q2 (n = 7 + 4k, q -> 3), Y (n = 8 + 5k, q >= 3.8 for k >= 363); exact lower bounds on the largest counterexample q at n = 8: 1709/1001, n = 9: 17691/10010, n = 10: 18418/10001, n = 47: 25000/10001. Conjecture: sup_n q_n* = infinity. | PROVED, refereed |
| 8 | **Exact thresholds.** q_6* = root of 797q^6 - 1686q^5 + 2159q^4 - 3004q^3 + 2499q^2 - 2158q + 1521 in (561/481, 7/6): at n = 6 the conjecture holds for every legal parameter with q >= q_6* and fails for every q in (1, q_6*). At n = 7 the largest limit q is the plastic number (6 graphs); the exact finite-eps supremum q_7* is in section 4. | CERTIFIED (q_6* final; q_7* candidates) |
| 9 | Not scooped as of 2026-10-14; an unrefereed Zenodo preprint independently confirms the n <= 5 exclusion and the value of q_6* (section 9). | checked |
| 10 | Every witness (the three n = 6 graphs and all 124 n = 7 candidates) is a **nondegenerate** CTLN at its witness parameters, certified exactly over all 2^n - 1 supports. | CERTIFIED |

## 1. Setting

For a digraph on n nodes with `A[i][j] = 1` iff `i -> j`, and parameters eps, delta > 0, the CTLN
has full-support Jacobian data `I - W = eps (C + tJ)`, `C = I + beta N`, with `N[i][j] = 1` iff
`j -/-> i`, `t = (1 - eps)/eps`, `beta = 1 + delta/eps`, `q = delta(1 - eps)/eps = (beta - 1) t/(t + 1)`;
legal iff `q > 1` iff `t(beta - 2) > 1`. The full support is a stable fixed point iff `K = C + tJ` is
positive stable (all eigenvalues in the open right half plane; decided by Hurwitz determinants of
`det(mu I + K) = phi(mu) + t omega(mu)`) and `x = adj(K) 1 / det K > 0`; `adj(K) 1 = a(beta)` is
t-free (Lemma A), and stability forces `det K > 0`, so positivity is `a_i(beta) > 0` for all i
(Lemma A.2). The conjecture says: full support stable iff the graph is the complete bidirected graph.
By Theorem F of `notes/lemma_F_reduction.md`, the conjecture at n holds iff no non-complete graph on
k <= n nodes has a positive stable full support, so the full-support question at every k <= n is the
whole question.

The eps -> 0 limit (Lemma I): limit-stable at beta iff the spectral abscissa `alpha(A|perp)` of the
adjacency compressed to the hyperplane orthogonal to the all-ones vector is `< -(1 - 1/beta)`; the
limit-stable set is `(2, beta_c)`, `beta_c = 1/(1 + alpha)`, nonempty iff `alpha < -1/2`.

## 2. The accounting chain

Every one of the 882,033,440 digraphs on 7 nodes (OEIS A000273(7), reproduced twice by nauty's
geng|directg) is dropped by a proved necessary condition valid at every legal parameter, or carries
an exact two-parameter certificate. Every stage was validated to diff zero against the n = 6 map
before use, and audited by independent agents (`notes/audit_n7_finite_workflow.json`).

| stage | condition (proved where) | dropped | remaining |
|---|---|---|---|
| 0 | enumeration | | 882,033,440 |
| 1 | dominated node (Curto Def 1.1; Lemma D) | 764,800,518 | 117,232,922 |
| 2 | no beta > 2 with all a_i > 0, exact Descartes in C (Lemma A.2) | 107,600,789 | 9,632,133 (8,845,537 certified positive + 786,596 undecided) |
| 3 | undecided lines resolved by exact root isolation | 723,010 | 8,909,123 positivity survivors |
| 4 | at most 21 edges (e_2 bound: legal stability forces more bidirected than non-adjacent pairs) | 5,135,775 | 3,773,348 |
| 5 | never stable on the legality curve (Lemma R, exact 1-D) | 3,773,183 | 165 |
| 6 | two-parameter CAD on the proven boundary set {a_i, D, Delta_6, L} (region2d, cross-checked, replayed) | 40 stable-but-not-positive-together | **125 = K7 + 124** |

Costs on a 4-core machine: stage 0-2 49 min, stage 3 1.5 h, stage 5 4.1 h (16.4 CPU-h, 6.6 ms/line),
stage 6 58 min. Validation: stage 2-3 reproduce the n = 6 strict positivity set exactly
(28,350, diff zero both ways); stage 5 on all 9,608 five-node digraphs (95 stable, every one with an
exact witness) and on the 28,350 six-node positivity survivors returns exactly the four stable
graphs; stage 6 is the engine that produced the audited n = 5 and n = 6 maps.

Files: `n7/filter7p.c`, `n7/resolve_U.py`, `n7/curve_test.py`, `n7/run_n7_curve.py`,
`n7/run_n7_region2d.py`, `n7/aggregate_n7.py`; certificates `certificates/n7_finite_summary.json`,
`certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl` (all 124 with full region certificates),
`certificates/n7_region2d/` (all 165 CAD records), `certificates/n7_curve_hits.jsonl`,
`n7/n7_pos_all_PU.txt.gz` (stage-2 output, 9,632,133 lines). Runbook: `n7/REPRODUCE_n7_finite.md`.

## 3. The 124 graphs

* Edge counts 23 to 31 of 42 possible arcs (23: 7, 24: 7, 25: 18, 26: 10, 27: 6, 28: 22, 29: 34, 30: 17, 31: 3);
  every one has at least 22 arcs, as the e_2 bound requires, and all are domination-free.
* 91 D-equivalence classes (same characteristic polynomial and in-walk counts, Lemma J): 71 singletons,
  16 pairs, two classes of 4, one of 5, one of 8. 64 of the 124 carry a bidirected in-twin pair, so
  many are one inert-edge move apart (bidirected in-twin switching preserves every fixed point, Lemma J.4).
* 123 are limit-stable with limit interval (2, beta_c) (or starting above 2 for 6 of them, whose
  positivity set starts strictly above beta = 2); 69 distinct minimal polynomials of the limit critical
  value q_c, the most frequent being q^5 + 2q^4 + q^3 - 2q^2 - 2q - 1 (17 graphs; the same quintic as
  the n = 6 hit idx 1387528), and the plastic number q^3 - q - 1 (6 graphs, the largest q_c).
* 40 of the 124 have a stable cell bounded above in t: at beta slightly above beta_c the full support
  is stable for a window of moderate t only (the "bulge"). Exact witnesses for two of them:
  `certificates/n7_limit_hits_finite_eps_bulge_witnesses.json`.
* 21 of the 123 contain an n = 6 counterexample as an induced subgraph; none contains K6.
* `&FNFL~[udo?` (28 arcs): `7 chi_perp = (l + 1)(l^2 + l + 1)(7l^3 + 14l^2 + 14l + 5)`, so
  alpha = -1/2 exactly; never limit-stable; stable positive at eps = 1/49, delta = 33/1568 (q = 99/98)
  and eps = 1/36, delta = 17/588; unstable for larger t at the same beta; region bounded in t at every
  beta (cell at beta = 53/26: t in (26, 53)). It is the n = 6 hit idx 1387528 with one node doubled
  (closed-twin blow-up), and the factor l^2 + l + 1 is a moving eigenvalue pair on the boundary line.
  `certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json`.

## 4. Thresholds

**n = 6 (final).** The supremum of q over all stable non-clique full supports at n = 6 is
`q_6* = root of 797q^6 - 1686q^5 + 2159q^4 - 3004q^3 + 2499q^2 - 2158q + 1521 in (561/481, 7/6)`
(about 1.1664), attained in the limit t -> infinity by Geneson's graph and its twin partner
(`certificates/qstar/qstar_n6.json`, tool `n7/qstar.py`). Together with the complete n <= 6 maps and
Theorem F (reduction): **for n <= 6 the conjecture holds at every legal parameter with q >= q_6* and
fails for every q in (1, q_6*)**; Geneson's bound n - 2 = 4 is replaced by 1.1664.

**n = 7 (candidates).** Lower bound from the limit: the plastic number rho = 1.3247 (root of
q^3 - q - 1, six graphs). The finite-eps supremum q_7* over the 124 regions is computed exactly by
`n7/qstar.py --sharpen` on the two-parameter certificates: see `certificates/qstar/qstar_n7_finite.json`
and the line below.

*(q_7* exact value: computation running at the time of this draft; this line is replaced when `certificates/qstar/qstar_n7_finite.json` exists.)*

## 5. Why the limit is not the whole story

**Lemma M is false.** The natural hope that stability at a legal finite t implies limit stability at
the same beta fails: for K = I + beta N + tJ the stable t-set at fixed beta need not be upward closed.
Exact certificates (`notes/lemma_M_refuted.md`): m = 3, A = 011/001/100, beta = 22/5, stable at
t = 11/25 (legal), unstable at t = 1 and 100; and with positivity at n = 7 (the graph above and the
bulging limit candidates). The general-matrix analogue is classical (conditional stability of a
root locus); legality is what makes m = 2 immune (the one-edge pair is stable exactly for illegal t).

**Identity (E) and Lemma AS.** `K = beta (cJ - A^T - gamma I)` with `c = 1 + t/beta`,
`gamma = 1 - 1/beta`, so stability is `nu_A(c) := min Re spec(cJ - A) > gamma`, legality is
`gamma > c - sqrt(c^2 - c)` (the stability boundary of a single one-way edge), and the large-t expansion
`nu_A(c) = -alpha + R*/c + O(c^-2)` with `R* = Re(1/s'(z0))`, `s(z) = 1^T (zI + A)^-1 1`, decides the
shape near (beta_c, infinity): R* > 0 means the stable region bulges beyond beta_c and is bounded in t
there. The finite-eps-only graphs are those with alpha = -1/2 exactly and R* > 1/8 (Lemma B).

**The e_2 bound.** The sum of the 2 x 2 principal minors of K is
`B(1 + 2t) + Nn(1 - beta)(1 + beta + 2t) + O(1 - t(beta - 2))` over bidirected, non-adjacent and
one-way pairs; at a legal point the last two terms are negative, so stability forces B > Nn, i.e.
more than half of all node pairs bidirected: at least m(m - 1)/2 + 1 arcs.

## 6. Lemma R: the finite map as a one-dimensional problem

`K(lambda beta, lambda t) = (1 - lambda) I + lambda K(beta, t)` (Remark F), so stability propagates
down every ray through the origin until the ray meets the legality curve t(beta - 2) = 1 (the
intermediate value theorem on `lambda t (lambda beta - 2) - 1`), and conversely a stable curve point
has stable legal neighbours. Hence **the full support is stable at some legal parameter iff
K(2 + 1/t, t) is positive stable for some t > 0**. Along the curve `tK = tI + (2t + 1) N + t^2 J` is an
integer polynomial matrix; the verdict changes only at positive roots of `det` and of the Hurwitz
minor Delta_{m-1} (Lemma B with t as the parameter), so one exact sample per root-free cell decides
it (`n7/curve_test.py`, 35-45 ms per 7-node graph; audited; validated on all 5-node digraphs and all
6-node positivity survivors). Since positivity is t-free, the stable positive regions are unions of
ray segments from the curve, and the two-parameter CAD is needed only for the curve-stable graphs.

## 7. Infinite families (`notes/research_wf3/theory/notes_theory.md`)

* **Closed-twin blow-ups (Lemmas W, X; refereed).** Replacing a node by a bidirected clique class of
  k copies with the same in- and out-neighbourhoods preserves the full-support fixed-point data by
  explicit characteristic-polynomial identities, with positivity independent of k when the class
  sends to everything.
* **Theorem F (refereed).** F_k = Geneson's graph + a bidirected clique class of size k (n = 6 + k):
  alpha < -1/2 for every k (Hurwitz minors 16(k + 1), 64(k + 1)(k + 16), 8192(k + 1)^2 (3k + 16)),
  positive on (2, 1 + rho). So counterexamples exist for every n >= 6. F_1 is one of the 123.
* **Theorem E (refereed).** E_a = blow-up of `&FNDL~[uds?` with class sizes (1, 1, 2, a, 1, 1, 1),
  n = 7 + a: `Phi = z (z^2 + z + 1)((a + 7) z (z^2 + z + 1) + 2)`, alpha = -1/2 exactly for all a,
  R* = 1/4; by Lemma B every E_a is a finite-eps-only counterexample (witnesses at q = 421/420 for
  a = 0..5). E_0 is `&FNFL~[udo?`.
* **Theorems Q, Q2, Y (refereed).** Sink-class blow-ups of three 7- and 8-node cores give limit-stable
  families with q -> 1.8442 (n = 7 + k), q -> 3 (n = 7 + 4k) and q >= 3.8 for n = 8 + 5k, k >= 363.
  So the largest counterexample q grows without any known bound; Geneson's safe region q >= n - 2 is
  the only upper bound. A catalogue of 9 (n = 6) + 617 (n = 7) clique-class families, each proved for
  all k >= 1 by exact root isolation in k, is in `families_catalogue_n6.json` / `_n7.json`
  (refereed: all 626 entries re-proved independently for all real k >= 1).
* **Refuted on the way.** "Every limit polynomial has a root on Re = -1/2" (55 of 135 do not), and
  Conjecture P (positivity of (D^T + rI)^-1 1 for r in (0, 1/2) forces an eigenvalue with Re <= -1/2):
  33 catalogue bases, e.g. `&FMAlW\d\y?`, have all eigenvalues of D strictly right of -1/2 and a positive
  vector at r = 49/100; the joins G v K_k are limit-stable for every k. **Open:** sup_n q_n* = infinity;
  and, at a FIXED legal (eps, delta), whether the stable non-clique support size is unbounded (the
  clique-class family F_k is stable at eps = 1/21, delta = 11/210 only for k <= 38).

## 8. Verification status

* Engine and n <= 6 map: audited (`notes/REPORT_n6.md` section on audits).
* n = 7 pipeline: two independent soundness audits of the prefilter and resolver (979,800 exact
  comparisons, all 786,596 undecided lines recomputed independently, byte-identical), one of Lemma R,
  the curve test and the edge bound (independent all-minors arrangement on 3,129 graphs, 6,000 random
  exact legal points on never-graphs), two exact re-verifications of the FNFL certificate; no blocking
  issue (`notes/audit_n7_finite_workflow.json`, artifacts `notes/audit_n7_finite/`).
* Theory: Lemma D1/D2, Lemmas W/X, Theorem F, Theorem E + Lemma B, Theorems Q/Q2/Y refereed
  (`notes/research_n7_parallel_workflow.json`, `notes/research_wf3/`); catalogue, claim 16 and the
  the 626-entry catalogue, claim 16 and the q* tool refereed as well (`notes/referee_n7_theory_workflow.json`,
  artifacts `notes/referee_wf4/`): every claim holds with high confidence; the referees' write-up notes are
  logged in ANOMALIES.md (2026-10-15). All family statements rest on Lemma I.3/I.4 (limit-to-finite-t
  transfer), refereed in the n = 6 phase.
* Nondegeneracy: `n7/nondegeneracy_check.py` certifies at a witness point of every hit that all 2^n - 1
  supports are nonsingular with no zero coordinate and no zero off-support drive
  (`certificates/nondegeneracy_witnesses.json`).
* Human sign-off: n = 6 hits signed off 2026-10-14; the 124 n = 7 graphs (and the family instances)
  await hand-verification. Suggested order: `&FNFL~[udo?` (two rational witnesses, section 3), then
  one representative per D-class from `certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit_summary.json`.

## 9. Literature position (as of 2026-10-14)

No published or posted work classifies n = 6, treats n = 7, exhibits finite-eps-only counterexamples,
or improves Geneson's threshold. Geneson arXiv:2607.21396 is v1 only. An anonymous, unrefereed
Zenodo/GitHub preprint ("Six Is Minimal", 2026-10-01) proves n <= 5 has no counterexample, uses
Geneson's graph as the 6-witness, and gives a family on 5r + 6 nodes; it reports counts that match
A000273(1) + ... + A000273(5) = 9,846 and sum 2^(m(m-1)) = 1,052,741, i.e. all digraphs on up to five
nodes. A claims diff (their repository read, not run; `notes/referee_wf4/lit_diff/`) finds the two works
consistent: their conventions equal ours up to a transpose of the nonedge matrix; their n <= 5 exclusion
equals our n = 5 map; their Geneson-class region gives 1 + u* = 1.16635502785, the same sextic as our
q_6*; their family on 5r + 6 nodes is a clique-class blow-up of the same seed with q -> 1 (a Lemma W
instance). They do not classify n = 6, treat n >= 7, exhibit finite-eps-only examples, or state Theorem F
or the growing-q families; they do give an explicit rational parameter schedule for their family and
global nondegeneracy certificates, which prompted our `nondegeneracy_check.py`. We cite it as an unrefereed
preprint, "independently obtained / consistent with", and claim no priority either way.

## 10. What is next

1. **Hand-verification** of the 124 (section 8); then the status labels in every certificate change.
2. **n = 8 limit map**: the overflow-guarded filter `n7/filter8.c` is validated byte-identical on all
   n <= 7 data; the full pass is 500-800 CPU-hours (runbook `n7/REPRODUCE_n8_limit.md`) and needs
   rented compute. An incidental n = 8 limit candidate was already found in a random dense sample.
3. **Theory**: prove or refute sup_n q_n* = infinity; characterise the limit-stable graphs (Conjecture P);
   Lean formalisation of Lemmas A, R and Theorem F (Lean is not installed here).
