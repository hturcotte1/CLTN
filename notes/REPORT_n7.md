# The target-free clique conjecture at n = 7: a complete exact map, and what it reveals for all n

Date: 2026-10-15. Companion to `notes/REPORT_n6.md` (definitions, conventions, the n = 6 theorem and
its audit). Conventions are CLAUDE.md section 2 throughout; everything below is exact (integers,
rationals, exact root isolation, Hurwitz determinants); decimals appear only as labelled approximations.

> **Status.** The 124 seven-node counterexamples are **verified**: hand-verified by the user on 2026-10-15
> and confirmed by four independent exact implementations written without the project's code
> (`notes/verify_124_workflow.json`; 124/124 at every one of 742 legal points, 0 undecided). The n = 6
> bases were signed off on 2026-10-14. Graphs on n >= 8 nodes from the infinite families remain
> candidates in the sense of CLAUDE.md section 3 (the theorems are refereed exact algebra; the instances
> were certified but not hand-verified).

## 0. Headline results

| # | result | status |
|---|---|---|
| 1 | **Complete finite-eps classification at n = 7.** Among all 882,033,440 digraphs on 7 nodes, exactly **124 non-clique graphs** have a full support that is a stable positive fixed point at some legal (eps, delta), each with an exact two-parameter region certificate; K7 is the only clique. 0 escalations. | CERTIFIED and VERIFIED |
| 2 | The 124 are **the 123 limit-stable graphs plus one graph, `&FNFL~[udo?`, invisible in the eps -> 0 limit** (spectral abscissa exactly -1/2; stable only for a bounded range of t at each beta). | CERTIFIED and VERIFIED |
| 3 | **The limit does not determine the finite-eps map**: "stable at finite eps => limit-stable" (Lemma M) is false, with exact certificates at m = 3 (pure stability) and n = 7 (with positivity). 40 of the 124 have stable regions that bulge past beta_c and are bounded above in t (Lemma AS: first-order coefficient R* > 0). | PROVED / CERTIFIED |
| 4 | **Lemma R (ray reduction).** The full support is stable at some legal parameter iff it is stable on the legality curve t(beta - 2) = 1; an exact one-dimensional test that makes the finite-eps map affordable (16 CPU-hours at n = 7). | PROVED, audited |
| 5 | **Theorem F.** Geneson's graph plus a bidirected clique class of k nodes is limit-stable for every k >= 0: **the conjecture fails for every n >= 6**, with sup q = 1.1664 (n = 6), 1.2747 (n = 7), the plastic number 1.3247 (n = 8..14), then decreasing to 1. | PROVED, refereed |
| 6 | **Theorem E.** An infinite family E_a (n = 7 + a) of finite-eps-only counterexamples: spectral abscissa exactly -1/2 for every a, stable positive at legal points with eps -> 0, beta -> 2+. `&FNFL~[udo?` is E_0 = the n = 6 hit idx 1387528 with one node doubled. | PROVED, refereed |
| 7 | **Growing thresholds.** Families Q (n = 7 + k, q -> 1.8442), Q2 (n = 7 + 4k, q -> 3), Y (n = 8 + 5k, q >= 3.8 for k >= 363); exact lower bounds at n = 8: 1709/1001, n = 9: 17691/10010, n = 10: 18418/10001, n = 47: 25000/10001. **Theorem U: q_n* -> infinity** (q_n* >= N - 1 for all n >= n0(N), every N). | PROVED, refereed |
| 8 | **Exact thresholds.** q_6* = root of 797q^6 - 1686q^5 + 2159q^4 - 3004q^3 + 2499q^2 - 2158q + 1521 in (561/481, 7/6): at n = 6 the conjecture holds for every legal parameter with q >= q_6* and fails for every q in (1, q_6*). **q_7* = the plastic number**, root of q^3 - q - 1 (1.3247): for n <= 7 the conjecture holds at every legal parameter with q >= rho and fails for every q in (1, rho) (section 4). | CERTIFIED (q_6* final; q_7* candidates) |
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
  is stable for a window of moderate t only (the "bulge"). Exactly characterised (Lemma AS, exact R*
  for all 124, `certificates/n7_bulge_Rstar_exact.json`): a graph bulges iff R* > 0 AND its fixed point
  is positive at beta_c; 51 have R* > 0 but 11 of those lose positivity before beta_c, and every bounded
  cell starts exactly at beta_c. The six plastic-number graphs do not bulge (their limit interval ends at
  a positivity root). Exact witnesses for two bulging graphs:
  `certificates/n7_limit_hits_finite_eps_bulge_witnesses.json`.
* 21 of the 123 contain an n = 6 counterexample as an induced subgraph; none contains K6.
* `&FNFL~[udo?` (28 arcs): `7 chi_perp = (l + 1)(l^2 + l + 1)(7l^3 + 14l^2 + 14l + 5)`, so
  alpha = -1/2 exactly; never limit-stable; stable positive at eps = 1/49, delta = 33/1568 (q = 99/98)
  and eps = 1/36, delta = 17/588; unstable for larger t at the same beta; region bounded in t at every
  beta (cell at beta = 53/26: t from 26 up to a root t_1 in (53, 54)). It is the n = 6 hit idx 1387528 with one node doubled
  (closed-twin blow-up), and the factor l^2 + l + 1 is a moving eigenvalue pair on the boundary line.
  `certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json`.

## 3b. All stable non-clique supports at n = 7, not only full ones

By Theorem F a proper support is stable iff its induced subgraph has a stable positive full support and
the off-support condition holds; the only non-clique induced subgraphs with a stable positive full support
have 6 nodes (the three n = 6 counterexamples). Enumerating each base plus a seventh node with any of the
4096 in/out patterns (the bases have trivial automorphism groups, so the 12,288 pairs are pairwise
non-isomorphic) and certifying the size-6 support with region2d including the off-support polynomial:
7,744 pairs are stable somewhere in the legal plane, 4,544 never, 0 escalations
(`certificates/n7_proper_supports_all.jsonl`, `_summary.json`). A seventh node that is a target of the
support is never stable (192 of 192), one receiving nothing always is (192 of 192), and no proper-support
region is bounded in t. Together with the 124 full supports this is the complete list of stable non-clique
supports at n = 7.

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

**Result (derived from the certificates and confirmed by the independent exact `qstar --sharpen` run,
`certificates/qstar/qstar_n7_finite.json`: upper bound = lower bound, exact, no undecided cells).**
`q_7* = rho`, the plastic number, the real root of `q^3 - q - 1` (about 1.3247), a supremum that is not
attained: six graphs approach it as t -> infinity with beta -> 1 + rho. Derivation: q < beta - 1 on
every cell; on the 258 stable cells unbounded in t the supremum is beta_hi - 1 exactly and the right
endpoints are the limit endpoints, whose maximum is 1 + rho (the limit stage, 0 mismatches); on the 40
bounded cells the outer isolating bound gives beta_hi <= 16/7, so q < 9/7 < rho. Four of the six
plastic-number graphs have limit intervals starting at q = 1, so the candidates' intervals cover (1, rho)
without a gap. Hence, conditional on hand-verification of the candidates and combined with the n <= 6
maps and Theorem F (reduction): **for n <= 7 the conjecture holds at every legal parameter with q >= rho
and fails for every q in (1, rho)** (the conditionality on hand-verification is now discharged); Geneson's bound n - 2 = 5 is replaced by 1.3247. For comparison the
infinite family F_k alone reaches 1.2747 at n = 7.

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
  **Theorem U (refereed).** q_n* -> infinity: for every N there is n0(N) with q_n* >= N - 1 for all
  n >= n0(N) (arms-forest cores realising prescribed limit polynomials plus one sink class; certified
  thresholds q_n* >= 2, 3, 4, 9/2, 23/5 for n >= 1812, 219, 657, 2475, 4909; explicit graphs at n = 28 and
  41). Geneson's q >= n - 2 remains the only upper bound; the growth rate of q_n* is open. A catalogue of 9 (n = 6) + 617 (n = 7) clique-class families, each proved for
  all k >= 1 by exact root isolation in k, is in `families_catalogue_n6.json` / `_n7.json`
  (refereed: all 626 entries re-proved independently for all real k >= 1).
* **Refuted on the way.** "Every limit polynomial has a root on Re = -1/2" (55 of 135 do not), and
  Conjecture P (positivity of (D^T + rI)^-1 1 for r in (0, 1/2) forces an eigenvalue with Re <= -1/2):
  33 catalogue bases, e.g. `&FMAlW\d\y?`, have all eigenvalues of D strictly right of -1/2 and a positive
  vector at r = 49/100; the joins G v K_k are limit-stable for every k. **Settled:** sup_n q_n* = infinity (Theorem U), and at a FIXED legal (eps, delta) the stable
  non-clique support size IS unbounded (Theorem A, refereed by an independent second proof): at
  eps = 1/21, delta = 11/210, where Geneson's clique-class family dies at k = 41, the family R_k built on
  the 7-node counterexample `&F^U|m^QTq?` is a stable positive non-clique full support for every k >= 0;
  every legal point with beta < 4.9607 carries such a family (unconditionally), and conditionally on the
  Theorem U lemmas every legal point does. **Open:** the growth rate of q_n* and of the smallest family
  member.

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
* Human sign-off: n = 6 hits 2026-10-14; the 124 n = 7 graphs hand-verified 2026-10-15 and machine-
  confirmed the same day (four independent verifiers, section above). One legal witness per graph for
  re-checking: `certificates/n7_candidates_legal_witnesses.json` (15 of the stored first witnesses in the
  limit file are illegal points, flagged as such; use the legal-witness file).

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
3. **Theory**: characterise the limit-stable graphs (Conjecture P is refuted, so a different invariant
   is needed); the growth rate of q_n*; extend the Lean formalisation (Theorems F, E, U, A).

## 11. Machine-checked statements (Lean 4 + Mathlib)

`notes/lean/CtlnCert/` (toolchain lean4 v4.35.0, Mathlib v4.35.0; `lake build` passes; axiom audit:
only propext, Classical.choice, Quot.sound; no `sorry`): 26 theorems, among them Remark F's affine
identity `K(l beta, l t) = (1 - l) I + l K(beta, t)` over any commutative ring; the characteristic-polynomial
scaling and the eigenvalue decomposition along rays; positive stability propagates down rays for
0 < l <= 1; the intermediate-value step with the explicit unique root `(t + sqrt(t^2 + beta t))/(beta t)`
in (2/beta, 1); the CLAUDE.md section 2 closed form `I - W = -delta I + (1 + delta) J - (eps + delta) A^T`
(with the transpose, as the code uses it); `I - W = eps K(J - I - A^T, beta, t)`; `q > 1 <-> t(beta - 2) > 1`;
Hurwitz stability of `-I + W` iff positive stability of K; and the forward direction of Lemma R in
eps/delta language. A second batch (31 more theorems) machine-checks the pipeline's prune lemmas: Lemma A.2
(positive stability forces det K > 0; under Hurwitz stability x > 0 iff adj(K) 1 > 0, so positivity is
t-free), the e_2 edge bound including its spectral step (via Mathlib's charpoly coefficients as sums of
principal minors and Gauss-Lucas: Hurwitz stability at a legal point forces more bidirected than
non-adjacent pairs and at least m(m - 1)/2 + 1 arcs), and Lemma D (a dominated node forces x_j = 0 at
every fixed point, also on induced subgraphs). A third batch adds Stodola's criterion (a monic real
polynomial with all roots in the open left half-plane has all coefficients positive; hence every
coefficient of det(mu I + K) is positive when K is positive stable, the statement behind the curve test's
fast path) and Remark J.4' (two non-adjacent in-twins in the support give -I + W the eigenvalue delta,
so it is never Hurwitz stable, for every eps and delta > 0). 87 theorems in all, no `sorry`, axiom audit
clean. The converse of Lemma R (openness of stability) is not formalised.
