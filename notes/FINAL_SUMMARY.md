# Final summary: the target-free clique conjecture, n <= 8 and beyond

Date: 2026-10-28 (first version 2026-10-15). Start here; details in `notes/REPORT_n6.md`,
`notes/REPORT_n7.md`, `notes/REPORT_n8.md`, the lemma notes and `certificates/INDEX.md`. Every number
below is either backed by a re-runnable exact file (named in the report that owns the number) or explicitly
an estimate or a display approximation.

## 1. The conjecture and what is now known

Curto-Morrison Open Question 10: in a CTLN, the supports of stable fixed points are exactly the
target-free cliques. Geneson (arXiv:2607.21396, July 2026) gave a 6-node counterexample and proved the
conjecture for q >= n - 2 - (n - 3) eps/2, where q = delta(1 - eps)/eps > 1 is the legality parameter.

| n | result | threshold (sup of q over counterexamples) | status |
|---|---|---|---|
| <= 5 | no counterexample at any legal parameter (also reproduced the published n = 5 classification at eps = 1/4, delta = 1/2 with zero mismatches) | none | certified |
| 6 | exactly 3 non-clique stable full supports among 1,540,944 digraphs (Geneson's graph, its bidirected in-twin partner, one new graph) | q_6* = root of 797q^6 - 1686q^5 + 2159q^4 - 3004q^3 + 2499q^2 - 2158q + 1521, about 1.1664 | certified, audited, hand-verified |
| 7 | exactly 124 among 882,033,440 digraphs, at every legal parameter; 123 visible in the eps -> 0 limit, one (`&FNFL~[udo?`) only at finite eps | q_7* = the plastic number, root of q^3 - q - 1, about 1.3247, not attained (two independent exact derivations) | certified, audited, hand-verified (a check the user recorded as not fully confident), confirmed by four independent implementations |
| 7, proper supports | every stable non-clique proper support classified: 7,744 of the 12,288 (graph, size-6 support) pairs built on the three n = 6 counterexamples are stable somewhere, 4,544 never | (inherits q_6*) | certified |
| 8, **limit map only** | exactly **9,133** non-complete graphs whose full support is limit-stable (eps -> 0), among the 801,558,826,511 digraphs with >= 29 arcs (every other digraph on 8 vertices is excluded by the Lean-checked e_2 bound); K8 the only complete one; 0 escalations. 61 are closed-twin lifts of n = 7 counterexamples, 9,072 are not; 6,720 contain no smaller counterexample. **The finite-eps map at n = 8 is not done** (the 417,603,634 boundary lines with alpha >= -1/2 were not stored, so finite-eps-only counterexamples of the `&FNFL~[udo?` type are outside the run; 15 such eight-node graphs are known from Lemma B, E_1 among them, with certified legal witnesses and no hand check) | q_8*(limit) = root of a degree-15 polynomial, 13109q^15 - 124967q^14 + ... - 9045013, in [459224485/268435456, 2755346911/1610612736], about 1.7107, not attained; all 9,133 upper ends below 2, 16 above 3/2; finite-eps q_8* >= this, not determined | certified; **verified by machine only** (two independent campaigns, 9,133/9,133, 0 fatal or major defects); **no human hand check** (declined by the user 2026-10-28) |
| all n >= 6 | counterexamples exist for every n (Theorem F); an infinite family invisible in the limit (Theorem E); the critical q is unbounded: q_n* -> infinity (Theorem U); unbounded support size at one fixed parameter (Theorem A) | q_n* >= q_8*(limit) about 1.7107 (n = 8, above the family bound 1709/1001), 17691/10010 (n = 9), 18418/10001 (n = 10), 25000/10001 (n = 47) | proved, internally reviewed (not peer-reviewed); instances certified (n >= 8 not hand-verified) |

So for n <= 7 the conjecture holds at every legal parameter with q >= rho (the plastic number, real root of
q^3 - q - 1, about 1.3247) and fails for every q in (1, rho); at n = 8 it fails at legal points with q
arbitrarily close to every value in (1, q_8*(limit)) (q_8*(limit) about 1.7107): for q < rho by the n <= 7
map, above that by the stable interval of `&GMx|\Aw^qz]w`, whose lower q end is below 5/4; whether it holds
above q_8*(limit) is open until the finite-eps map is done. Geneson's safe region q >= n - 2 is replaced
by a constant at n <= 7, while the families show the true threshold is unbounded in n (Theorem U); the only
upper bound is Geneson's n - 2, and the growth rate is open.

## 2. Structural facts

* **Stability is not monotone in the limit direction.** "Stable at finite eps => limit-stable at the
  same beta" is false (exact certificates), so the eps -> 0 classification is not the whole story:
  40 of the 124 seven-node regions bulge past beta_c and are bounded in t, exactly when R* > 0 and the
  fixed point is positive at beta_c (Lemma AS; exact R* for all 124), and one graph exists only at
  finite eps (spectral abscissa exactly -1/2, a moving eigenvalue pair on the boundary line). This is
  why the n = 8 limit map is explicitly not a finite-eps map.
* **Ray reduction (Lemma R).** Stability at any legal parameter forces stability on the legality
  curve t(beta - 2) = 1, so "stable somewhere" is an exact one-dimensional test; this made the
  complete finite-eps map at n = 7 a 16-CPU-hour computation instead of thousands.
* **Necessary conditions valid at every parameter** (all machine-checked in Lean): a dominated node
  forbids a positive full-support fixed point (Lemma D); stability forces det K > 0 and makes
  positivity t-free (Lemma A.2); a stable full support has more bidirected than non-adjacent pairs,
  hence more than half of all arcs (e_2 bound; at n = 8 this cuts the space to 44.7% of all digraphs, and
  the bound is attained by two 29-arc counterexamples); every coefficient of det(mu I + K) is positive
  (Stodola); two non-adjacent in-twins force the eigenvalue delta (never stable).
* **Closed-twin blow-ups** preserve the full-support fixed-point data of the weighted quotient (Lemmas
  W, X), which is how every infinite family is built; they do not preserve the region (46 of the 61
  n = 8 lifts have a different interval from the parent, one disjoint), most blow-ups of counterexamples
  are not limit counterexamples (797 of 858 at n = 8; their finite-eps status is open for those neither
  dominated nor never-positive), and blow-ups of non-counterexamples can be (125 at n = 8).
* **What the n = 8 counterexamples share** (`notes/REPORT_n8.md` section 5, `notes/n8_structure.md`):
  11 to 27 of 56 arcs missing (never near-cliques); largest bidirected clique at most n - 3 = 5 (61 attain
  it; the bound n - 3 is attained at n = 6, 7, 8 and in 38 proved families for every class size); all
  strongly connected, Hamiltonian, with an induced oriented 3-cycle. These features are shared by the
  338,505 strict survivors that are not counterexamples, so they do not distinguish counterexamples (whether
  they follow from the necessary conditions is open: S1+ is checked for n <= 9 only, and U1/U2 are not
  proved). **Conjecture S** (clique number <= n - 3; Hamiltonian; an induced oriented
  3-cycle) is proposed; S1 (for limit counterexamples) is verified exhaustively for every n <= 9 on the only
  class where it can fail.
  The q record at n = 8 comes from a sink-free core; its non-edge spectrum (minimum real part about -0.3689)
  is closer to -1/2 than the k -> infinity limit of Theorem Q (about -0.3516), so the sink-class family
  overtakes it, but only for very large class size (k >= 286 for q close to 9/5); two graphs have the rational critical value q = 5/3.
* **Refuted on the way:** Conjecture P (positivity of the resolvent vector forces an eigenvalue at
  Re <= -1/2), "every limit polynomial has a root on Re = -1/2", and heredity (5,902 of the 9,133 n = 8
  counterexamples have no vertex deletion that even passes the n = 7 prefilter).

## 3. Verification chain

1. Validation gate: Geneson's counterexample reproduced (permanent test); the published n = 5
   classification reproduced with diff zero; the enumerations reproduce OEIS A000273 through n = 7, and
   at n = 8 the restricted space is counted exactly by Burnside (801,558,826,511) and reproduced by the
   prefilter totals, job by job and shard by shard.
2. Exact engine: Hurwitz determinants on exact characteristic polynomials, three-valued (a vanishing
   minor is undecided, never a verdict); cylindrical decomposition on a proved boundary set; every root
   isolation replayed; pointwise cross-checks at every cell. 0 escalations at n = 5, 6, 7, 8.
3. Independent audits (agents without access to the project's conclusions): engine and n = 6 map;
   the n = 7 prefilter, resolver, curve test, edge bound; every theorem (Theorems F, E, Q, Q2, Y, the
   catalogue, the q* tool) reviewed internally (not external peer review), all holding with high confidence; the two n = 8 enumeration swaps
   (K8 by bidirected subgraph, three graphs by orbits) each refereed by four agents, with per-B ledgers
   and a content check on K8 - e (143,266 canonical strict lines identical under both schemes).
4. Human sign-off: n = 6 (2026-10-14) and all 124 n = 7 graphs (2026-10-15; the user recorded this hand check
   as not fully confident and conditioned sign-off on the machine re-verification), the latter additionally
   confirmed by four separately written exact verifiers at 742 points (0 undecided). **At n = 8 the user
   declined hand-verification (2026-10-28)**; in its place two independent machine campaigns: four
   verifiers written from the definitions (9,133/9,133 at 45,665 + 18,282 + 45,974 + 13,388 points, 0
   refuted, 0 undecided), then six further lenses (completeness via the 858 n = 7 blow-ups and the
   29,116-switch Lemma J closure, stability for all large t, the never side at 874,996 points, the boundary
   rule on 36.9 M regenerated digraphs, a float dynamics lens, an end-to-end referee of the enumeration
   chain): 0 fatal, 0 major, 1 minor documentation defect. The n = 8 result is therefore **verified by
   machine, not by hand**, and every statement of it carries the limit-only scope.
5. Lean 4 + Mathlib: 87 theorems, no `sorry`, standard axioms only (`notes/lean/`); the e_2 bound and
   Lemma D among them are the two prunes of the n = 8 chain. Lemma I (the spectral limit) and the
   bidirected-subgraph decomposition lemma are not formalised.
6. Nondegeneracy certified at a witness point of every n <= 7 hit (not run on the n = 8 list);
   literature checked on 2026-10-28 (no overlapping work found in a limited arXiv search, index to 27 Oct 2026;
   nothing on n >= 7 or n = 8 was found; an unrefereed
   preprint independently confirms q_6* and the n <= 5 exclusion).

## 4. Open and pinned

* **Done since the first version:** the n = 8 limit map, run locally at low priority on this four-core
  machine (2026-10-15 to 2026-10-27, eleven container restarts: four before the per-job driver, seven during it;
  resumable per-job driver), not on rented
  compute; the full account is `notes/REPORT_n8.md`.
* **Open: the finite-eps map at n = 8.** Route B (regenerate the 417,603,634 boundary lines, curve test
  on them plus the 7,078 `never_full` lines): about 2,300-2,400 CPU-hours (regeneration plus curve test) extrapolated from the measured n = 7
  and n = 8 rates, which lists the finite-eps-only counterexamples among those lines (15 are already known,
  E_1 among them); a complete finite-eps map additionally needs a proof that alpha > -1/2 excludes stability
  at every legal point; Route A (the full
  n = 7 pipeline on the restricted space): on the order of 40,000 CPU-hours if the n = 7 positivity ratio
  holds, i.e. rented compute (`notes/REPORT_n8.md` section 8).
* **Out of reach: n = 9 by enumeration.** The e_2-restricted space at n = 9 has 5,903,691,602,470,371
  digraphs up to isomorphism (`notes/n9_feasibility_count.json`, exact Burnside, validated at n <= 8),
  7,365 times the n = 8 space, at least about 4 million CPU-hours at the n = 8 per-digraph rate (an optimistic
  extrapolation). Feasible at n = 9: the
  598,540,288 one-node extensions of the 9,133; the S3 falsification class; S1 (limit sense) is already verified there.
* **Question A, settled (Theorem A, internally reviewed, with an independent second proof):** yes. At the single
  legal point eps = 1/21, delta = 11/210 (where Geneson's clique-class family is stable only up to 46 nodes) the
  family R_k built on the 7-node counterexample `&F^U|m^QTq?` is a stable positive non-clique full
  support for every k >= 0 (n = 7 + k); every legal point with beta < 4.9607 carries such a family,
  and conditionally on the Theorem U lemmas every legal point does.
* **Question B, settled (Theorem U, internally reviewed):** q_n* -> infinity. For every N there is n0(N) with
  q_n* >= N - 1 for all n >= n0(N); certified thresholds q_n* >= 2, 3, 4, 9/2, 23/5 for n >= 1812, 219,
  657, 2475, 4909. The growth rate (Geneson's n - 2 versus these super-exponential n0) is open; the
  complete maps give q_6* 1.1664, q_7* 1.3247, q_8*(limit) 1.7107 (approximations).
* **Conjecture S** (`notes/n8_structure.md` section 3): prove S1 via its finite reduction; test S2, S3
  at n = 9 on the one-node extensions.
* Formalisation gaps: the converse of Lemma R (openness of stability), the spectral limit lemma
  (Lemma I) and the decomposition lemma of `notes/k8_decomposition.md` are not in Lean.
* Beyond n = 8, the structural route (families, Lemma R, the blow-up calculus, Conjecture S) is the way
  forward; the n = 7 and n = 8 data say most counterexamples at those sizes are new rather than inherited
  (primitive share 3/3, 102/123, 6,720/9,133).

## 5. Where everything is

`certificates/INDEX.md` lists every certificate and the script that regenerates it (n <= 7);
`ANOMALIES.md` is the append-only log of every surprise and correction (read it before trusting any
number in a note that it supersedes); `notes/REPORT_n6.md`, `notes/REPORT_n7.md` and `notes/REPORT_n8.md`
are the full reports; `n7/REPRODUCE_n7_finite.md` and `n7/REPRODUCE_n8_limit.md` are the runbooks (the
n = 8 run as executed is documented in ANOMALIES.md 2026-10-16 to 2026-10-28 and in `notes/k8_decomposition.md`);
`n7/handcheck.py` prints one candidate at one parameter point for checking by hand. n = 8 files:
`n8/n8_limit_HITS.jsonl` (the 9,133 with exact intervals and witnesses), `n8/n8_limit_results.jsonl`
(all 347,639 strict survivors), `n8/n8_limit_tally.json`, `n8/n8_accounting.json`,
`certificates/n8_limit_summary.json`, `certificates/n8_summary.json`, `certificates/n8_counterexamples_classes.json`,
`certificates/n8_candidates_legal_witnesses.json` (use `witness_in_first_stable_cell` for a hand check),
`notes/verify_n8/SYNTHESIS.json` and `notes/verify_n8_2/SYNTHESIS.json` (the two machine campaigns),
`notes/verify_n8_2/handcheck_kit/HANDCHECK.md` (six fully worked candidates), `notes/n8_structure.md` with
`notes/n8_structure_tools/` (structure analysis), `notes/literature_position_2026-10-28.md`,
`notes/n9_feasibility_count.json`.
