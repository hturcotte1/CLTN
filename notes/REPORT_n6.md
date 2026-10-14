# The target-free clique conjecture at n = 6: a certified map of the gap

Status: COMPLETE for the full-support problem at n = 6 (and hence, by
Theorem F, for the conjecture at n = 6). Every number below is backed by a
file in `certificates/` and every lemma by a proof in `notes/`. The two new
counterexamples were independently verified by adversarial agent workflows
and were hand-verified and signed off by the project owner on 2026-10-14
(CLAUDE.md §3 gate closed; `human_sign_off` field in both certificates). A
replay of the limit sweep with the final code and an adversarial audit of the
whole chain (section 7) are complete.

Definitions used throughout: "stable" means linearly (asymptotically) stable:
every eigenvalue of `-I + W_sigma` has strictly negative real part, which for a
strictly positive full-support fixed point is asymptotic stability of the TLN
(the vector field is linear there); points where `det(I - W_sigma) = 0` or a
Hurwitz minor vanishes are non-hyperbolic and are not counted as stable.
"Legal" means `q = delta(1-eps)/eps > 1`, i.e. `0 < eps < delta/(1+delta)`;
the boundary `q = 1` is outside the CTLN definition. `theta` only rescales the
fixed point (Lemma G). The extension to proper supports through Theorem F uses
the strict off-support inequality (nondegeneracy in the sense of
Curto-Geneson-Morrison and Geneson).

## 1. Validation gates (both reproduce)

* **Stage 0, Geneson's counterexample** (arXiv:2607.21396, Thm 3.2). The paper
  gives the graph as a receiver-first non-edge matrix; our adjacency is its
  complement and transpose; the weight convention is identical to ours. The
  printed integer certificate `C v = 179298949 * 1` reproduces exactly, the
  verifier's fixed point is proportional to the paper's `v`, and the full
  support is a stable non-clique fixed point at `eps = 1/100, delta = 29/2500`
  (and `1/1000`) while being an unstable positive fixed point at `eps = 1/10,
  1/20, 1/50`. Test block 8; `certificates/stage0_geneson_2607.21396.json`.
* **Stage 1, the published n = 5 classification** (`ccurto/n5-graphs-package`).
  The package stores the transpose of our adjacency (`graph2net.m` has no
  transpose). With `A = sA^T`, every recorded fixed-point support matches the
  exact verifier on all 9,608 graphs at the package defaults `eps = 1/4,
  delta = 1/2`: diff 0 (24,442 supports, 14,488 stable, 0 violations, 0
  boundary cases). The package records no stability data. Test block 9 is a
  sampled version that also fails under the wrong convention (27 of 30).
  `certificates/stage1_n5_diff.json`.

## 2. Reformulation (proved in notes/lemmas.md, lemma_F_reduction.md, lemma_I_spectral.md)

With `beta = 1 + delta/eps`, `t = (1-eps)/eps` (legal iff `t(beta-2) > 1`,
i.e. `q > 1`) and `N` the receiver-first non-edge matrix of the support,
`I - W_sigma = eps (I + beta N + t J)`. Because `J` has rank one:

* the fixed point is positive iff all entries of `a = adj(I + beta N) 1` are
  nonzero with a common sign (Lemma H), and a positive AND stable fixed point
  forces all `a_i > 0` because `det K > 0` for a Hurwitz matrix (Lemma A.2)
  -- conditions in `beta` alone;
* at fixed `t` every verdict change happens at a real root of `a_i`, of
  `D = det(I + beta N + tJ)` (linear in `t`), of the off-support numerators
  `o_k`, or of the Hurwitz minor `Delta_{m-1}` (Lemma B, via Orlando);
* as `eps -> 0` at fixed `beta` the full support is stable iff the spectral
  abscissa of the adjacency compressed to `1^perp` is below `-(1 - 1/beta)`
  (Lemma I.3); the limit-stable set is an interval `(2, beta_c)`, so limit
  counterexamples need more than half of all ordered pairs to be edges (I.5);
* **Theorem F**: the conjecture for all n-node graphs at a parameter point is
  equivalent to "no non-complete k-node graph (k <= n) has a positive stable
  full-support fixed point" there; proper supports reduce to the induced
  subgraph plus the clique rule. Hence the n = 6 problem is exactly the
  full-support question for 5- and 6-node graphs.
* A dominated node (Curto 2510.05098 Def 1.1) forbids any full-support fixed
  point (Lemma D; proof needs uniform theta but not eps < 1), and removing it
  preserves stable supports (Lemma E). `graphs.dominates` was unsound and is
  replaced by the verbatim definition.

## 3. Engines and their validation

* `fastregion.py`: exact q-region at fixed `eps` or in the limit (integer
  evaluation + interpolation, VAS isolation returning irreducible minimal
  polynomials, cell verdicts by the validated Hurwitz routine, every sample
  cross-checked against `ctln.analyze_support`, root isolation replayed by
  Sturm counts). Agrees with the independent symbolic path `symbolic.q_region`
  on all 7,194 (graph, support, eps) cases at n = 4 and 900 at n = 5.
* `region2d.py`: cylindrical algebraic decomposition of the legal `(beta, t)`
  plane on the proven boundary set, restricted to the positivity cells, with a
  Descartes-bisection replay of every root count. Reproduces the known k <= 4
  theorem on all 237 graphs (only complete graphs ever stable) and agrees with
  the earlier Delta_5 implementation on all 706 n = 5 survivors.

## 4. n = 5: settled for all legal parameters

Of 9,608 digraphs on 5 nodes, 390 have a legal `beta` with all `a_i > 0`
(final filter; an earlier common-sign filter gave a superset of 706, which
was also certified in full); the two-parameter decomposition certifies 389 as
never stable and only `K5` as stable (0 escalations, every cell cross-checked,
regenerated with the committed engine). With Theorem F this also disposes of
every size-5 support inside every 6-node graph. This agrees with two
unrefereed GitHub claims (Aug and Oct 2026) by an independent exact method.
`certificates/region2d_n5_shard0of1_tally.json`, `certificates/n5_summary.json`.

## 5. n = 6 full support

* **Limit sweep** (`eps -> 0`), all 1,540,944 digraphs, final code (replay):
  1,512,594 graphs have no legal `beta` with all `a_i > 0`, hence (Lemma A.2)
  no stable positive full-support fixed point at any parameter -- they may
  still have positive but unstable fixed points (the empty graph does);
  28,350 strict survivors; exactly four graphs are limit-stable: `K6` and
  three non-clique graphs; 0 escalations. The first run of the sweep used a
  weaker (common-sign) filter on part of the index range and produced a
  superset of 34,071 survivors; the two-parameter certification was run on
  that superset. 204,378 graphs are domination-free and every survivor is
  among them.
* **Slices on the survivors** (0 escalations, 34,070 pointwise cross-checks):
  `eps = 1/1000, 1/100`: `K6` + the three counterexamples; `eps = 1/10`: `K6`
  + Geneson's class only; `eps = 1/4, 1/2`: `K6` only.
* **Two-parameter certification** of all 34,071 superset survivors (cylindrical
  algebraic decomposition of the legal `(beta, t)` plane, every cell decided
  exactly and cross-checked pointwise, every root count replayed): 34,067
  never stable anywhere, 4 stable somewhere (`K6` and the three graphs below),
  0 escalations, 0 errors, 65,000 CPU-seconds. Consistent with the limit
  sweep and all five slices; no survivor is missing. `certificates/n6_summary.json`,
  `notes/n6_map.md`.

**Theorem (computer-certified; human sign-off on the two new certificates recorded 2026-10-14).**
Among the 1,540,944 digraphs on 6 nodes, exactly three non-complete graphs
have a positive, linearly stable full-support fixed point at some legal
parameter: Geneson's graph, its bidirected-twin-switch partner (index
1,424,629), and the 18-edge graph of index 1,387,528. Every other 6-node
graph is certified never stable at every legal `(eps, delta)`. Combined with
the n = 5 result and Theorem F: a 6-node CTLN violates the target-free clique
conjecture at a legal parameter point if and only if it contains one of these
three graphs (or, for the 5-node part, nothing) as the full vertex set -- i.e.
the conjecture fails at n = 6 only for full supports, and only for these two
twin-equivalence classes.

### The counterexamples (full support, n = 6)

| graph | edges | limit-stable q | status |
|---|---|---|---|
| Geneson's (idx 1000255) | 19 | `1 < q < q_c`, `q_c` the sextic root ~1.16636 | known; stable for ALL legal `t` when `beta < beta_1` (degree-14 algebraic, in (17/8, 15/7)): counterexamples exist arbitrarily close to `q = 1` |
| idx 1424629 (new graph) | 19 | identical to Geneson's | Geneson's graph with one out-edge moved between bidirected in-twins (Lemma J.4: identical fixed points for every support and parameter); verified by a 7-agent workflow |
| idx 1387528 (new graph) | 18 | `1 < q < q_c'`, `q_c' = beta_1' - 1` the root of `q^5+2q^4+q^3-2q^2-2q-1` in `(1, 13/12)`, `q_c' ~ 1.0830`; at `eps = 1/101` the window is `1 < q <~ 1.0723` | genuinely new, singleton equivalence class, no twins; verified by a 7-agent workflow and by an independent from-scratch symbolic re-derivation of its whole region; witness `eps = 1/101, delta = 27/2626, q = 1350/1313` |

(Reading the cell data: a stable cell is a beta-interval together with the
stable t-range on the slice at the sample beta; `t_lo = 1/(beta_sample - 2)`
means the whole legal t-range over that beta-cell is stable. The tuples are
sectors of the decomposition, not boxes.)

Modulo bidirected twin switching there are therefore exactly two
counterexample graphs at n = 6 with a stable non-clique full support. Both
regions have the same shape: `{2 < beta < beta_c, t > max(1/(beta-2), T(beta))}`
with `T(beta) < 1/(beta-2)` on an initial interval `(2, beta_1)`, so both are
stable at EVERY legal `t` when `delta/eps` is close enough to 1. For the
18-edge graph the two thresholds coincide: stable iff `2 < beta < beta_1'`
at every legal `t`, where `beta_1'` is the root of
`beta^5 - 3beta^4 + 3beta^3 - 3beta^2 + 2beta - 1` in `(152/73, 25/12)`,
`beta_1' ~ 2.08302` (it is the positivity boundary: `a_5 = 0`).

### Structure (Lemma J)

Two graphs have identical full-support data for all parameters iff, after
relabelling, they have the same characteristic polynomial and the same in-walk
counts `(A^T)^k 1`; the 34,071 survivors form 21,800 such classes (17,307 on
the strict survivor set). Only bidirected in-twin switching is CTLN-natural
(preserves all supports); other class members can differ in their proper
supports.

## 6. What this says about the gap

For n = 6 there is no `q* > 1` below which counterexamples start: both
counterexample classes are stable at every legal point with `delta/eps` close
enough to 1, so counterexamples exist for every `q > 1` arbitrarily close to
1. The meaningful quantity is the upper reach of the counterexample region:
`q_c = 1.16635...` (the sextic root) for Geneson's class and
`q_c' = beta_1' - 1 ~ 1.0830` (the root of `q^5+2q^4+q^3-2q^2-2q-1`; at
`eps = 1/101` the window ends near 1.0723) for the other, against Geneson's
proven safe bound `q >= 4 - 1.5 eps`. The gap
`[~1.17, 4 - 1.5 eps)` at n = 6 contains NO counterexample for any graph:
the conjecture holds there for every 6-node CTLN. Every counterexample at
n = 6 needs `eps` small: none exists at `eps = 1/4` or `1/2`; at `eps = 1/10`
only Geneson's class survives (for `1 < q < ~1.03`).

## 7. Soundness audit and trust base

An adversarial audit workflow (`notes/audit_n6_workflow_full.json`: CAD
soundness review, independent from-scratch symbolic re-derivation of the
idx 1387528 region -- same set `{2 < beta < beta_1', all legal t}` and the
same minimal polynomial -- spot-checks of never verdicts, lemma referee,
completeness critic) returned "major gaps, none blocking, no verdict error".
The major items were: the human sign-off is outstanding (by design); the
never-side records are not self-contained certificates (they store the
critical values and cell counts, not the full decomposition) and
verification means re-running the committed engine; several misstatements
in an earlier draft of this report (fixed in this version); and the hit
files were not tracked in git (fixed: every HITS file, gzipped bulk records
and a sha256 manifest are now committed). A prover-written two-parameter
version of Lemma B (`notes/lemma_Bprime_2d.md`: sign-invariance on connected
components of the complement of the boundary curves, correctness of deciding
only full-dimensional cells via openness, and elementary delineability for
the projection set used) closed the proof-coverage gap. What the certificates still
trust: sympy 1.14's exact polynomial algebra (Berkowitz determinant,
factorisation, gcd, resultant, discriminant) for the projection polynomials;
the projection identities were re-verified by independent Bareiss
determinants on a sample of records (`notes/prover_scripts_B/`). The Descartes
replay is independent of sympy's VAS isolator but uses sympy's squarefree
part. Routh-Hurwitz and Orlando are cited, not re-proved.

## 8. Literature status

Scoop checks on 2026-10-12 and 10-13: nothing on arXiv maps the gap; two
unrefereed repos claim n <= 5 clean and the seed's region (consistent with our
results); Geneson's recent papers are elsewhere. J.4' is a special case of
Curto-Geneson-Morrison Lemma 3.13 / Cor 2.7; J.4 and J.2 were not found in
print. See `notes/literature_twins_check.md`, `notes/overlap_report_github_repos.md`.
