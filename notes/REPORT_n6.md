# The target-free clique conjecture at n = 6: a certified map of the gap

Status: DRAFT written while the two-parameter certification of the n = 6
positivity survivors is still running (section 5 will be completed from
`certificates/n6_summary.json`). Everything else below is final and backed by
files in `certificates/` and proofs in `notes/`. Two claimed counterexamples
have been independently verified by adversarial agent workflows but still
await the human hand-verification that CLAUDE.md §3 requires.

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
  positive (Lemmas A, A.2, H) -- a condition in `beta` alone;
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

Of 9,608 digraphs on 5 nodes, 706 admit a positive full-support fixed point
for some `beta`; the two-parameter decomposition certifies 705 as never stable
and only `K5` as stable (0 escalations, every cell cross-checked). With
Theorem F this also disposes of every size-5 support inside every 6-node
graph. This agrees with two unrefereed GitHub claims (Aug and Oct 2026) by an
independent exact method. `certificates/region2d_n5_shard0of1_tally.json`.

## 5. n = 6 full support

* **Limit sweep** (`eps -> 0`), all 1,540,944 digraphs: 1,506,873 cannot have
  a positive fixed point at any legal `beta`; 34,071 survivors; exactly four
  graphs are limit-stable: `K6` and three non-clique graphs; 0 escalations.
  204,378 graphs are domination-free and every survivor is among them.
* **Slices on the survivors** (0 escalations, 34,070 pointwise cross-checks):
  `eps = 1/1000, 1/100`: `K6` + the three counterexamples; `eps = 1/10`: `K6`
  + Geneson's class only; `eps = 1/4, 1/2`: `K6` only.
* **Two-parameter certification** of all 34,071 survivors: PENDING (fill in
  from `certificates/n6_summary.json`: totals, hits, consistency with limit
  and slices).

### The counterexamples (full support, n = 6)

| graph | edges | limit-stable q | status |
|---|---|---|---|
| Geneson's (idx 1000255) | 19 | `1 < q < q_c`, `q_c` the sextic root ~1.16636 | known; stable for ALL legal `t` when `beta < beta_1` (degree-14 algebraic, in (17/8, 15/7)): counterexamples exist arbitrarily close to `q = 1` |
| idx 1424629 (new graph) | 19 | identical to Geneson's | Geneson's graph with one out-edge moved between bidirected in-twins (Lemma J.4: identical fixed points for every support and parameter); verified by a 7-agent workflow |
| idx 1387528 (new graph) | 18 | `1 < q < root of q^5+2q^4+q^3-2q^2-2q-1 in (1, 13/12)` (~1.0723 at eps=1/101) | genuinely new, singleton equivalence class, no twins; verified by a 7-agent workflow; witness `eps = 1/101, delta = 27/2626, q = 1350/1313` |

Modulo bidirected twin switching there are therefore (so far) exactly two
counterexample graphs at n = 6 with a stable non-clique full support.

### Structure (Lemma J)

Two graphs have identical full-support data for all parameters iff, after
relabelling, they have the same characteristic polynomial and the same in-walk
counts `(A^T)^k 1`; the 34,071 survivors form 21,800 such classes (17,307 on
the strict survivor set). Only bidirected in-twin switching is CTLN-natural
(preserves all supports); other class members can differ in their proper
supports.

## 6. What this says about the gap

For n = 6 there is no `q* > 1` below which counterexamples start: Geneson's
graph is stable at every legal point with `delta/eps` close enough to 1. The
meaningful quantity is the upper reach of the counterexample region, which
for the two classes is `q_c ~ 1.166` and `~1.08`, against Geneson's proven
bound `q >= 4 - 1.5 eps`. (Final statement after section 5 is complete.)

## 7. Literature status

Scoop checks on 2026-10-12 and 10-13: nothing on arXiv maps the gap; two
unrefereed repos claim n <= 5 clean and the seed's region (consistent with our
results); Geneson's recent papers are elsewhere. J.4' is a special case of
Curto-Geneson-Morrison Lemma 3.13 / Cor 2.7; J.4 and J.2 were not found in
print. See `notes/literature_twins_check.md`, `notes/overlap_report_github_repos.md`.
