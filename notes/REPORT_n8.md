# The target-free clique conjecture at n = 8: the complete eps -> 0 limit map, and what it changes

Date: 2026-10-28. Companion to `notes/REPORT_n7.md` (the n = 7 finite-eps map, Lemma R, the infinite
families, q_7*) and `notes/REPORT_n6.md` (definitions, conventions, the n = 6 theorem). Conventions are
those of the project instructions file, section 2, throughout; everything below is exact (integers, rationals, exact root isolation,
Hurwitz determinants); decimals appear only as labelled approximations ("about", "display").

> **Status and scope.** This report covers the **eps -> 0 limit map of the full support at n = 8**: every
> digraph on 8 vertices with at least 29 arcs (801,558,826,511 of them) was read exactly once and the
> 347,639 strict spectral survivors were classified exactly, giving K8 plus **9,133 non-complete graphs
> whose full support is limit-stable**: 9,133 counterexamples to the conjecture in the limit sense. The
> list is **VERIFIED BY MACHINE** by two independent campaigns (first: four verifiers written from the
> definitions, 9,133/9,133; second: six lenses including completeness from the n = 7 side, stability for
> all large t, the never side and an end-to-end referee of the enumeration chain, 0 fatal or major
> defects) and has **no human hand check**: the user declined hand-verification on 2026-10-28 and ordered
> the second machine campaign instead (`certificates/n8_limit_summary.json`, `status`,
> `second_verification.human_decision`; ANOMALIES.md 2026-10-28). **Not done:** the finite-eps map at
> n = 8. The 417,603,634 boundary (B) lines of the prefilter (alpha(A|perp) >= -1/2, never limit-stable)
> were not stored, so finite-eps-only counterexamples of the `&FNFL~[udo?` type (alpha = -1/2 exactly,
> Theorem E) that live among them are outside this run; 15 such eight-node graphs are already known (the
> Lemma B graphs of `notes/research_wf3/theory/cube_n8_witnesses.json`, E_1 of Theorem E among them, each
> with an exact legal witness re-decided in `notes/referee_finish_overclaims_checks.json`, block B; none
> hand-verified), and the complete finite-eps list is not determined. Proper supports at n = 8 are not covered either
> (Theorem F reduces them to induced subgraphs on <= 7 nodes, mapped in `notes/REPORT_n7.md`, plus an
> off-support condition for the eighth node that was not enumerated). Every statement below carries
> these caveats (`certificates/n8_limit_summary.json`, `scope`; `certificates/n8_summary.json`, `scope`).

## 0. Headline results

Sources: `certificates/n8_summary.json` unless stated; `n8/n8_limit_tally.json`; `n8/n8_accounting.json`;
`notes/n8_structure.md`; `notes/verify_n8/SYNTHESIS.json`; `notes/verify_n8_2/SYNTHESIS.json`.

| # | result | status |
|---|---|---|
| 1 | **Complete limit map at n = 8.** Among the 801,558,826,511 digraphs on 8 vertices with >= 29 arcs (the e_2 edge bound, Lean-checked, excludes the other 991,800,366,337 of the 1,793,359,192,848 digraphs on 8 vertices at every legal parameter), exactly **9,133 non-complete graphs** have a full support that is a limit-stable positive fixed point (stable at all large t on some legal ray); K8 is the only complete one; 0 escalations, 0 overflows. Each carries an exact certificate: one open stable beta-interval with algebraic endpoints (minimal polynomial + isolating interval), four stored pointwise checks (all four stable and positive for 9,051 records; section 3, "Witness hygiene"), and a certified legal interior witness. | CERTIFIED; VERIFIED BY MACHINE (two campaigns); no hand check |
| 2 | **Scale.** 3 (n = 6), 123 (n = 7), 9,133 (n = 8) limit counterexamples. **9,072 of the 9,133 are not closed-twin lifts** of n = 7 counterexamples (61 are); 6,720 contain no n = 6 or n = 7 counterexample as an induced subgraph; for 5,902 no single-vertex deletion even passes the n = 7 prefilter. The n = 8 list is overwhelmingly new, not inherited. | CERTIFIED (lifts recomputed twice) |
| 3 | **Threshold in the limit.** `q_8*(limit)`, the supremum of q over the 9,133 stable intervals, is the root of the degree-15 polynomial `13109 q^15 - 124967 q^14 + ... - 9045013` in `[459224485/268435456, 2755346911/1610612736]` (about 1.71074), a supremum approached by `&GMx|\Aw^qz]w` and `&GVrm]A{^pZ]w` (one bidirected in-twin switch apart) and **not attained** (the Hurwitz determinant vanishes there). It exceeds q_7* = the plastic number (about 1.3247) and the family lower bound 1709/1001 (about 1.7073) of REPORT_n7. All 9,133 upper ends are below 2; 16 exceed 3/2; the smallest is about 1.000122. The finite-eps q_8* is not determined by this run (it is >= q_8*(limit)). | CERTIFIED (ordering certified exactly) |
| 4 | **The record comes from a sink-free core.** The two extremal graphs have no node sending to all others; their q is limited by the spectrum of their own non-edge digraph (minimum real part about -0.3689, display), which lies to the left of (closer to -1/2 than) the k -> infinity limit of Theorem Q (about -0.3516, display, i.e. q -> 1.844); the sink-class families overtake the n = 8 record only for large k (q close to 9/5 needs k >= 286, n >= 293), and at n = 8 itself the record is attained with no sink class at all (exact comparison: `notes/referee_finish_overclaims_checks.json`, block C). Two other graphs have the **rational** critical value beta_c = 8/3, q = 5/3 exactly, from the factor `4z^2 + 3z + 1`. | CERTIFIED (`notes/n8_structure.md` 1.8) |
| 5 | **Structure.** Every counterexample misses 11 to 27 of the 56 arcs; the largest bidirected clique has at most 5 = n - 3 vertices (61 attain it); every one is strongly connected, Hamiltonian, and contains an induced oriented 3-cycle. These universal features do not single out counterexamples: the 338,505 strict survivors that are not counterexamples share them (strong connectivity fails for 18). **Conjecture S** (clique number <= n - 3; strongly connected; an induced oriented 3-cycle) is proposed; its part S1 is verified, in the limit sense, for every n <= 9 by an exhaustive test of the bidirected-K_{n-2} class (682,176 configurations at n = 9, only K9 survives) and holds in all 631 proved infinite families for every class size. | DATA + CONJECTURE |
| 6 | **Completeness cross-checks from the proved side.** All 626 n = 8 members of the clique-class catalogue (593 distinct graphs) are in the list; Theorem F's F_2 is in the list with its interval ending at the plastic number as proved; of the 858 distinct closed-twin blow-ups of the 123 n = 7 counterexamples, 61 are in the list and the other 797 are not limit counterexamples (their finite-eps status is open where they are neither dominated nor never-positive); the Lemma J in-twin-switch closure of the list (29,116 switches) stays inside the list with identical intervals. | CERTIFIED |
| 7 | **The enumeration chain** reads the restricted space exactly once: the Burnside count 801,558,826,511 is reproduced by the prefilter totals (accounting ok, 9 shards, 3,226 jobs, 201 split graphs checked), including two swaps of enumeration scheme forced by high-symmetry graphs (K8 by bidirected subgraph, three graphs by orbits), each with a proved lemma, a per-B ledger and a content check (K8 - e: 143,266 canonical strict lines identical under both schemes). | CERTIFIED; chain refereed |
| 8 | As of 2026-10-28 a limited search (arXiv, index to 27 Oct 2026; no other preprint servers) found no posted work on n = 8, n = 7, finite-eps-only counterexamples or an improvement of Geneson's threshold (section 7). | searched, limits in section 7 |

## 1. Setting

The notation is that of REPORT_n7 section 1: `I - W = eps (C + tJ)`, `C = I + beta N`, `N[i][j] = 1` iff
`j -/-> i`, `t = (1 - eps)/eps`, `beta = 1 + delta/eps`, `q = delta(1 - eps)/eps = (beta - 1) t/(t + 1)`,
legal iff `q > 1` iff `t(beta - 2) > 1`. The full support is a stable fixed point iff `K = C + tJ` is positive
stable (Hurwitz determinants, three-valued) and `adj(K) 1 = a(beta) > 0` (t-free, Lemma A.2). In the
**eps -> 0 limit** (Lemma I; `notes/lemma_I_spectral.md`): the full support is limit-stable at beta iff the
spectral abscissa `alpha(A|perp)` of the adjacency compressed to the hyperplane orthogonal to the all-ones
vector is `< -(1 - 1/beta)`; the limit-stable set is `(2, beta_c)`, `beta_c = 1/(1 + alpha)`, nonempty iff
`alpha < -1/2`; and limit stability at beta gives stability for all large t on the ray through beta
(Lemma I.3/I.4), so a limit counterexample is a counterexample at legal parameters with eps -> 0 and
q -> beta - 1. In the limit, `q = beta - 1`.

**The restricted space.** The e_2 bound (REPORT_n7 section 5; `notes/lean/CtlnCert/CtlnCert/EdgeBound.lean`)
says a stable full support at any legal parameter has more bidirected than non-adjacent pairs, hence at
least `m(m - 1)/2 + 1 = 29` arcs at m = 8, hence an underlying simple graph with 15 to 28 edges. So the
whole question at n = 8 (limit or finite) lives on the digraphs with >= 29 arcs: 5,350 underlying graphs and
**801,558,826,511** digraphs up to isomorphism, 44.7% of A000273(8) = 1,793,359,192,848
(`n8/n8_restricted_counts.json`; exact Burnside count over Aut(G) per underlying graph, validated with no
arc bound against A000273(6) and A000273(7); ANOMALIES.md 2026-10-16).

**What the prefilter decides.** `n7/filter8.c` (the overflow-guarded n = 8 version of the n = 7 filter,
byte-identical to that filter on all 1,540,944 n = 6 digraphs and on the 638,042 n = 7 spectral survivors;
ANOMALIES.md 2026-10-15) reads each digraph once and (i) drops it if a node is dominated (Lemma D,
Lean-checked: a dominated node forces a zero coordinate at every fixed point), (ii) runs the exact Hurwitz
test of the shifted compressed characteristic polynomial and prints **S** (strict: `alpha < -1/2`, every
minor positive), **B** (a zero coefficient or minor met before any negative one: `alpha >= -1/2`), or
nothing (a negative one met first: `alpha >= -1/2` is what is proved, any non-positive minor giving a root
with Re >= 0; in the 300 sampled non-printed lines alpha > -1/2 strictly,
`notes/verify_n8_2/never_side/check_ab_summary.json`); an integer overflow would print **O** (none occurred).
Only `alpha < -1/2` can be limit-stable (Lemma I.3(c2)), so B and dropped lines are non-survivors of the
limit map, and this exclusion is sound regardless of the B/dropped split (`notes/verify_n8_2/SYNTHESIS.json`,
`never_side.own_rederivation`). B lines were discarded at the shell level and not stored; they are where
the finite-eps-only graphs of Theorem E (alpha = -1/2 exactly) sit, which is why the finite-eps map is out
of scope (section 8). The second campaign corrected an earlier description of the B lines as
"alpha = -1/2 exactly": 230 of 300 sampled B lines have alpha > -1/2 strictly
(`notes/verify_n8_2/never_side/check_ab_summary.json`, `sample_B`; `notes/verify_n8_2/SYNTHESIS.json`), and
the chain referee's separate cross-check found a later negative minor (alpha > -1/2) in 23 of 122 B lines
(ANOMALIES.md 2026-10-28, "Corrections from the chain referee").

**What the limit stage certifies.** For every S line, `fastregion.region(A, range(8), INF, replay=True)`
(`n8/run_n8_limit.py`) computes the exact critical set in beta (roots of the positivity polynomials
`a_i(beta)`, of `det K` and of the Hurwitz determinant Delta of the limit characteristic polynomial), isolates
the roots, and decides each open cell by an exact sample with replay-checked isolation. A record with verdict
`stable_somewhere` stores its stable cells as open intervals with endpoints given as minimal polynomial plus
isolating interval; the endpoint itself is never certified stable (a vanishing minor is three-valued
undecided). `never_positivity` means no beta > 2 with all `a_i > 0`; `never_full` means positivity holds
somewhere but never together with stability.

## 2. The accounting chain

Every one of the 801,558,826,511 digraphs of the restricted space is dropped by a proved necessary
condition valid for the limit map, or carries an exact limit certificate. Sources:
`n8/n8_accounting.json`, `n8/n8_limit_tally.json`, `certificates/n8_limit_summary.json`,
`certificates/n8_summary.json` (`chain`, `counts`).

| stage | condition (proved where) | dropped | remaining |
|---|---|---|---|
| 0 | all digraphs on 8 vertices, A000273(8), for reference | | 1,793,359,192,848 |
| 0' | e_2 edge bound: >= 29 arcs, underlying graph with 15..28 edges (`edge_bound_ctln` in `notes/lean/CtlnCert/CtlnCert/EdgeBound.lean`, Lean-checked; hypotheses: legal point, m >= 2, full support Hurwitz stable) | 991,800,366,337 | **801,558,826,511** on 5,350 underlying graphs (exact Burnside; read exactly once) |
| 1 | dominated node (Curto Def 1.1; Lemma D, Lean-checked) | 687,860,834,577 | 113,697,991,934 domination-free |
| 2 | spectral test `alpha(A|perp) < -1/2` (Lemma I.3), exact Hurwitz in `filter8` | 113,280,040,661 printed nothing | 417,951,273 spectral survivors = 417,603,634 **B** (alpha >= -1/2, never limit-stable, not stored) + 347,639 **S** (strict); 0 overflow |
| 3 | exact limit stage on every S line (`fastregion.region`, replay) | 331,427 `never_positivity` + 7,078 `never_full` | 9,134 `stable_somewhere` = **K8 + 9,133**; 0 escalations |

The identity `strict == never_positivity + never_full + stable_somewhere + escalate` holds
(`certificates/n8_summary.json`, `counts.accounting_identity_holds`), and the second campaign recounted
every verdict from `n8/n8_limit_results.jsonl` (338,505 never = 331,427 + 7,078; 9,134 stable; HITS equal to
the non-complete stable records, 9,133/9,133; `notes/verify_n8_2/SYNTHESIS.json`, `result_under_test`).

**Enumeration schemes and the two swaps.** The space was enumerated as `geng 8 15:28 | directg -e29:56`
in nine geng shards (`n8/parts`, shards 0..8, 5,129,155,711 digraphs) plus one job per underlying graph or
per `directg -s p/P` part for the remaining graphs (`n8/jobs`), with the exact Burnside count per graph as
the acceptance test of every unsplit job and of every split graph's part sum (`n8/check_accounting.py`:
ok, 0 problems, 9 shards, 3,226 jobs, 201 split graphs checked). Two graphs classes could not be enumerated
by `directg` in acceptable time because its isomorph rejection slows with the automorphism group
(measured 0.051 M/s at |Aut| = 192 down to 0.006 M/s at 1,440; ANOMALIES.md 2026-10-27):

| swap | what | lemma | ledger | content check |
|---|---|---|---|---|
| 1 | **K8** (575,009,339 digraphs, |Aut| = 40,320): 348 `directg` parts (estimated 3-8 days on four cores) replaced by **40 batches enumerated by bidirected subgraph B** (`n8/k8decomp`): for each of the 12,345 graphs B on 8 vertices with 1..28 edges, `directg -o` on the complement of B, then `addarcs B`, then the unchanged `filter8` | `notes/k8_decomposition.md`: isomorphism classes of digraphs with underlying graph K_n <-> pairs (B up to isomorphism, orientation of K_n - B up to Aut(B)); a re-enumeration of exactly the same classes, not a prune; four adversarial referees SOUND (`notes/k8_decomposition_referees.json`) | `n8/check_k8decomp.py`: 40 batches, 12,345 graphs B, grand sum 575,009,339 = exact 575,009,339, OK (`certificates/n8_summary.json`, `chain.swaps[0].ledger`); independent count OEIS A001174(8) - A000568(8) = 575,016,219 - 6,880 | validated at n = 4..7 by canonical-form multiset equality with `directg` (38 / 570 / 21,424 / 2,141,832 classes) and identical filter output; `addarcs` unit-tested on 3,245 cases; filter relabelling invariance tested on 4,907,680 n = 8 digraphs |
| 2 | **K8 - K3, K8 - K_{1,5}, K8 - (K2 + K_{1,5})**: 99 `directg` parts replaced by **14 orbit-enumerator parts** (`n8/orbitenum`: Aut(G) by brute force, least-element representatives of B under Aut(G) and of the orientation of G - B under the stabiliser of B) | relative form of the same lemma (`notes/k8_decomposition.md`, "Relative form"); four adversarial referees SOUND | `n8/check_orbit.py`: K8 - K3 8 parts, 63,271 representatives, 1,179,317,832 = exact; K8 - K_{1,5} 4 parts, 50,139 representatives, 330,741,720 = exact; K8 - (K2 + K_{1,5}) 2 parts, 24,369 representatives, 84,866,528 = exact (`chain.swaps[1].ledger`) | eight graphs at n = 5..7 (up to 10,840,382 classes) and n = 8 graphs (128,067 classes) set-equal to `directg` with identical filter output; **post hoc on K8 - e** (5,339,129,738 digraphs, the richest graph): the 16 orbit parts give exactly the 143,266 canonical strict-survivor lines of the 321 `directg` parts the run used (`n8/orbitenum/k8e_check`; ANOMALIES.md 2026-10-28) |

Both swaps are recorded in ANOMALIES.md 2026-10-27 ("K8 enumerated by its bidirected subgraph", "Three
more high-symmetry graphs enumerated by orbits"). The chain referee of the second campaign re-derived the
edge bound, the Burnside counts (30 graphs, two shards, 20 dense graphs), coverage exactly once, filter
agreement on 3,391 digraphs and the accounting, and found no gap; its one "misattribution" finding
(38,215 lines) was a pre-canonicalisation artefact, resolved to 0 after canonical labelling
(`notes/verify_n8_2/SYNTHESIS.json`, `chain_referee`). All 347,639 strict lines come from 17 high-symmetry
underlying graphs (K8: 166,042; K8 - e: 143,266; K8 - P3 (graph6 `GV~~~{`): 19,701; K8 - 2K2 (`G]~~~{`): 16,123; ...;
recount of every strict line's underlying graph with `labelg` in `notes/referee_finish_numbers.json`; the label
"K8 - 2K2" for 19,701 in ANOMALIES.md 2026-10-27/28 is a mislabel of `GV~~~{`, which ANOMALIES.md 2026-10-22
identifies as K8 minus P3), none from the 135 packed jobs
(55,077,003,264 digraphs) or from the nine geng shards (ANOMALIES.md 2026-10-27/28), so the per-B ledgers
of the swapped schemes are load-bearing, not a formality.

**The sharded limit stage.** The single-threaded limit stage was killed by a container restart at 59,610 of
347,639 lines and finished as four shards of the unmodified `n8/run_n8_limit.py` over four disjoint
symlinked subsets of the job files, merged by `n8/limit_sharding/merge_finalize.py` (disjoint sets whose
union is exactly the S/O line set; tally recomputed from the merged records equals the sum of the shard
tallies; every one of the 64,196 complete records of the interrupted single run equals the shard record
for the same graph). Shard tallies (`certificates/n8_summary.json`, `chain.limit_stage.shard_tallies`;
`n8/limit_sharding/shard{0..3}_tally.json`):

| shard | S | never_positivity | never_full | stable_somewhere | complete | seconds |
|---|---|---|---|---|---|---|
| 0 | 86,910 | 82,828 | 1,770 | 2,312 | 0 | 2,733.6 |
| 1 | 86,910 | 82,867 | 1,760 | 2,283 | 1 (K8) | 2,675.9 |
| 2 | 86,910 | 82,873 | 1,727 | 2,310 | 0 | 2,697.9 |
| 3 | 86,909 | 82,859 | 1,821 | 2,229 | 0 | 2,669.2 |
| sum | 347,639 | 331,427 | 7,078 | 9,134 | 1 | 10,776.6 |

**Costs and timeline.** The prefilter ran locally on this four-core machine at low priority from
2026-10-15 to 2026-10-27 (ANOMALIES.md entries 2026-10-15 to 2026-10-27), through eleven container
restarts (four on 15-16 October during the geng-shard phase, seven from 19 to 27 October during the per-job
phase; ANOMALIES.md 2026-10-15, 2026-10-16, 2026-10-22, 2026-10-23, 2026-10-27) with a resumable per-job driver; measured aggregate rate about 1.5-1.7 M digraphs/s with four
pipelines (`filter8` the hot process at 0.4 M/s per pipeline; ANOMALIES.md 2026-10-16), i.e. about 5.5-6
days of clean wall time for the restricted space, the rest being restarts and the `directg` tail on
high-symmetry graphs that the two swaps removed. The limit stage took 10,776.6 CPU-seconds. The aggregation
(`n8/aggregate_n8.py`, exact, 2 workers) took 117 s.

Files: `n8/n8_limit_HITS.jsonl` (one record per counterexample: digraph6, arcs, the exact stable interval
with minimal polynomials and isolating intervals in beta and q, the vanishing critical polynomials at each
end, four pointwise witnesses), `n8/n8_limit_results.jsonl` (all 347,639 verdicts), `n8/n8_limit_tally.json`,
`n8/n8_accounting.json`, `n8/n8_prefilter_TOTAL.txt`, `n8/n8_restricted_counts.json`;
`certificates/n8_limit_summary.json` (scope, enumeration, limit stage, status, both verification blocks),
`certificates/n8_summary.json` (the full chain, counts, D-classes, q data per record, q_8*, provenance),
`certificates/n8_counterexamples_classes.json` (6,260 D-classes), `certificates/n8_candidates_legal_witnesses.json`
(one legal witness per graph). Runbook of the original plan: `n7/REPRODUCE_n8_limit.md`; the sharding by
underlying graph, the job driver and the swaps are documented in ANOMALIES.md 2026-10-16 and 2026-10-27.

## 3. The 9,133 graphs

Sources: `certificates/n8_summary.json`, `certificates/n8_limit_summary.json`, `notes/n8_structure.md`.

* **Arc counts** 29 to 45 of 56 (29: 2, 30: 27, 31: 120, 32: 269, 33: 562, 34: 1,079, 35: 1,485, 36: 1,406,
  37: 1,305, 38: 768, 39: 510, 40: 476, 41: 581, 42: 353, 43: 153, 44: 34, 45: 3), i.e. 11 to 27 missing arcs.
  **The e_2 bound is attained**: the two 29-arc graphs `&GM?WUZeEie^O` and `&GMBO]At]Xp]?` are
  near-tournaments (3 and 2 bidirected pairs, 2 and 1 non-adjacent pairs, no twins, trivial automorphism
  group, no node sending to all others) with q windows just above 1 (upper ends about 1.01309 and 1.00116,
  display). At n = 7 the minimum of 22 arcs was allowed but not attained (23 was the minimum).
  All 9,133 are domination-free, non-complete, non-clique and pairwise non-isomorphic (own canonical form
  and nauty `labelg` agree; the stored `edges` field equals the arc count in every record).
* **Non-adjacent pairs** 0: 6,037, 1: 2,895, 2: 198, 3: 3; bidirected pairs 2 to 17 of 28; largest
  bidirected clique 2: 2,202, 3: 5,847, 4: 1,023, 5: 61, never 6 (`notes/n8_structure.md` 1.1-1.3).
* **D-equivalence classes** (Lemma J key: integer characteristic polynomial of A plus the sorted multiset
  of per-vertex in-walk sequences `((A^T)^k 1)_i`, k = 1..8; Theorem J.2 of `notes/lemma_J_equivalence.md`:
  equal keys iff identical full-support region at every (eps, delta)): **6,260 classes**, sizes 1: 4,939,
  2: 824, 3: 2, 4: 374, 5: 16, 6: 3, 7: 1, 8: 84, 10: 2, 12: 1, 14: 4, 16: 8, 19: 1, 32: 1. Within every
  class the stored lower-end kind, the upper-end minimal polynomial and root, the lower-end root where
  present, and the arc count are identical (0 exceptions); the only within-class differences are 262
  members whose `vanishing` labels name a different vertex index of a positivity polynomial `a_i`, the same
  datum up to the vertex permutation that D-equivalence allows (`certificates/n8_counterexamples_classes.json`).
* **Twins.** 4,271 graphs carry at least one bidirected in-twin pair (5,202 pairs in all); the Lemma J
  in-twin-switch closure of the list (29,116 switches of these 4,271) stays inside the list with identical
  intervals (0 exceptions, 0 missed counterexamples; `lemma_J_switch_closure_second_campaign`). 186 graphs
  have closed twins (same in- and out-neighbourhoods plus the bidirected pair); 188 have a non-trivial
  automorphism group (186 from closed twins, 2 others). No graph has non-adjacent in-twins (Remark J.4',
  Lean-checked, never stable).
* **Lifts versus new.** Closed-twin blow-ups (Lemma W) of the 123 n = 7 limit counterexamples: 861 built,
  858 distinct, **61 in the list** (from 52 parents), 797 not in the list; those 797 are not limit
  counterexamples (prefilter flags: 690 dropped, 47 B, 60 never in the limit stage, i.e. 43 `never_positivity`
  and 17 `never_full`; `notes/verify_n8_2/completeness_twins/part2_summary.json`), and 0 of 4,782 sampled legal
  points inside the parents' intervals is stable and positive (a sample, not an exact decision of the whole
  interval; the finite-eps status of the B, `never_full` and alpha > -1/2 ones is open); so
  **9,072 of the 9,133 are not single closed-twin lifts**. Lemma W preserves the fixed-point data of the
  weighted quotient, not the region: 46 of the 61 listed blow-ups have a different interval from the
  parent, one disjoint (`&G^z~Vwn?zBRw`; `notes/verify_n8_2/SYNTHESIS.json`). Conversely 125 of the 186
  closed-twin graphs are blow-ups of graphs that are **not** n = 7 counterexamples (114 of them not even
  strict n = 7 survivors): the blow-up calculus is not monotone in either direction. By induced
  subgraphs: 2,323 contain an n = 7 counterexample (all 123 occur; 5 contain `&FNFL~[udo?`), 308 contain an
  n = 6 counterexample (all three occur), **6,720 contain neither**; for 5,902 no single-vertex deletion is
  even a strict n = 7 survivor (`notes/n8_structure.md` 1.6).
* **Stable intervals.** Every record has exactly one stable interval, bounded above (0 unbounded). Lower
  end: beta = 2 (q -> 1) for 8,162, an algebraic positivity root above 2 for 971 (never a spectral root:
  limit stability holds down to beta = 2 whenever positivity does, as Lemma I.3 says). Upper end: the
  spectral boundary beta_c (Delta vanishes) for 7,191, a positivity root for 1,942; the exact positivity
  sets recomputed in `notes/n8_structure.md` 1.7 agree with every stored interval (9,133/9,133, 0
  inconsistencies; the mechanism table there: 4,690 + 1,974 + 1,498 from beta = 2, 504 + 444 + 23 from a
  root). 4,298 distinct upper-end minimal polynomials and 147 distinct lower-end ones, all irreducible
  (sympy); degrees of the upper ends 1: 2, 3: 208, 4: 34, 5: 947, 6: 874, 7: 307, 10: 535, 12: 8,
  15: 4,022, 21: 2,196. The most frequent upper-end polynomial is `q^5 + 2q^4 + q^3 - 2q^2 - 2q - 1` (408
  graphs), the same quintic as the n = 6 counterexample idx 1387528 and the most frequent at n = 7 (17
  graphs there); the plastic number `q^3 - q - 1` ends 178 intervals (6 at n = 7); every stored q
  polynomial equals the beta polynomial composed with beta = 1 + q, and every isolating interval contains
  exactly one root (Sturm; `q_data.stored_minpolys_checked`).
* **Narrow intervals.** 37 intervals starting at beta = 2 are narrower than 1/1000 in beta (plus 7 narrower
  ones starting above 2; 0 undecided at width 1e-9). All 37 end at the spectral boundary just above 2
  (largest q end among them about 1.00098, display): `alpha(A|perp)` is barely below -1/2. None is a
  lift or a catalogue member, none has twins (`notes/n8_structure.md` 1.9).
* **Witness hygiene** (candidates unaffected). 6,435 records have an illegal first stored witness (t = 100,
  q <= 1), the n = 7 pattern again; the 37 narrow-interval records have all four stored witnesses outside
  the interval (producer fallback beta = lo + 1/1000) and store `stable = false` there; 16 graphs have only
  the t = 100 witness unstable (a finite-eps effect). 9,051 records have all four witnesses stable and
  positive, 9,096 have a legal stable positive stored witness, and
  `certificates/n8_candidates_legal_witnesses.json` gives a certified legal interior point
  (`witness_in_first_stable_cell`, q > 1 exact) for every one of the 9,133: use that field for any check by
  hand (`notes/verify_n8/SYNTHESIS.json`, `things_the_human_must_know`). Six candidates are worked in full
  with exact fractions in `notes/verify_n8_2/handcheck_kit/HANDCHECK.md`.

### 3b. Comparison n = 6 / 7 / 8

Sources: `notes/REPORT_n7.md` (n = 6, 7), `certificates/n7_finite_summary.json`, `certificates/INDEX.md`,
`certificates/n8_summary.json`, `notes/n8_structure.md` section 2.

| | n = 6 | n = 7 | n = 8 |
|---|---|---|---|
| digraphs enumerated | 1,540,944 (all) | 882,033,440 (all) | 801,558,826,511 (>= 29 arcs; of 1,793,359,192,848) |
| spectral survivors / strict | 28,350 strict (replay) | 638,042 / 1,889 | 417,951,273 / 347,639 |
| limit counterexamples | 3 | 123 | 9,133 |
| finite-eps-only counterexamples | 0 | 1 (`&FNFL~[udo?`) | at least 15 known (the Lemma B graphs, E_1 among them; certified witnesses, not hand-verified); complete list not determined (B lines not stored) |
| complete finite-eps map | yes | yes | **no** |
| D-equivalence classes | 2 | 91 (124 graphs) | 6,260 |
| primitive (no smaller counterexample induced) | 3 | 102 | 6,720 |
| closed-twin lifts of the previous size | - | FNFL = n = 6 hit idx 1387528 with one node doubled | 61 (of 858 blow-ups) |
| largest bidirected clique (n - 3), graphs attaining it | 3 (2) | 4 (16) | 5 (61) |
| arcs: allowed minimum / attained range | 16 / - | 22 / 23..31 | 29 / **29**..45 |
| q* (sup over counterexamples) | sextic `797q^6 - ... + 1521` in (561/481, 7/6), about 1.1664 | plastic number, root of `q^3 - q - 1`, about 1.3247, not attained | limit: degree 15, in [459224485/268435456, 2755346911/1610612736], about 1.7107, not attained; finite-eps: >= this, not determined |
| verification | hand (2026-10-14) + audit | hand (2026-10-15, recorded by the user as not fully confident) + four independent verifiers | two independent machine campaigns; **no hand check** |

## 4. Thresholds

**n <= 7 (final, REPORT_n7 section 4).** For n <= 6 the conjecture holds at every legal parameter with
q >= q_6* (about 1.1664) and fails for every q in (1, q_6*); for n <= 7 it holds for q >= rho (the plastic
number, about 1.3247) and fails for every q in (1, rho).

**n = 8, limit map** (`certificates/n8_summary.json`, `q_star_limit`, `q_data`). Define `q_8*(limit)` as the
supremum over the 9,133 records of the upper end of the stable q-interval (in the limit q = beta - 1).

* `q_8*(limit)` is the root of
  `13109 q^15 - 124967 q^14 + 538385 q^13 - 1420107 q^12 + 2727809 q^11 - 4550635 q^10 + 7298813 q^9 - 10746735 q^8 + 13594735 q^7 - 15127269 q^6 + 16652387 q^5 - 18092785 q^4 + 16098203 q^3 - 9885513 q^2 + 8301263 q - 9045013`
  in `[459224485/268435456, 2755346911/1610612736]` (width below 1e-9; about 1.71074452, display), the
  upper end of the intervals of `&GMx|\Aw^qz]w` and `&GVrm]A{^pZ]w` (same D-class, 40 arcs; stored beta
  interval [8/3, 11/4]). The **ordering is certified**: every other upper end has an exact isolating
  interval lying entirely below the champion's, and the tie is the same root (same irreducible polynomial,
  Sturm count 1 on the hull). It agrees with the first campaign's independent refinement
  (`notes/verify_n8/v_structure/q_upper_ends_refined.json`: `[21526147/12582912, 5381537/3145728]`).
* **Not attained**: each record certifies the open cell below the root only; at the root the Hurwitz
  determinant Delta of the attaining graphs vanishes (three-valued undecided), and no other record is
  certified stable at or above it.
* Exactly `q_8*(limit) > q_7* = rho` (the plastic number's isolating interval [49/37, 53/40] lies below),
  `q_8*(limit) > 1709/1001` (the n = 8 family lower bound of REPORT_n7 section 0, about 1.7073), and
  `17/10 < q_8*(limit) < 7/4`.
* **Distribution of the 9,133 upper ends** (all strictly above 1 and strictly below 2; none equal to any of
  the thresholds): above 1001/1000: 9,096; above 101/100: 8,555; above 21/20: 5,778; above 11/10: 3,190;
  above 6/5: 1,012; above 13/10: 382; above 7/5: 60; above 3/2: 16; above 8/5: 7; above 17/10: 3; above
  7/4: 0. By bin: (10001/10000, 1001/1000]: 37; (1001/1000, 101/100]: 541; (101/100, 21/20]: 2,777;
  (21/20, 11/10]: 2,588; (11/10, 6/5]: 2,178; (6/5, 13/10]: 630; (13/10, 7/5]: 322; (7/5, 3/2]: 44;
  (3/2, 8/5]: 9; (8/5, 17/10]: 4; (17/10, 7/4]: 3.
* **Smallest upper end:** `&GMBQMHtMH~]W`, root of a degree-21 polynomial in
  `[77079731/77070336, 1079116235/1078984704]` (about 1.000122, display), Delta vanishing; confirmed by an
  independent exact check (`notes/n8_structure_tools/check_min_q.json`: limit-stable at
  beta = 20001/10000, not at beta = 10001/5000). The value "about 1.000276" in
  `notes/verify_n8/SYNTHESIS.json` and `certificates/n8_limit_summary.json`
  (`machine_confirmation.findings.q_upper_end_min_approx`) is a display error in those two files: their own
  per-graph file `notes/verify_n8/v_structure/q_upper_ends_refined.json` has the value above (to be logged in
  ANOMALIES.md; no certificate depends on it; `certificates/n8_summary.json`, `smallest_upper_end_discrepancy_note`).
* **A rational critical value.** `&GMw|Vr]Rb~Nw` and `&GUo~EX~Yxz^w` have `beta_c = 8/3`, `q = 5/3` exactly
  (stored isolating interval [5/3, 5/3], degree 1): their limit polynomial has the factor `4z^2 + 3z + 1`
  with roots of real part exactly -3/8 (`notes/n8_structure.md` 1.8).
* **Lower ends.** 8,162 intervals start at q = 1 (beta = 2) and 971 at an algebraic q > 1; in the limit the
  stable q-intervals of the 9,133 are open intervals in (1, q_8*(limit)) and 8,162 of them start at 1.

**What this does and does not say about the threshold theorem at n = 8.** By Theorem F (reduction,
`notes/lemma_F_reduction.md`) the conjecture at n = 8 fails at a legal parameter iff some non-complete
graph on k <= 8 nodes has a stable positive full support there. Every one of the 9,133 is stable at legal
points with q arbitrarily close to any value in its interval (Lemma I.3/I.4), so **the conjecture fails at
n = 8 for q in a set whose closure contains (1, q_8*(limit))**: for q in (1, rho) by the n <= 7 map (the
failure set is monotone in n, Theorem F), and for q in [rho, q_8*(limit)) by the limit interval of
`&GMx|\Aw^qz]w`, whose lower q end lies in [1, 5/4] (`n8/n8_limit_HITS.jsonl`) with 5/4 < 49/37 < rho
(`notes/referee_finish_overclaims_checks.json`, block D; no file records whether the union of all 9,133
intervals covers (1, q_8*(limit)), which the paper leaves as a question); whether it holds for q >= q_8*(limit) is
**open**, because the finite-eps map is not done: at n = 7 the finite-eps supremum equalled the limit one
(bounded cells gave q < 9/7 < rho), but 40 of the 124 n = 7 regions bulge past beta_c at finite t
(Lemma AS), and the n = 8 bulges have not been computed. The finite-eps `q_8*` satisfies
`q_8* >= q_8*(limit)`.

## 5. Structure: what the 9,133 share (`notes/n8_structure.md`, exact tools in `notes/n8_structure_tools/`)

Facts checked on all 9,133 records; statements marked *proved* have proofs cited in that note.

* **Not near-cliques, and they cannot be.** 11 to 27 missing arcs; the missing-arc digraph D = J - I - A is
  never a matching, a star, a star forest or a forest, and always has a directed cycle (*proved*: a D-source
  with a D-out-arc forces a negative resolvent vector, Lemma D2 consequence). The largest bidirected clique
  has at most 5 = n - 3 vertices (61 attain it); equivalently the missing arcs are never covered by two
  vertices. No graph contains an induced K6 or K7.
* **Universal features (hold on all 9,133):** strongly connected, indeed Hamiltonian (U1); contain an
  induced oriented 3-cycle, the only induced subgraph on 3 to 7 vertices common to all (U2); clique bound
  n - 3 (U3); a node receiving from all others also sends to all others (U4, *proved*); no non-adjacent
  in-twins (U5, *proved*); more bidirected than non-adjacent pairs (U6, *proved*, tight for 2 graphs);
  D has a directed cycle (U7, *proved*). **They are not a signature**: the same tests on the 338,505 strict
  spectral survivors that are not counterexamples hold for all of them except U1, which fails for 18
  (all `never_positivity`). What separates the 9,133 is the joint condition, a common beta at which
  positivity and `Re spec(D|perp) > -1/beta` hold together; no graph-local feature tested captures it
  (`discriminate.json`, `ham_never.json`).
* **Universal statements that fail**, with exact exception counts: has a node sending to all others (5,276
  exceptions); underlying graph complete (3,096); contains an n = 6 or n = 7 counterexample (6,720); stable
  interval starts at beta = 2 (971); ends at the spectral boundary (1,942); trivial automorphism group
  (188); every non-sink node on an induced oriented 3-cycle (2,311); one-way arcs among non-sink nodes
  strongly connected (682); minimum out-degree >= 2 (59) / in-degree >= 2 (1); some vertex deletion a strict
  n = 7 survivor (5,902); `psi_D(0) = 0` (2,951); D a matching / star / star forest / forest (9,133 each).
* **Sink classes (nodes sending to all others; Lemma X, Theorems Q, Q2).** None in 5,276 graphs, 1 in
  3,185, 2 in 562, 3 in 110 (class patterns [] 5,276, [1] 3,185, [1,1] 513, [1,1,1] 71, [2] 49, [2,1] 36,
  [3] 3). For all 3,857 graphs with a sink node, Lemma X2' holds at the certified interior beta: the core
  (graph minus all sink nodes) has a positive fixed point and every sink class has in-set weight above 1
  (3,857/3,857). In 1,885 of them some sink deletion leaves an n = 7 limit counterexample; in the other 1,972 no
  sink deletion leaves an n = 7 limit counterexample (in 2 of these, `&GNbR^x~[zR[?` and `&GNbR^z~[rR[_`, a
  sink deletion leaves the finite-eps-only `&FNFL~[udo?`; 1,970 have no sink deletion that is any n = 7
  counterexample; `notes/referee_finish_overclaims_checks.json`, block A). Cores: 596 distinct 7-node cores (1,356 graphs have
  an n = 7 counterexample as core, 1,829 do not), 40 six-node cores, 5 five-node cores.
* **Instances of the proved constructions.** F_2 of Theorem F (`&GZy~VwuKW~^w`, canonical hit
  `&GNbr^zz\z^[?`) is in the list with its interval ending at the plastic number (stored critical interval
  [86/37, 7/3], minimal polynomial `x^3 - 3x^2 + 2x - 1` in beta), as Theorem F predicts for k = 2. All 626
  n = 8 members of the clique-class catalogue (`families_catalogue_n6/_n7.json`: n = 6 bases with k = 2,
  n = 7 bases with k = 1) are in the list as 593 distinct graphs, 0 missing. The 15 "cube" finite-eps-only n = 8 graphs of the theory notes (E_1 of Theorem E is one of them,
  isomorphic to `&GNQf^zr\IsNw`; alpha = -1/2 exactly; each with an exact legal witness) are **absent**, as
  they must be. Theorem Y's 8-node core `&G^ZZZr|KjqTo` is absent (it needs its sink classes). Geneson's
  double closed-twin blow-ups are absent, as expected (its single blow-ups are limit-never at n = 7). The
  incidental n = 8 candidate `&GLjlNZ}NzH\o` of REPORT_n7 section 10 is present.
* **The 16 graphs with q above 3/2** fall into 10 D-classes; all are dense (37..42 arcs), rigid (trivial
  Aut, no closed twins), 14 with no non-adjacent pair; 10 end at the spectral boundary and 6 at a
  positivity root; only 2 contain an n = 7 counterexample. **The extremal pair** `&GMx|\Aw^qz]w`,
  `&GVrm]A{^pZ]w` (40 arcs, 12 bidirected pairs, no non-adjacent pair, no sink node, D-equivalent, one
  in-twin switch apart) is positive for every beta above the root of `b^4 - 3b^3 + 3b^2 - 3b + 1` in (2, 3)
  (about 2.15372, display) and spectrally limited: `psi_D = z (8z^6 + 16z^5 + 34z^4 + 33z^3 + 30z^2 + 13z + 5)`
  with minimum real part of the roots about -0.36890 (display), so beta_c about 2.71074. Theorems Q and
  Q2 obtain large q from a sink class whose in-set is "all but one node", reaching -0.3516 (Q) and -1/4
  (Q2) only as k -> infinity; both limits lie to the right of -0.3689 (further from -1/2, larger q), so the
  n = 8 record does **not** beat them: at n = 8 the record is attained by a graph with no sink class at all,
  and Theorems Q and Q2 exceed it only in the limit of large class size (Q needs k >= 286 to reach q close
  to 9/5).
  The q-maximisers at n = 7 (the six plastic-number graphs) were positivity-limited; at n = 8 the top 8 of
  the 16 are spectrally limited.
* **Conjecture S** (labelled as such; not proved). Every non-complete limit counterexample on n vertices
  has (S1) no bidirected clique on n - 2 vertices (`omega_bi <= n - 3`); (S2) a Hamiltonian directed cycle,
  in particular is strongly connected; (S3) an induced oriented 3-cycle. Evidence: the complete limit maps
  at n = 6, 7, 8 (S1 sharp at each n: 2, 16, 61 graphs attain n - 3); all 631 proved infinite families for
  every class size (S1 holds, with equality for all k in 38 of them, so n - 3 cannot be lowered to n - 4);
  104 certified instances from n = 7 to n = 73; and an **exhaustive test of S1 at n = 9** on the whole class
  where it could fail, the digraphs with a bidirected K7 (682,176 configurations: 600,440 dominated, 81,735
  fail the strict spectral test, exactly one survives, K9), and likewise at n = 3..8
  (`clique_class_test_n{3..9}.json`; decisions by `graphs.is_domination_free`, three-valued
  `ctln.hurwitz_stable` and `fastregion.region`). So S1 holds for every n <= 9 in the limit sense. The note
  gives a finite reduction that may make S1 provable (with a bidirected K_{n-2}, D has rank at most 4 and
  everything depends on 16 pattern counts and one relation, the kind of finite integer family the
  Hurwitz-minors-in-k proofs of Theorems F, Q, A handle). S2 and S3 are untested at n = 9.

## 6. Verification status

* **Engine, n <= 7 maps, theory:** as in REPORT_n7 section 8 (audited; the 124 n = 7 graphs hand-verified
  2026-10-15, a check the user recorded as made without full confidence (ANOMALIES.md 2026-10-15), and
  confirmed by four independent exact implementations).
* **n = 8, first campaign** (workflow wf_79a996db-7e2, `notes/verify_n8/SYNTHESIS.json`; 2026-10-28;
  ANOMALIES.md 2026-10-28 "ALL 9,133 ... MACHINE-CONFIRMED"): four verifiers written from the project instructions file's
  definitions without the project's code, each with its own digraph6 decoder. `v_fraction` (pure
  Fractions, Faddeev-LeVerrier and Berkowitz characteristic polynomials asserted equal, Hurwitz minors
  three-valued cross-checked by Routh; all four stored witnesses plus an own interior point per graph,
  45,665 points); `v_roots` (exact Routh table and Cauchy/Sturm index on 18,282 points, sympy complex
  root isolation on every 6th graph and the 37 narrow ones, calibrated on 600 random polynomials);
  `v_structure` (own canonical form and `labelg`: pairwise non-isomorphic, none complete, missing arcs
  11..27, all domination-free, witness fields consistent, an own legal stable positive point per graph, q
  upper ends refined); `v_sympy` (sympy exact: first legal witness on all, all witnesses on every 10th,
  endpoint sharpness on every 50th, 183/183 sharp, Sturm count 1 on every stored isolating interval).
  **Result: 9,133/9,133 confirmed by each of the four, 0 refuted, 0 undecided, 0 vanishing Hurwitz
  minors, 0 disagreements between verifiers or with the 36,532 stored witness flags.** Controls: K8
  stable, Geneson's graph confirmed, 20 random never-survivors unstable or non-positive at 280 legal ray
  points. Not certified by that campaign: completeness of the list (the search was not re-run), the whole
  intervals beyond the certified isolating intervals, and Lemma I itself.
* **Human decision** (ANOMALIES.md 2026-10-28): the user declined hand-verification of the n = 8 list
  ("I don't want to hand verify it, just run another verification yourself") and ordered a second
  independent machine campaign in its place. **There is no human hand check at n = 8.**
* **n = 8, second campaign** (workflow wf_6abae796-7df, `notes/verify_n8_2/SYNTHESIS.json`; agent
  directories `completeness_twins`, `interval_all_t`, `never_side`, `dynamics_float`, `chain_referee`,
  `handcheck_kit`), lenses the first campaign did not use:
  * completeness from the proved side: closed-twin blow-ups of the 123 n = 7 counterexamples (858
    distinct: 61 in the list, 797 not limit counterexamples; 4,782 sampled legal points, 0 stable and positive) and the
    Lemma J in-twin-switch closure (29,116 switches of 4,271 records, all in the list with identical
    intervals);
  * stability for **all** large t: for every one of the 9,133, the exact set of t on the limit ray where
    all Hurwitz minors, positivity and legality hold is a single interval (t*, infinity) with finite t*
    (maximum about 8,565.5, display); 36,532 stored witness points agree;
  * the never side: 300 `never_positivity` and 300 `never_full` records at 874,996 legal ray points, 0
    stable and positive; an own all-beta replay reproduces every verdict; the `filter8` boundary rule
    re-derived on 36.9 M regenerated digraphs (S/B/non-printed samples of 300 each, 0 disagreements on S
    versus non-S); the first five `never_full` certificates re-evaluated cell by cell (62/62 agree);
  * a float dynamics lens, explicitly non-decisional: 434 candidates integrate to their exact positive
    fixed point where the decay is measurable (302/303), 200 never graphs never settle at a positive
    full-support point; 14 candidates have basins smaller than a 0.3% perturbation (a remark, not a
    defect: local asymptotic stability is what is claimed and is decided exactly);
  * an end-to-end chain referee (section 2);
  * a hand-check kit with six fully worked candidates (`notes/verify_n8_2/handcheck_kit/HANDCHECK.md`).
  **Synthesis: 0 fatal, 0 major, 1 minor documentation defect** (the B lines described as "alpha = -1/2
  exactly" instead of "alpha >= -1/2", corrected in ANOMALIES.md 2026-10-28), 13 explained non-issues;
  the chain referee's one MAJOR was an incomplete verifier run of its own, covered twice elsewhere.
* **Aggregation** (`n8/aggregate_n8.py`, `certificates/n8_summary.json`): recounts every verdict from the
  results file and asserts the tally, the prefilter totals and the accounting identity; decodes all 9,133
  (non-complete, non-clique, domination-free, pairwise non-isomorphic under `labelg`, arcs = stored
  `edges`); re-verifies that each stored q polynomial is the beta polynomial at 1 + q, that every
  isolating interval holds exactly one root, that all 4,298 + 147 distinct minimal polynomials are
  irreducible; refines every upper end to width < 1e-9 by exact bisection and cross-checks it against
  the first campaign's refinement; recomputes the lifts and the twin counts and matches them to the second
  campaign; certifies the ordering of the upper ends.
* **Known wording and display issues, no certificate affected:** (i) the smallest upper end "about
  1.000276" in `notes/verify_n8/SYNTHESIS.json` and `certificates/n8_limit_summary.json` should read about
  1.000122 (section 4); (ii) the gloss `what_an_endpoint_of_kind_root_certifies` in
  `certificates/n8_summary.json` calls `a_k` a coefficient of the characteristic polynomial, whereas in the
  records, in REPORT_n7 section 1 and in Lemma A.2 the `a_i` are the components of `adj(K) 1`, the
  positivity polynomials, which is how the structure note's exact positivity sets read them (0
  inconsistencies on 9,133); (iii) the stored rational q upper bounds of the records are coarse (two are
  exactly 2 although the algebraic endpoints are about 1.0742) and the refined values are the ones to
  quote. The first campaign also found repeated characteristic-polynomial roots off
  the imaginary axis at 404 of 3,112 sampled points (202 graphs; twin structure), harmless for stability.
* **Nondegeneracy** over all 2^8 - 1 supports at a witness point (REPORT_n7 section 0, item 10, for n <= 7)
  has **not** been run on the n = 8 list.

## 7. Literature position (as of 2026-10-28; `notes/literature_position_2026-10-28.md`)

As of 2026-10-28 we know of no published or posted work that classifies n = 8, treats n = 7, exhibits
finite-eps-only counterexamples, or improves Geneson's threshold (arXiv:2607.21396, v1 only, 23 Jul 2026;
his later 2026 postings are extremal combinatorics). A search of arXiv (topic searches "threshold-linear
networks", "combinatorial threshold-linear", "CTLN", "target-free", and author searches for Geneson, Curto,
Morrison, Shiu; the search index reached entries up to 27 Oct 2026) found no CTLN item dated after 23 Jul
2026 apart from a revision of Londono Alvarez, arXiv:2606.07336 (v2, 28 Sep 2026), which concerns
fixed-point composition under low-rank gluings and not target-free cliques; it was read only at the
abstract level, so no overlap claim is made beyond that it is the nearest relative of the twin and blow-up
lemmas. The unrefereed preprint "Six Is Minimal" (Zenodo 10.5281/zenodo.23089682, v1.1.0-candidate,
2026-10-01) and its repository (HEAD e28a309, re-checked by `git ls-remote`), and the ChickenNungets/CTLN
repository (HEAD 52ad000), were unchanged; the comparison with "Six Is Minimal" in REPORT_n7 section 9
stands unchanged (it covers n <= 5, Geneson's graph and a family on 5r + 6 nodes; we claim no priority
either way). Limits of the check: arXiv's index freshness (27 Oct), HTTP errors on the arXiv API and some
author pages, no search of other preprint servers; repeat before any posting. Nothing was contacted,
submitted or posted.

## 8. What is next

1. **The finite-eps map at n = 8.** The n = 7 pipeline (`n7/REPRODUCE_n7_finite.md`: exact positivity
   prefilter, resolution of undecided lines, e_2 bound, the Lemma R curve test, two-parameter CAD on the
   curve-stable graphs) transplants to n = 8 on the same restricted space, since the e_2 bound holds at every
   legal parameter. Cost, **extrapolated** from the measured n = 7 and n = 8 figures (every input number from
   the file named; the per-line costs at n = 8 will be somewhat higher than at n = 7 because the matrices are
   8 x 8):
   * *Route A, complete map.* (a) A positivity prefilter pass over the 801,558,826,511 restricted digraphs:
     the enumeration itself is the cost, as for the limit prefilter: about 5.5-6 days on this four-core
     machine at the measured aggregate rate of 1.5-1.7 M digraphs/s (ANOMALIES.md 2026-10-16), roughly
     525-595 CPU-hours (131-148 wall-hours times four cores), with the K8 and high-symmetry graphs again through `n8/k8decomp` and
     `n8/orbitenum`. (b) The curve test on the positivity survivors. At n = 7, 8,909,123 of the 117,232,922
     domination-free digraphs survived positivity (`certificates/n7_finite_summary.json`), a ratio of about
     7.6%; **if** the ratio were the same at n = 8 the survivors would number about 8.6 e9 (an assumption,
     not a measurement). The n = 7 curve stage cost 16.4 CPU-hours for the 3,773,348 lines it actually evaluated (REPORT_n7
     section 2), about 15.6 ms per evaluated line (REPORT_n7 quotes the same 16.4 CPU-hours as 6.6 ms per line
     over all 8,909,123 positivity survivors, 5,135,775 of which the edge bound dropped at no cost; the
     evaluated-line figure is the right denominator for this extrapolation); at that rate 8.6 e9 lines are about 37,000 CPU-hours, i.e. about a year on
     four cores or about two to three weeks on a 96-vCPU machine, before the n = 8 slowdown. (c) The
     two-parameter CAD on the curve-stable graphs at about 40 s each (ANOMALIES.md 2026-10-15, "m = 8
     CAD ~40 s each"); their number is unknown (at n = 7 it was 165 of 3,773,348 curve lines).
   * *Route B, the FNFL question only.* Regenerate the prefilter keeping the B lines (cost (a), about 2 GB
     per `n7/REPRODUCE_n8_limit.md`), then run the curve test on the 417,603,634 B lines plus the 7,078
     `never_full` S lines: about 1,800 CPU-hours at 15.6 ms per line (about 19 days on four cores, about
     a day on 96 vCPUs), plus CAD on the curve-stable ones. This lists the finite-eps-only
     counterexamples of the Theorem E type (alpha = -1/2 exactly) among the B and `never_full` lines (15 such
     graphs are already known, E_1 among them; which others exist is open), which
     is where the n = 7 experience (`&FNFL~[udo?`, family E) puts them; it is a complete finite-eps map
     only together with a proof that `alpha > -1/2` strictly excludes stability at every legal point, which
     the project has as data at n = 7 (none of the 124 has alpha > -1/2) and not as a theorem (Lemma
     I.3(c2) gives instability for all large t only). Route B first, then A if the budget allows.
   * Route A also yields the finite-eps `q_8*` (bulges, Lemma AS), and Route B does so only under the same
     proviso on the `alpha > -1/2` lines; either gives the n = 8 analogue of the
     bounded-cell count (40 of 124 at n = 7).
2. **n = 9 feasibility.** The e_2 bound restricts n = 9 to digraphs with at least 37 arcs. By the same
   Burnside count that reproduces A000273(5..8) and the n = 8 restricted count exactly
   (`notes/n9_feasibility_count.py`, `notes/n9_feasibility_count.json`: 9,608; 1,540,944; 882,033,440;
   1,793,359,192,848; 801,558,826,511, all exact), there are **5,903,691,602,470,371** such digraphs up to
   isomorphism at n = 9 (of 13,027,956,824,399,552 in all), 7,365 times the n = 8 restricted space. At the
   measured n = 8 rate of 0.4 M digraphs/s per pipeline (optimistic for n = 9, whose matrices are larger) that
   is at least about 4.1 million CPU-hours (about 43,000 days on four
   cores at 1.6 M/s): **a complete n = 9 limit map by enumeration is out of reach** by this method, even
   rented; it would need a prune of several orders of magnitude (Conjecture S1, if proved, removes only the
   bidirected-K7 class) or a different enumeration (by underlying graph with per-graph Burnside acceptance,
   as at n = 8, plus a much faster filter). What is feasible at n = 9 (`notes/n8_structure.md` section 3):
   (a) all one-node extensions of the 9,133: 9,133 x 2^16 = 598,540,288 labelled graphs, prefiltered by
   domination and the strict spectral test in C and classified by `fastregion.region`; this finds the
   non-primitive n = 9 counterexamples and tests S1-S3 on them, but by the heredity numbers (5,902 of 9,133
   n = 8 graphs are spectrally primitive) we expect, but do not know, that most n = 9 counterexamples are not extensions; (b) the
   cheapest falsification attempt of S3: enumerate the digraphs whose one-way part has no directed
   triangle and prefilter them (at n = 8 no strict survivor at all lacks an induced oriented 3-cycle);
   (c) S1 at n = 9 is already done (the bidirected-K7 class, 682,176 configurations).
3. **Theory.** Prove S1 via the finite reduction (rank <= 4 non-edge digraph, 16 pattern counts); decide
   whether the growth 3, 123, 9,133 and the primitive share 3/3, 102/123, 6,720/9,133 have a structural
   explanation; characterise the limit-stable graphs (Conjecture P was refuted at n = 7, and the n = 8 data
   say the universal features U1-U7 are shared by almost all strict survivors and so are not a signature;
   that U1-U3 follow from the necessary conditions is not proved); the
   growth rate of q_n* (the n = 8 record comes from a sink-free core, the asymptotic families from sink
   classes); nondegeneracy certificates for the 9,133 (`n7/nondegeneracy_check.py`).
4. **Formalisation.** Lemma I (the spectral limit) and the k8 decomposition lemma are the two steps of the
   n = 8 chain that are not in Lean (section 9).
5. **Housekeeping for the orchestrator** (files this report may not modify): log in ANOMALIES.md the
   smallest-upper-end display error and the `a_k` gloss (section 6), the rational critical value 8/3, the
   heredity failure (5,902), the 125 blow-ups of non-counterexamples, and the sink-free record; fix the two
   display strings.

## 9. Machine-checked statements (Lean 4 + Mathlib)

`notes/lean/CtlnCert/` is unchanged by the n = 8 work: 87 theorems, no `sorry`, axiom audit clean
(REPORT_n7 section 11). The n = 8 chain uses two of its lemmas as prunes: the **e_2 edge bound**
(`notes/lean/CtlnCert/CtlnCert/EdgeBound.lean`: Hurwitz stability at a legal point forces more bidirected
than non-adjacent pairs and at least m(m - 1)/2 + 1 arcs; this is what restricts n = 8 to the
801,558,826,511 digraphs with >= 29 arcs) and **Lemma D** (a dominated node forces a zero coordinate at
every fixed point). Also machine-checked and used: Lemma A.2 (positivity is t-free), the closed form of
`I - W`, `q > 1 <-> t(beta - 2) > 1`, Hurwitz stability of `-I + W` iff positive stability of K, Stodola,
Remark J.4'. **Not formalised** and relied on at n = 8: Lemma I (the spectral limit criterion
`alpha(A|perp) < -1/2`, including I.3(c2), which makes the B exclusion sound; refereed in the n = 6 phase
and re-derived by the second campaign), the bidirected-subgraph decomposition lemma of
`notes/k8_decomposition.md` (refereed by four agents, validated by multiset equality at n = 4..7 and the
K8 - e content check), the Lemma J equivalence used for the D-classes, and the converse of Lemma R.
