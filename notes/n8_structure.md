# What the n = 8 limit counterexamples share: structure of the 9,133

Date: 2026-10-28. Conventions are those of the project instructions file, section 2: `A[i][j] = 1` iff `i -> j`. The **non-edge
(missing-arc) digraph** is `D = J - I - A`, so `D[i][j] = 1` iff `i != j` and `i -/-> j`. This is the `D` of
`notes/research_wf3/theory/notes_theory.md` section 1. A node that **sends to all** others is a D-sink: its
D-row is zero, it is a class node with `T = all` in Lemma X, and it is what the theory notes call a "sink class"
node. A node that **receives from all** others is a D-source. All arithmetic in the tools is exact (integers,
`fractions.Fraction`, sympy exact root isolation, nauty `labelg` canonical labels). Decimals appear only for
display and are marked "(display)".

**Scope and status of the input.** The input is the eps -> 0 **limit** map at n = 8: 9,133 non-complete
limit-stable full supports (`n8/n8_limit_HITS.jsonl`). It is VERIFIED BY MACHINE by two independent campaigns
and has **no human hand check**. Two kinds of counterexample are outside it: those visible only at finite eps,
and any that sit among the 417,603,634 unstored boundary lines (`certificates/n8_limit_summary.json`, `scope`).
Every statement below is about this list. Where something is proved, the proof is cited; otherwise it is a
count.

**Tools.** Everything is in `notes/n8_structure_tools/` and runs from that directory with
`nice -n 10 python3 <script>` (at most 2 worker processes). Section 5 gives the order. Each count below names
the JSON file it comes from.

---

## 0. Summary

1. **A large share of the n = 8 counterexamples are new.** 6,720 of the 9,133 (73.6%) contain no n = 6 or
   n = 7 counterexample as an induced subgraph. For 5,902 of them, no single-vertex deletion even passes
   the two necessary conditions at n = 7 (domination-freeness and `alpha(A|perp) < -1/2`). Only 61 are
   closed-twin lifts of n = 7 counterexamples. Only 593 are members of the proved clique-class catalogue
   (`motifs.json`, `heredity.json`, `families.json`).
2. **They are far from cliques and also cannot be close to one.** Between 11 and 27 of the 56 arcs are
   missing. The missing-arc digraph is never a matching, a star, a star forest or a forest. The largest
   bidirected clique has at most **5 = n - 3** vertices, and the bound is attained (61 graphs). Equivalently,
   the missing arcs can never be covered by two vertices (`structure_summary.json`).
3. **Universal features hold on all 9,133.** Every candidate is strongly connected, and in fact Hamiltonian.
   Every candidate contains an induced oriented 3-cycle, which is the only induced subgraph on 3 to 7 vertices
   common to all of them. Any node that receives from all others also sends to all others. There are no
   non-adjacent in-twins. However, these features are equally universal among the 338,505 strict spectral
   survivors that are *not* counterexamples (strong connectivity fails for 18 of them), so they are shared by
   the graphs that pass the necessary conditions (whether they are implied by them is not proved for U1-U3)
   and do not distinguish counterexamples (`discriminate.json`, `ham_never.json`).
4. **Sink classes behave as Lemma X2' predicts.** 3,857 candidates have a node that sends to all others.
   For all 3,857, at a certified interior beta, the core is positive and every sink class is activated
   (`sigma_S > 1`) (`analysis2.json`). In 1,972 of them, deleting any single sink node leaves a graph that is not
   an n = 7 limit counterexample (in 2 of these, `&GNbR^x~[zR[?` and `&GNbR^z~[rR[_`, a sink deletion is the
   finite-eps-only `&FNFL~[udo?`; `notes/referee_finish_overclaims_checks.json`, block A): the sinks create
   the stability (`families.json`).
5. **The largest q at n = 8 comes from a graph with no sink class.** The two extremal graphs `&GMx|\Aw^qz]w` and
   `&GVrm]A{^pZ]w` are one bidirected in-twin switch apart. Neither has a node that sends to all others.
   Their positivity set is the whole half-line beyond a quartic root. Their q is limited only by the spectrum
   of `D|perp`, whose minimum real part is about -0.3689 (display). Theorem Q reaches about -0.3516 (display), which
   lies further from -1/2 (larger q), only in the limit k -> infinity of a sink class; at n = 8 itself the
   record is set without any sink class. Two graphs have the **rational** spectral end beta_c = 8/3
   (q = 5/3), from the factor `4z^2 + 3z + 1` of the limit polynomial (section 1.8).
6. **Conjecture (section 3).** Every non-complete limit counterexample on n vertices has bidirected clique
   number at most n - 3, is strongly connected (Hamiltonian), and contains an induced oriented 3-cycle. Part (i)
   holds at n = 6, 7, 8 and in every proved family for all k. It is now checked **exhaustively at n = 9** on the
   whole class it would have to fail on (all 682,176 configurations with a bidirected K7: only K9 survives). In
   its stronger form it is also checked at n = 3..8.

---

## 1. (a) Facts verified on the data (all 9,133, exact)

### 1.1 Arc and pair counts; the e_2 bound is attained

Source: `structure_summary.json` (`n8`).

* Arcs: 29: 2, 30: 27, 31: 120, 32: 269, 33: 562, 34: 1079, 35: 1485, 36: 1406, 37: 1305, 38: 768, 39: 510,
  40: 476, 41: 581, 42: 353, 43: 153, 44: 34, 45: 3. Missing arcs (|D|) therefore range from 11 to 27.
* Unordered pairs (28 of them), split by type:
  * bidirected pairs `b`: 2..17 (2: 8, 3: 36, 4: 202, 5: 433, 6: 856, 7: 1396, 8: 1606, 9: 1388, 10: 934,
    11: 575, 12: 388, 13: 553, 14: 463, 15: 214, 16: 73, 17: 8);
  * single-arc pairs `s`: 10..26;
  * non-adjacent pairs `z`: 0: 6037, 1: 2895, 2: 198, 3: 3.
  So 6,037 candidates have a complete underlying graph, and no candidate has more than 3 non-adjacent pairs.
* The e_2 bound (`b > z`, Lean-checked) holds with every margin `b - z` from 1 to 17. The margin
  distribution equals the arc distribution shifted, because `arcs = 28 + b - z`. **Margin 1 means arcs = 29:
  the bound is attained.**
* **The two 29-arc candidates** (`special_sets.json`, `arcs29`; positivity sets from `positivity.jsonl`):

  | | `&GM?WUZeEie^O` | `&GMBO]At]Xp]?` |
  |---|---|---|
  | (b, s, z) | (3, 23, 2) | (2, 25, 1) |
  | largest bidirected clique | 2 | 2 |
  | out-degrees | 3,2,3,5,3,4,4,5 | 3,3,4,3,5,4,4,3 |
  | in-degrees | 3,3,5,4,4,4,4,2 | 3,3,3,5,4,4,3,4 |
  | nodes sending to all | none | none |
  | twins of any kind | none | none |
  | Aut | trivial | trivial |
  | contains a smaller counterexample | no | no |
  | stable beta-interval | (2, beta_c), beta_c in [2, 123/61] | (2, beta_c), beta_c in [2, 607/303] |
  | how the interval ends | spectral end (Delta) | spectral end (Delta) |
  | refined q upper end | [64800179/63963136, 16200045/15990784] (about 1.01309, display) | [39760951/39714816, 106029203/105906176] (about 1.00116, display) |
  | positivity set (`positivity.jsonl`) | (2, root of b^5 - 3b^4 + 4b^3 - 5b^2 + 2b - 1), about 2.0633 (display) | two components: (2, root of 2b^3 - 5b^2 + 2b - 1) and (root of 2b^5 - 8b^4 + 11b^3 - 8b^2 + 3b - 1, root of b^3 - 3b^2 + 2b - 1 = 1 + plastic number); the stable interval lies in the first |

  Both are **almost oriented**. They have only 2 or 3 bidirected pairs, and the minimum number allowed with
  their 1 or 2 non-adjacent pairs. Neither is a lift, a catalogue member, or contains a smaller
  counterexample. Both live in a q-window just above 1 (q < 1.014, display). The e_2 extreme is reached by
  sparse, rigid, near-tournament graphs. It is not reached by thinned cliques.

### 1.2 The non-edge (missing-arc) digraph D

Source: `structure_summary.json`.

* **Shape.** D is a matching in 0 of 9,133 candidates, a star in 0, a union of stars in 0, and a forest (as an
  undirected graph) in 0. The maximum degree of D's underlying graph is 4: 2, 5: 652, 6: 3689, 7: 4790. D
  touches 8 vertices in 9,002 candidates, 7 in 129 and 6 in 2. The untouched vertices are exactly the
  bidirected-universal nodes (`uncovered_vertices_are_exactly_bidirected_universal: true`).
* **Not coverable by two vertices.** No two vertices of D's underlying graph cover all its edges in any of the
  9,133. This is the same statement as `omega_bi <= 5` (section 1.3): a vertex cover of size c leaves a
  bidirected clique on 8 - c vertices.
* **Cycles.** D always has a directed cycle (9,133/9,133). This is forced: a D-source with a D-out-arc is
  impossible, because `u = (D^T + rI)^{-1} 1` would be negative (theory notes section 1, consequence of Lemma
  D2). In an acyclic non-empty D such a source always exists. The one-way part of A (equivalently, the reversed
  one-way part of D) always has a directed cycle, and always a directed 3-cycle, i.e. an **induced oriented
  3-cycle** of A (9,133/9,133).
* **Connectivity.** D is strongly connected in exactly the 5,276 candidates with no node sending to all
  (`D_strong_iff_no_sink: true`), since a sink node has D-out-degree 0. D restricted to the non-sink nodes is
  strongly connected in 8,616 and not in 517.

### 1.3 The bidirected part

* Largest bidirected clique `omega_bi`: 2: 2202, 3: 5847, 4: 1023, **5: 61, never 6 or more**
  (`structure_summary.json`). No candidate contains an induced K6 or K7 (`motifs.json`: `contains_K6 = 0`,
  `contains_K7 = 0`).
* The answer to "is every candidate a bidirected clique minus a few arcs?" is **no**. Every candidate lacks at
  least 11 arcs. Bidirected pairs range from 2 to 17 of 28, and the 29-arc candidates have only 2 or 3.

### 1.4 Sources, sinks, sink classes, cores

Sources: `structure_summary.json`, `families.json`, `analysis2.json`.

* No candidate has an A-source (in-degree 0) or an A-sink (out-degree 0). Minimum out-degree is 1: 59,
  2: 2119, 3: 5746, 4: 1208, 5: 1. Minimum in-degree is 1: 1, 2: 114, 3: 3788, 4: 5182, 5: 48.
* **Nodes receiving from all others:** 129 candidates have one such node and 2 have two. Each such node also
  sends to all others, i.e. it is bidirected-universal (133 nodes; 0 exceptions). *Proved* by the D-source
  argument of 1.2. These 131 graphs are joins `H v K1` or `H v K2`.
* **Nodes sending to all others (sink-class nodes):** none in 5,276 candidates, 1 in 3,185, 2 in 562, 3 in 110.
  Grouping sinks with the same in-set from the core into classes gives class-size patterns
  [] 5276, [1] 3185, [1,1] 513, [1,1,1] 71, [2] 49, [2,1] 36, [3] 3.
* **Number of core nodes that do NOT send to a sink**, counted per sink node: 0: 133 (the
  bidirected-universal nodes), **1: 1089** (in-set "all but one node", the Theorem Q2 shape), 2: 2225,
  3: 1039, 4: 153.
* **Core = G minus all sink nodes.** Its size is 8: 5276, 7: 3185, 6: 562, 5: 110. There are 596 distinct
  7-node cores, 40 distinct 6-node cores and 5 distinct 5-node cores. Membership:
  * 7-node cores: 1,356 are n = 7 counterexamples, 1,829 are not;
  * 6-node cores: 62 are n = 6 counterexamples, 500 are not;
  * 5-node cores: none are counterexamples (there are none at n = 5).
* **Lemma X2' checked on the data.** For all 3,857 candidates with sinks, at the certified interior beta of the
  stable cell (`certificates/n8_candidates_legal_witnesses.json`): `w_G > 0` (9,133/9,133), `w_core > 0`
  (3,857/3,857), and `sum_{i in S} w_core,i > 1` for every sink class (3,857/3,857). No exceptions.
* **Deleting one sink node** (`delete_one_sink`): 1,885 candidates have a sink whose deletion leaves an n = 7
  limit counterexample (one of the 123; `families.py` compares with the limit list only). In the other 1,972
  with sinks, no sink deletion does (in 2 of them a sink deletion is the finite-eps-only `&FNFL~[udo?`, so
  1,970 have no sink deletion that is any n = 7 counterexample): 1,829 with one sink, 130 with two,
  13 with three. Iterated peeling of sink nodes ends with a residual core of 8: 5276, 7: 2824, 6: 729, 5: 287,
  4: 17 vertices (`structure_summary.json`, `peeling`).
* **The limit polynomial has a zero root.** Write `psi_D = chi_D - chi_{D+J}` (= 8 chi_{D|perp}). A candidate
  with k sink nodes always has multiplicity of z = 0 at least k - 1, as Lemma X1' predicts with the factor
  `z^{p-1}`. Counts: 2 sinks -> multiplicity 1: 128, 2: 296, 3: 138; 3 sinks -> 2: 60, 3: 50
  (`analysis2.json`). 2,951 candidates have `psi_D(0) != 0`.

### 1.5 Twins, D-classes, symmetry

Sources: `structure_summary.json`, `analysis2.json`, `consistency_extra.json`.

* **Non-adjacent in-twins:** 0 of 9,133 (proved impossible, Lemma J.4').
* **Bidirected in-twin pairs:** none in 4,862 candidates, 1: 3676, 2: 354, 3: 182, 4: 41, 6: 18. The 4,271
  candidates with at least one pair are exactly the 4,271 used in the Lemma J switch closure of
  `notes/verify_n8_2`.
* **Closed twins** (bidirected, same in- and out-neighbourhoods): 186 candidates (181 with one pair, 2 with two
  pairs, 3 with a class of three). Only 61 of these are lifts of n = 7 counterexamples. For the other 125,
  G - u is not an n = 7 counterexample: in 114 it is not even a strict n = 7 survivor, and in 11 it is a
  'never' graph. **So a closed-twin blow-up of a non-counterexample can be a counterexample.**
* **D-equivalence classes** (Lemma J.2 key: characteristic polynomial plus sorted per-vertex in-walk counts):
  6,260 classes. Sizes: 1: 4939, 2: 824, 3: 2, 4: 374, 5: 16, 6: 3, 7: 1, 8: 84, 10: 2, 12: 1, 14: 4, 16: 8,
  19: 1, 32: 1.
* **Automorphism group:** trivial in 8,945, order 2 in 182, 4 in 3, 6 in 3. 186 of the 188 non-rigid
  candidates have closed twins. The other two, `&GNi{NXLMhnNw` and `&GVQ}Pp~Jy}Nw`, have a non-twin
  automorphism.

### 1.6 Instances of the proven constructions; lifts versus new graphs

Sources: `families.json`, `motifs.json`, `heredity.json`, `special_sets.json`.

* **Theorem F.** F_2 (`&GZy~VwuKW~^w`, n = 8) is in the list as `&GNbr^zz\z^[?`. Its stable interval ends at a
  positivity root with q-minimal polynomial `q^3 - q - 1`: the plastic number, as Theorem F (section 3.2)
  predicts for k = 2.
* **Clique-class catalogue** (Lemma X1/X2, `families_catalogue_n6/_n7.json`). All 626 n = 8 members (the
  n = 6 bases with k = 2 and the n = 7 bases with k = 1) are in the list, as 593 distinct graphs. 0 are
  missing. This is an independent recall check of the limit map on the proved families (it cannot detect graphs
  outside them).
* **Lemma W lifts.** I rebuilt all closed-twin blow-ups of the 123 n = 7 counterexamples myself: 858 distinct
  graphs, of which 61 are in the list, from 52 parents. This is identical to the set found by
  `notes/verify_n8_2/completeness_twins`. The 61 are exactly the candidates with a closed twin u such that
  G - u is an n = 7 counterexample.
* **G[S,T;2] over n = 6:** 12 candidates have a closed-twin pair whose deletion leaves an n = 6
  counterexample.
* **Induced smaller counterexamples.**
  * 2,323 contain an n = 7 counterexample, and all 123 n = 7 counterexamples occur;
  * 5 contain `&FNFL~[udo?`;
  * 308 contain an n = 6 counterexample (all three occur);
  * 2,413 contain one or the other;
  * **6,720 contain none.**
  * In the 2,643 embeddings of an n = 7 counterexample, the eighth node sends to all seven in 2,076 and
    misses an out-arc in 567.
* **Heredity fails badly.** For each candidate G, I classified every G - v against the complete n = 7
  limit map (`n7/n7_limit_results.jsonl`):
  * **5,902 have no vertex whose deletion is even a strict n = 7 survivor** (domination-free with
    `alpha < -1/2`);
  * for 797, the best deletion is a 'never_positivity' survivor;
  * for 111, the best deletion is a 'never_full' survivor;
  * for 2,323, the best deletion is a counterexample.
* **Lifts (61) versus the other 9,072** (`special_sets.json`):

  | | 61 lifts | 9,072 others |
  |---|---|---|
  | arcs | 34..44 (40..42: 37 of 61) | 29..45 |
  | sink nodes | 0: 9, 1: 3, 2: 28, 3: 21 (49 have 2 or more) | 0: 5267, 1: 3182, 2: 534, 3: 89 |
  | omega_bi | 3: 11, 4: 36, 5: 14 | 2: 2202, 3: 5836, 4: 987, 5: 47 |
  | closed twins | all (by definition) | 125 |
  | Aut nontrivial | 61 | 127 |
  | contains an n = 6 counterexample | 21 | 287 |
  | catalogue members | 30 | 563 |
  | interval starts at beta = 2 | 56 | 8,106 |
  | largest refined q lower bound | 7402057/5242880 (about 1.4118, display) | 21526147/12582912 (about 1.7107, display) |
  | median refined q lower bound | 84188887/75497472 (about 1.1151, display) | 51201791/47710208 (about 1.0732, display) |

  The lifts are dense, symmetric, sink-heavy and bidirected-heavy. They carry 14 of the 61 graphs with
  `omega_bi = 5`. The new graphs hold both extremes: the sparsest graphs (29 arcs) and the largest q.

### 1.7 The stable interval: what limits it

Sources: own exact positivity sets in `positivity.jsonl` / `positivity.json`; consistency in
`positivity_check.json`. Each set is obtained by integer interpolation of `det M_b` and `adj(M_b) 1`, then exact
root isolation.

* The positivity set P on (2, infinity) is one interval in 8,984 candidates and two in 149. It is unbounded
  above in 2,133.
* **Consistency.** The stored stable interval lies in one component of P in 9,133/9,133 candidates. Its lower
  end is 2 or a root of P's boundary. Its upper end is a boundary root of P exactly when the stored end is a
  positivity root. 0 inconsistencies.
* **Mechanism of each stable interval:**

  | lower end | upper end | count |
  |---|---|---|
  | 2 | spectral (beta_c), positivity continues, P bounded | 4,690 |
  | 2 | spectral, P unbounded | 1,974 |
  | 2 | positivity root | 1,498 |
  | positivity root > 2 | spectral, P bounded | 504 |
  | positivity root > 2 | positivity root | 444 |
  | positivity root > 2 | spectral, P unbounded | 23 |

  The lower end is never spectral: limit stability always holds right down to beta = 2 when positivity does.
  This is consistent with Lemma I.3: the limit-stable set is `(2, beta_c)`.

### 1.8 The 16 candidates with q upper end above 3/2, and the two extremal ones

Source: `special_sets.json` (`q_gt_3_2`, `q_gt_3_2_list`, `q_gt_3_2_D_classes`); refined ends from
`notes/verify_n8/v_structure/q_upper_ends_refined.json`, where no isolating interval contains 3/2. Exact
anatomy from `extremal.py`.

* The 16 fall into 10 D-classes: four singletons and pairs and one class of four:
  * the extremal pair {`&GMx|\Aw^qz]w`, `&GVrm]A{^pZ]w`};
  * {`&GMw|Vr]Rb~Nw`, `&GUo~EX~Yxz^w`} (q end exactly 5/3);
  * {`&GMw|TrmVb~Nw`, `&GZ_|Vo~Yx^^w`};
  * a class of four with the common end in [22932153/14680064, 1638011/1048576] (about 1.5621, display);
  * six singletons, among them `&GZptUbY\B}^w` (= the record graph `&G^ZnRba\ba^w` of theory notes section 5.1,
    refined end about 1.7094, display) and its neighbour `&GZptVbI\B}^w`.
* Common features of the 16:
  * dense: 37..42 arcs, 9..14 bidirected pairs, 14 have no non-adjacent pair;
  * rigid: all have trivial Aut and no closed twins;
  * 12 carry a bidirected in-twin pair;
  * 11 have exactly one sink node and 5 have none;
  * only 2 contain an n = 7 counterexample, and only 2 are catalogue members.
* **How the interval ends.** 10 end at the spectral boundary and 6 end at a positivity root (the class of
  four and the pair at about 1.5621, display). Only the extremal pair has a positivity set that is unbounded
  above (`special_sets.json`, `positive_for_all_large_beta = 2`; the other 8 spectral ones have bounded P with
  the spectral end inside it). The top 8 of the 16 by q are all spectrally limited.
* **The extremal pair.** Both have 40 arcs, 12 bidirected pairs, no non-adjacent pair, `omega_bi = 4`, no node
  sending to all, and trivial Aut. They are D-equivalent (same Lemma J key) and not isomorphic. B is obtained
  from A by the bidirected in-twin switch of twins {4, 7} at node 1 (or 2) (`consistency_extra.json`).
  Exactly:
  * `det M_b = -4b^7 + 11b^6 - 25b^5 + 29b^4 - 30b^3 + 16b^2 - 8b + 1`;
  * `adj(M_b) 1` has components `-(b^2 - b + 1)^3`, `-(b^2+1)(b^2-b+1)`, ...;
  * the only boundary point of P beyond 2 is the root of `b^4 - 3b^3 + 3b^2 - 3b + 1` in (2, 3), about
    2.15372 (display), so P = (2.15372..., infinity);
  * `psi_D = z (8z^6 + 16z^5 + 34z^4 + 33z^3 + 30z^2 + 13z + 5)`, with minimum real part of the roots about
    -0.36890 (display), so beta_c about 2.71074 (display).

  The stable interval is therefore (2.15372, 2.71074) and q ranges over (about 1.1537, about 1.7107), both
  display. The refined upper end is `[21526147/12582912, 5381537/3145728]`, degree 15.
* **Why q is large, compared with Theorem Q2.** Theorems Q and Q2 split the problem in two:
  1. a core that is positive at every beta > 2, so the ceiling is not set by positivity;
  2. sink classes with in-sets "all but one node", which move the moving spectrum toward
     `spec(D - bbar 1^T)`, whose minimum real part is -0.3516 for Q and -1/4 for Q2 (display), as k -> infinity.

  The extremal n = 8 pair has property 1 in the form "positive for all beta beyond 2.1537". It has property 2
  **with no sink class at all**: the 8-node core's own `D|perp` already has minimum real part -0.3689
  (display). That is further LEFT (closer to -1/2) than the k -> infinity limit of Theorem Q's 7-node core
  (-0.3516); Q overtakes the n = 8 record only for large k (q close to 9/5 needs k >= 286). The large-q
  graphs at n = 8 therefore come from cores alone, not from sink classes, although the cores themselves do
  not beat the asymptotic sink-class limits. This is the "growth comes from new
  non-twin structure" alternative of theory notes 5.3. The one-sink members of the 16 (11 of them) are mixed
  cases. For example, in `&GMw|Vr]Rb~Nw` the sink has in-set "all but two core nodes" (in-set
  {0,1,2,3,5} of core {0..6}).
* **The rational critical value.** In `&GMw|Vr]Rb~Nw` and `&GUo~EX~Yxz^w`,
  `psi_D = z (4z^2 + 3z + 1)(2z^4 + 2z^3 + 4z^2 + 2z + 1)`. The quadratic factor has roots
  `(-3 +- i sqrt 7)/8`, real part exactly -3/8, so **beta_c = 8/3 and q = 5/3 exactly** (stored isolating
  interval [5/3, 5/3], degree 1). Their positivity set is (root of `b^3 - 2b^2 - 1`, about 2.2056 (display);
  root of `b^3 - 3b^2 + b - 1`, about 2.7693 (display)). The upper root 2.7693 is the same positivity ceiling
  as the `&F^V\K[epa?` class family of theory notes 5.1. The stable interval (2.2056..., 8/3) is open, so
  q = 5/3 is not attained.

### 1.9 The 37 narrow-interval graphs

Source: `special_sets.json` (`narrow37`), flag `in_37_no_legal_stored` of `q_upper_ends_refined.json`.

* Every one starts at beta = 2 and ends at the spectral boundary. The largest refined q upper end among them is
  8790397/8781824 (about 1.00098, display): `alpha(A|perp)` is barely below -1/2.
* None is a lift, a catalogue member, or has closed twins, and all are rigid.
* 22 have no sink node and 15 have one. 15 carry a bidirected in-twin pair. 5 contain an n = 7 counterexample
  and none contains an n = 6 one. Arcs range from 31 to 41 (35: 13). Non-adjacent pairs: 0: 27, 1: 10.
* The smallest q upper end of the whole list belongs to this set: `&GMBQMHtMH~]W`, with refined end
  [77079731/77070336, 1079116235/1078984704] (about 1.000122, display). I confirmed it exactly in
  `check_min_q.json`: the graph is limit-stable at beta = 20001/10000 (q = 1.0001) and not at
  beta = 10001/5000 (q = 1.0002). The positivity set reaches 2.3247 (display), so the end is spectral. **This
  contradicts the value "smallest about 1.000276" recorded in `certificates/n8_limit_summary.json`
  (`q_upper_end_min_approx`) and in `notes/verify_n8/SYNTHESIS.json`** (section 4).

### 1.10 Universal statements: which hold and which fail

"Holds" means for all 9,133 candidates. The last column is the same test on the 338,505 non-complete strict
spectral survivors that are NOT counterexamples (`discriminate.json`, `ham_never.json`).

| # | statement P ("every n = 8 counterexample has P") | n = 8 result | n = 7 (123 + FNFL) | n = 6 (3) | on the 338,505 'never' survivors |
|---|---|---|---|---|---|
| U1 | strongly connected; indeed Hamiltonian | **holds** | holds | holds | fails for 18 (all 'never_positivity'); Hamiltonian iff strongly connected there |
| U2 | contains an induced oriented 3-cycle (the only universal induced 3..7-vertex subgraph, `motifs.json`) | **holds** | holds | holds | holds (all) |
| U3 | largest bidirected clique <= 5 (= n - 3) | **holds** (61 attain 5) | holds, <= 4 (16 attain) | holds, <= 3 (2 attain) | holds (all) |
| U4 | a node receiving from all others also sends to all others | **holds** (*proved*) | holds | holds | holds |
| U5 | no non-adjacent in-twins | **holds** (*proved*, J.4') | holds | holds | holds |
| U6 | more bidirected than non-adjacent pairs | **holds** (*proved*, e_2; tight: 2 graphs) | holds | holds | not tested |
| U7 | D has a directed cycle | **holds** (*proved*, 1.2) | holds | holds | holds |
| F1 | has a node sending to all others | fails: **5,276** exceptions | fails: 62 | fails: 1 | |
| F2 | underlying graph complete (no non-adjacent pair) | fails: **3,096** | fails: 57 (+FNFL) | fails: 3 | |
| F3 | contains an n = 6 or n = 7 counterexample | fails: **6,720** | fails: 102 (contains n = 6) | | |
| F4 | stable interval starts at beta = 2 | fails: **971** | | | |
| F5 | stable interval ends at the spectral boundary | fails: **1,942** | | | |
| F6 | trivial automorphism group | fails: **188** | fails: 3 (+FNFL) | holds | |
| F7 | every non-sink node lies on an induced oriented 3-cycle | fails: **2,311** | fails: 33 | holds | |
| F8 | one-way arcs among non-sink nodes form a strongly connected digraph | fails: **682** | fails: 11 | holds | |
| F9 | minimum out-degree >= 2 / minimum in-degree >= 2 | fails: **59** / **1** | fails: 1 / holds | holds | |
| F10 | some vertex deletion is a strict n = 7 survivor | fails: **5,902** | | | |
| F11 | `psi_D(0) = 0` (D|perp singular) | fails: **2,951** | | | |
| F12 | missing-arc digraph is a matching / star / star forest / forest | fails: **9,133** each | fails: all | fails: all | |

Sources: `structure_summary.json` (including `more_universals`), `discriminate.json`, `ham_never.json`,
`n7_strict_props.json`, `motifs.json`, `heredity.json`, `analysis2.json`. For n = 7, all 1,765 'never' strict
survivors also satisfy U1, U2, U3, U4 and Hamiltonicity (`n7_strict_props.json`).

**Reading.** U1–U7 are true but are not the signature of a counterexample. With the single exception of U1
(18 graphs), they already hold for every graph that passes domination-freeness and the strict spectral
condition. What separates the 9,133 from the 338,505 is the **joint** condition: a common beta at which
positivity and `Re spec(D|perp) > -1/beta` hold together. No graph-local feature I tested captures it.

---

## 2. (b) Comparison with n = 7 and with the theory

* **Non-edge picture** (`notes/research_wf3/theory/notes_theory.md` section 1). Every consequence of Lemmas D1
  and D2 that I tested holds on all 9,133:
  * D-sources are D-isolated (U4);
  * D is cyclic (U7);
  * the positivity sets match the stored stable intervals (1.7).

  The non-adjacent pairs grow from at most 1 (n = 6, 7) to at most 3. At n = 7, 66 of 123 had no
  non-adjacent pair; at n = 8, 6,037 of 9,133. The arc range moved from 23..31 of 42 (minimum 22 allowed, not
  attained) to 29..45 of 56 (minimum 29 allowed, **attained**).
* **Bidirected clique bound.** It is n - 3 and attained at every n: max 3 at n = 6 (2 graphs), 4 at n = 7
  (16), 5 at n = 8 (61). REPORT_n7 section 3 noted "none contains K6". The n = 8 analogue is "none contains
  K6", and K6 is now the critical size.
* **Sink classes** (Theorems Q, Q2, Lemma X1'/X2').
  * Share of candidates with at least one node sending to all others: 2 of 3 at n = 6, 61 of 123 at n = 7,
    3,857 of 9,133 at n = 8.
  * Classes of size two or more: 2 at n = 7 (plus FNFL), 88 at n = 8 (49 + 36 + 3).
  * X2' holds on every sink candidate (1.4).
  * In 2,076 of 2,643 embeddings of an n = 7 counterexample the extra node sends to all. At n = 7 it was 21 of
    23 embeddings of an n = 6 counterexample. The "T = all" preference persists but is weaker.
* **Lemma W (closed-twin blow-ups).** 61 of the 858 blow-ups of n = 7 counterexamples survive. Lemma W
  preserves the fixed-point data of the weighted quotient, not the region, and most blow-ups are not limit
  counterexamples. Conversely,
  125 counterexamples are blow-ups of non-counterexamples (1.5). The blow-up calculus is not monotone in
  either direction.
* **Lemma X1/X2 catalogue.** 626 of 626 n = 8 members are present (593 distinct graphs). Theorem F: F_2 is
  present with the plastic-number end, as proved.
* **Theorem E and the finite-eps-only candidates.** The 15 'cube' n = 8 finite-eps-only graphs of theory
  notes 4.4 (E_1 of Theorem E is one of them; exact legal witnesses) have alpha = -1/2 exactly, and **none of
  them is in the limit list**, as it must be (`consistency_extra.json`). They are n = 8 analogues of
  `&FNFL~[udo?` and live among the unstored B lines, as would any others. This analysis does not cover them.
* **Theorem Y.** Its 8-node core `&G^ZZZr|KjqTo` is not itself in the list. It is a core, positive at every
  beta > 2, that needs its sink classes.
* **Theorems Q / Q2 versus the n = 8 record.** See 1.8. At n = 8 the record q comes from a sink-free core
  whose own moving spectrum (minimum real part -0.3689) is closer to -1/2 than Theorem Q's k -> infinity
  limit (-0.3516). At n = 7 the
  q-maximisers (the six plastic-number graphs) were positivity-limited: 3 with a sink node, 3 without
  (computed from `certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl`, `q_minpoly = q^3 - q - 1`). At n = 8
  the top 8 of the 16 are spectrally limited.
* **Scale.** There are 3, 123 and 9,133 limit counterexamples at n = 6, 7, 8. Primitive ones (no smaller
  counterexample induced) number 3, 102 and 6,720. D-classes number 2, 91 (with FNFL) and 6,260. The share of
  candidates with a bidirected in-twin pair is 64/124 at n = 7 and 4,271/9,133 at n = 8.

---

## 3. (c) Conjectures suggested by the data (labelled as such; none is proved)

**Structural Conjecture S.** Let n >= 3 and let G be a non-complete digraph on n vertices whose full support is
limit-stable at some beta > 2. Then:

* **(S1) clique bound:** G has no bidirected clique on n - 2 vertices (`omega_bi <= n - 3`). Equivalently, the
  missing arcs cannot be covered by two vertices.
  **(S1+), stronger:** every non-complete digraph with a bidirected clique on n - 2 vertices is dominated or has
  `alpha(A|perp) >= -1/2`.
* **(S2) connectivity:** G is strongly connected. Even stronger: G has a Hamiltonian directed cycle.
* **(S3) oriented triangle:** G contains an induced oriented 3-cycle.

**Evidence.**

1. Complete limit maps.
   * n = 6, 7, 8: S1, S2 (Hamiltonian form) and S3 hold for all 3 + 123 + 9,133 counterexamples, and for
     `&FNFL~[udo?`.
   * S1 is **sharp** at each n: n - 3 is attained by 2, 16 and 61 graphs.
   * n <= 5: there are no counterexamples.
2. **Every proved infinite family, for every k** (`family_test.json`, symbolic part). For each of the 631
   families (the 626-entry catalogue, Theorems F, Q, A, Q2, Y), `omega_bi` of the member with class sizes k_j is
   `max_J (sum_{j in J} k_j + omega(G[cap_J S_j]))`. The bound `omega_bi <= n - 3` holds for all k with
   slack >= 0. It holds **with equality for all k** in 38 families, so S1 cannot be improved to n - 4.
3. **Certified instances up to n = 73** (`family_test.json`, instances). 104 graphs pass S1, S2, S3 and U4:
   * Theorem F witnesses;
   * the high-q graphs at n = 8..11;
   * the Q, Q2 and Y instances at n = 19, 32, 47, 73;
   * the Theorem U explicit graphs at n = 28, 41;
   * Theorem A's R_k for k = 0..40;
   * the self-joins of R for k = 1..4.
4. **Exhaustive test of S1+ for n = 3..9** (`clique_class_test_n{3..9}.json`). Every digraph with a bidirected
   K_{n-2} is a multiset of n - 2 vertex patterns (16 types relative to the two outside vertices) plus one of 4
   relations between those two vertices. So the enumeration is complete; there are 682,176 configurations at
   n = 9. Results:
   * n = 9: 600,440 are dominated, 81,735 fail the strict spectral test, and **exactly one survives: K9**
     (limit verdict stable, complete);
   * n = 3..8: the same, with only K_n surviving.
   * The pipeline matches the complete maps at n = 6, 7, 8 (no non-complete survivor, as the maps require).
   * The strict spectral test passed all 1,889 n = 7 and 22,719 sampled n = 8 strict survivors
     (`clique_class_control.json`).
   * Decisions use `graphs.is_domination_free`, `ctln.hurwitz_stable` (three-valued; a vanishing minor means
     `alpha >= -1/2`, which is never limit-stable by Lemma I.3(c2)) and `fastregion.region`.
   * **S1 is therefore verified for n <= 9 (limit sense), unconditionally on any n = 9 enumeration.**
5. Against S2, S3: nothing found. S2 at the positivity level: the 18 non-strongly-connected strict survivors
   at n = 8 all fail positivity.

**Why S1 may be provable.** With a bidirected K_{n-2}, D has rank at most 4. Its moving limit spectrum and the
positivity vector depend only on the 16 pattern counts and the x–y relation, through
`Q P = (counts)` (Schur complement over the two outside vertices). Domination-freeness already forbids the
pattern (i <- x one-way, i <- y one-way) and forces strong restrictions on the rest (600,440 of 682,176
configurations die at n = 9). So S1+ is a statement about a finite-type integer family. It is of the kind the
"Hurwitz minors as polynomials in k" proofs of Theorems F, Q and A handle. S2 is related to Curto–Morrison
Q13 ("proper sources die") but is not implied by it: in a CTLN, absent arcs are strong inhibition, so an
upstream set is not dynamically decoupled.

**Weaker observations that the data suggest but I do not elevate to conjecture.**

* The share of primitive counterexamples grows: 3/3, 102/123, 6,720/9,133. Most counterexamples at size n+1
  will not be one-node extensions of those at size n.
* The q record at fixed n is set by sink-free or one-sink cores whose own `D|perp` spectrum is far from -1/2.
  Sink classes are not needed at fixed small n; they matter asymptotically (Theorems Q, Q2, Y, U).

**The obvious tests at n = 9.**

* (a) **Done here:** S1+ on the whole bidirected-K7 class (above).
* (b) **One-node extensions.** Enumerate all one-node extensions of the 9,133: 9,133 x 2^16 = 598,540,288
  labelled graphs. Prefilter with domination and strict spectral tests in C (the n = 8 filter8 code does
  this per graph). Classify the survivors with `fastregion.region`. Check S1, S2 and S3 on every n = 9 limit
  counterexample found.
  * This does not cover primitive n = 9 counterexamples. By the heredity numbers (5,902 of 9,133
    'spectrally primitive' at n = 8) we expect them to be the majority; this is not known.
  * A complete n = 9 test needs the n = 9 limit map, restricted by the e_2 bound to digraphs with at least
    37 arcs.
* (c) **The cheapest full falsification attempt of S3 at n = 9.** Enumerate digraphs without an induced
  oriented 3-cycle, i.e. those whose one-way part has no directed triangle, and run the prefilter. At n = 8,
  no strict survivor at all lacks one, among 347,638.

---

## 4. Discrepancies and items for the anomaly log

These are not appended to `ANOMALIES.md`, which I was told not to modify.

1. **Wrong minimum q upper end in two certificate files.** `certificates/n8_limit_summary.json`
   (`machine_confirmation.findings.q_upper_end_min_approx = "1.000276"`) and `notes/verify_n8/SYNTHESIS.json`
   say "smallest ~1.000276". Their own refined data
   (`notes/verify_n8/v_structure/q_upper_ends_refined.json`) and my exact check (`check_min_q.json`) give
   `&GMBQMHtMH~]W` with end in (10001/10000, 5001/5000): it is stable at q = 1.0001 and unstable at
   q = 1.0002. The candidate list is unaffected. The quoted number is wrong.
2. **A rational critical value:** beta_c = 8/3, q = 5/3 for `&GMw|Vr]Rb~Nw` and `&GUo~EX~Yxz^w`, from the
   factor `4z^2 + 3z + 1` of psi_D (real part -3/8).
3. **Heredity fails.** 5,902 of 9,133 counterexamples have no vertex deletion that passes even the strict
   n = 7 prefilter.
4. **Blow-ups of non-counterexamples are counterexamples:** 125 closed-twin candidates, 114 over bases that
   are not even strict n = 7 survivors.
5. **The q record at n = 8 has no sink class.** Its moving spectrum, minimum real part -0.3689 (display),
   is attained without any sink class, but does NOT beat Theorem Q's k -> infinity limit of -0.3516 (display):
   -0.3689 lies to its left.
6. **Format note.** The n = 73 Theorem Y witness in `notes/research_wf3/theory/family_Y_certificate.json`
   uses a non-standard single-byte digraph6 header (chr(136)). nauty's `N(n)` for n >= 63 is
   `126` plus three bytes. My decoder accepts both and length-checks.

---

## 5. Files and re-running

All in `notes/n8_structure_tools/`. Run in this order (2 workers, `nice -n 10`); total about 25 CPU-minutes
plus about 12 minutes for n = 9.

| script | output | content |
|---|---|---|
| `lib.py` | | own digraph6 codec (standard + long header), nauty `labelg` wrapper, pair types, D, twins, exact charpoly, D-class key, exact limit fixed point, automorphism count |
| `small_lists.py` | | the n = 6 (3) and n = 7 (123 + FNFL) lists from the certificates |
| `features.py` | `features.jsonl` | per-candidate exact features (with the stored interval data) |
| `small_features.py` | `small_features.json` | the same features for n = 6, 7 |
| `summary.py` | `structure_summary.json` | all distributions of sections 1.1–1.5 and 1.10 |
| `discriminate.py`, `ham_never.py` | `discriminate.json`, `ham_never.json` | the universal statements on all 347,638 non-complete strict survivors |
| `motifs.py` | `motifs.json`, `motifs_per_hit.jsonl` | induced-subgraph census k = 3..7, containment of n = 6 / 7 counterexamples, K6, K7 |
| `families.py` | `families.json` | lifts, F_2, catalogue members, G[S,T;2], cores, sink deletions, sink in-set deficiencies |
| `analysis2.py` | `analysis2.json` | D-classes, zero root of psi_D, closed-twin non-lifts, Aut, Lemma X2' check |
| `heredity.py` | `heredity.json` | vertex deletions against the complete n = 7 limit map |
| `more_universals.py` | `more_universals.json` | further universal statements at n = 6, 7, 8 |
| `n7_strict_props.py` | `n7_strict_props.json` | the same on all 1,889 strict n = 7 survivors |
| `positivity.py`, `positivity_check.py` | `positivity.jsonl/.json`, `positivity_check.json` | exact positivity sets and interval mechanisms |
| `special_sets.py` | `special_sets.json` | 29-arc pair, the 16 with q > 3/2, the 37 narrow graphs, lifts vs new |
| `extremal.py d6 ...` | stdout (JSON) | exact anatomy (det, adj, P, psi_D) of named graphs |
| `consistency_extra.py` | `consistency_extra.json` | E_1 / cube candidates absent; the extremal in-twin switch |
| `check_min_q.py` | `check_min_q.json` | exact check of the smallest q upper end |
| `family_test.py` | `family_test.json` | S1–S3 on 104 certified instances (n = 7..73) and on all 631 proved families for all k |
| `clique_class_test.py n` | `clique_class_test_n{n}.json` | exhaustive bidirected-K_{n-2} class test (n = 3..9) |
| `clique_class_control.py` | `clique_class_control.json` | positive control of the strict spectral test |
