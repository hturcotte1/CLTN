# interval_all_t -- exact all-large-t certification of the 9,133 n=8 limit candidates

Agent: interval_all_t (second verification campaign). Sample rule: ALL 9,133 candidates
(i % 1 == 0), i.e. the mandated i%4==0 set + 37 null-witness graphs, extended to every graph
because the measured cost was 0.45 s/graph (2 workers, nice 10, ~34 min wall).

Method (engine.py, own code; run.py driver; tally.py; controls.py):
* digraph6 decoded by own parser; edge count checked against the record.
* beta = exact rational of certificates/n8_candidates_legal_witnesses.json
  'witness_in_first_stable_cell'; certified strictly inside the record's first stable interval by
  exact bisection of the stored minimal polynomial on the isolating interval (beta_in_interval).
* eps = 1/(t+1), delta = (beta-1) eps; (t+1) d (I-W) == K'(t) = d I + b N + d t J, N = J-I-A^T,
  checked SYMBOLICALLY against the W definition for every graph (W_definition_check).
* M = -I+W = -K'/(d(t+1)); char poly of -K' by Faddeev-LeVerrier over Q[t]; Hurwitz minors H_k by
  fraction-free Bareiss over Q[t]; H_k(M) = H_k(-K')/(d(t+1))^{k(k+1)/2}, same sign for t>0.
  Identities a_8 = det K' and H_8 = a_8 H_7 asserted for every graph.
* x(t) = (t+1) d adj(K') 1 / det K' (Cramer numerators N_i, D = det K').
* Conditions: H_1..H_7 > 0, D > 0, N_i > 0 (i=0..7), L = (b-2d)t - d > 0 (q > 1).
  Every condition polynomial factored over Q (sympy), positive real roots of every distinct
  irreducible factor isolated exactly (sympy intervals, refined until factors' intervals are
  disjoint), sign of every condition evaluated exactly (Fraction) in every gap. S = union of
  gaps with all conditions > 0. t* = largest positive root of any condition polynomial.
  A root is never counted as stable (open set; vanishing minor = undecided).
* Cross-checks: cell witness t in S and re-decided pointwise by an independent Fraction engine
  (W built from the definition); Hurwitz-minor signs and x signs of the symbolic polynomials
  vs the pointwise engine at t = 100, 1000, 10000, 100000 at beta_cell; the 4 stored witnesses
  (their own beta) re-decided pointwise and compared flag by flag; the same symbolic t-analysis
  at the stored beta, with the stored stable&positive flags compared with membership of t in the
  stability+positivity set; stored beta inside the limit interval <=> that set contains all large t.

Results (tally.json, table.tsv, results.jsonl):
* 9,133/9,133: beta_cell strictly inside the interval; S contains (t*, infinity); t* finite;
  S is a single interval (t*, inf) for every graph; cell witness t in S; 0 problems, 0 errors,
  0 undecided points. Binding condition at t*: legality L for 8,388 graphs, H7 for 745.
  t* <= 100 for 8,139 graphs, <= 1,000 for 9,054, max t* ~ 8565.50 (idx 5189, degree-6 factor).
* 36,532/36,532 stored witness points agree flag-for-flag with the independent exact re-decision;
  36,532/36,532 consistent with the exact t-set at the stored beta; 9,096 stored betas inside the
  limit interval all give unbounded sets; the 37 null-witness betas are outside and give bounded sets.
* Controls (controls.json): K8 stable for all t>0 at 3 betas; empty graph never; 80 'never' graphs
  from n8_limit_results.jsonl never limit-stable at 4 betas; 12 n=7 counterexamples unbounded at a
  stable-cell beta, bounded at an unstable-cell beta. NOTE: control graph &GMA[KPT^B|Ig
  (verdict never) is stable+positive for t in ~(0.155,0.695) at beta=21/10 -- illegal (q<1),
  consistent with 'never' (limit) verdict; see control_finite_t_note.json.
