# Final summary: the target-free clique conjecture, n <= 7 and beyond

Date: 2026-10-15. Start here; details in `notes/REPORT_n6.md`, `notes/REPORT_n7.md`, the lemma notes
and `certificates/INDEX.md`. Every number below is backed by a re-runnable exact certificate.

## 1. The conjecture and what is now known

Curto-Morrison Open Question 10: in a CTLN, the supports of stable fixed points are exactly the
target-free cliques. Geneson (arXiv:2607.21396, July 2026) gave a 6-node counterexample and proved the
conjecture for q >= n - 2 - (n - 3) eps/2, where q = delta(1 - eps)/eps > 1 is the legality parameter.

| n | result | threshold (sup of q over counterexamples) | status |
|---|---|---|---|
| <= 5 | no counterexample at any legal parameter (also reproduced the published n = 5 classification) | none | certified |
| 6 | exactly 3 non-clique stable full supports among 1,540,944 digraphs (Geneson's graph, its bidirected in-twin partner, one new graph) | q_6* = root of 797q^6 - 1686q^5 + 2159q^4 - 3004q^3 + 2499q^2 - 2158q + 1521, about 1.1664 | certified, audited, hand-verified |
| 7 | exactly 124 among 882,033,440 digraphs, at every legal parameter; 123 visible in the eps -> 0 limit, one (`&FNFL~[udo?`) only at finite eps | q_7* = the plastic number, root of q^3 - q - 1, about 1.3247, not attained (two independent exact derivations) | certified, audited, hand-verified, confirmed by four independent implementations |
| all n >= 6 | counterexamples exist for every n (Theorem F); an infinite family invisible in the limit (Theorem E); the critical q is unbounded: q_n* -> infinity (Theorem U) | q_n* >= 1709/1001 (n = 8), 17691/10010 (n = 9), 18418/10001 (n = 10), 25000/10001 (n = 47) | proved and refereed; instances certified (n >= 8 not hand-verified) |

So for n <= 7 the conjecture holds at every legal parameter with q >= 1.3247 and fails for every
q in (1, 1.3247); Geneson's safe region q >= n - 2 is replaced by a constant at these n, while the
families show the true threshold grows with n without any known bound.

## 2. Structural facts

* **Stability is not monotone in the limit direction.** "Stable at finite eps => limit-stable at the
  same beta" is false (exact certificates), so the eps -> 0 classification is not the whole story:
  40 of the 124 seven-node regions bulge past beta_c and are bounded in t, exactly when R* > 0 and the
  fixed point is positive at beta_c (Lemma AS; exact R* for all 124), and one graph exists only at
  finite eps (spectral abscissa exactly -1/2, a moving eigenvalue pair on the boundary line).
* **Ray reduction (Lemma R).** Stability at any legal parameter forces stability on the legality
  curve t(beta - 2) = 1, so "stable somewhere" is an exact one-dimensional test; this made the
  complete finite-eps map at n = 7 a 16-CPU-hour computation instead of thousands.
* **Necessary conditions valid at every parameter** (all machine-checked in Lean): a dominated node
  forbids a positive full-support fixed point (Lemma D); stability forces det K > 0 and makes
  positivity t-free (Lemma A.2); a stable full support has more bidirected than non-adjacent pairs,
  hence more than half of all arcs (e_2 bound); every coefficient of det(mu I + K) is positive
  (Stodola); two non-adjacent in-twins force the eigenvalue delta (never stable).
* **Closed-twin blow-ups** preserve the full-support fixed-point data (Lemmas W, X), which is how
  every infinite family is built and why E_1 at n = 8 has exactly the same region as its 7-node base.
* **Refuted on the way:** Conjecture P (positivity of the resolvent vector forces an eigenvalue at
  Re <= -1/2), and "every limit polynomial has a root on Re = -1/2".

## 3. Verification chain

1. Validation gate: Geneson's counterexample reproduced (permanent test); the published n = 5
   classification reproduced with diff zero; the enumerations reproduce OEIS A000273 through n = 7.
2. Exact engine: Hurwitz determinants on exact characteristic polynomials, three-valued (a vanishing
   minor is undecided, never a verdict); cylindrical decomposition on a proved boundary set; every root
   isolation replayed; pointwise cross-checks at every cell. 0 escalations at n = 5, 6, 7.
3. Independent audits (agents without access to the project's conclusions): engine and n = 6 map;
   the n = 7 prefilter, resolver, curve test, edge bound; every theorem (Theorems F, E, Q, Q2, Y, the
   catalogue, the q* tool) refereed, all holding with high confidence.
4. Human sign-off: n = 6 (2026-10-14) and all 124 n = 7 graphs (2026-10-15), the latter additionally
   confirmed by four separately written exact verifiers at 742 points (0 undecided).
5. Lean 4 + Mathlib: 87 theorems, no `sorry`, standard axioms only (`notes/lean/`).
6. Nondegeneracy certified at a witness point of every hit; literature checked (not scooped; an
   unrefereed preprint independently confirms q_6* and the n <= 5 exclusion).

## 4. Open and pinned

* **Pinned (needs rented compute):** the n = 8 limit map, 500-800 CPU-hours, runbook
  `n7/REPRODUCE_n8_limit.md`, overflow-safe filter validated; one n = 8 limit candidate already found.
* **Open question A** (shared with the preprint): at one FIXED legal (eps, delta), are there stable
  non-clique supports of unbounded size? All known families use size-dependent parameters.
  (prover + referee workflow in progress at the time of writing; result appended below when it lands)
* **Question B, settled (Theorem U, refereed):** q_n* -> infinity. For every N there is n0(N) with
  q_n* >= N - 1 for all n >= n0(N); certified thresholds q_n* >= 2, 3, 4, 9/2, 23/5 for n >= 1812, 219,
  657, 2475, 4909. The growth rate (Geneson's n - 2 versus these super-exponential n0) is open.
* Formalisation gaps: the converse of Lemma R (openness of stability) and the spectral limit
  lemma (Lemma I) are not in Lean.
* Next sizes at finite eps are out of reach by enumeration (n = 8 has 1.79 trillion digraphs); the
  structural route (families, Lemma R, the blow-up calculus) is the way forward.

## 5. Where everything is

`certificates/INDEX.md` lists every certificate and the script that regenerates it;
`ANOMALIES.md` is the append-only log of every surprise and correction (read it before trusting any
number in a note that it supersedes); `notes/REPORT_n6.md` and `notes/REPORT_n7.md` are the full
reports; `n7/REPRODUCE_n7_finite.md` and `n7/REPRODUCE_n8_limit.md` are the runbooks;
`n7/handcheck.py` prints one candidate at one parameter point for checking by hand.
