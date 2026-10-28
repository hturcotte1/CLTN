# Geneson, arXiv:2607.21396 -- retrieval and verification

Retrieved 2026-10-12/13 by the literature agent. Nothing below is reconstructed from memory or citations.

## 0. Provenance (how the text was obtained)

- https://arxiv.org/abs/2607.21396 : resolved. Title "A six-neuron counterexample to the target-free clique conjecture",
  author Jesse Geneson only, "Submitted on 23 Jul 2026", 23 Jul 2026 14:58:30 UTC, 12 KB, subjects math.CO; cs.DM.
  ONLY v1 is listed (no v2 on the abstract page at fetch time). 12 pages.
- WebFetch of the PDF returned the raw PDF bytes (WebFetch could not read them, but saved the file). I extracted text locally with
  `pdftotext -layout` from that saved file. File: .../scratchpad/pdfwork/geneson.pdf ; text: .../scratchpad/pdfwork/geneson.txt .
  The PDF page-1 margin stamp reads "arXiv:2607.21396v1 [math.CO] 23 Jul 2026". ALL VERBATIM QUOTES BELOW COME FROM THIS pdftotext OUTPUT.
  Caveat: pdftotext flattens displayed math (fractions, matrix brackets, subscripts); I re-flowed those displays by hand and say so where it matters.
  I did not fetch the TeX source (arxiv.org/src/...), so I have not seen the LaTeX.
- WebFetch of https://arxiv.org/html/2607.21396v1 also returned content, but it passes through a summarizing model, so I treated it as a cross-check only.
  Its key facts (matrix N, delta = 29 eps/25, beta = 54/25, Thm 4.1 inequality, Cor 4.2) agree with the PDF text.
- Another WebFetch/WebSearch layer produced a paraphrase of the abstract; I discarded it in favour of the PDF text.

## 1. The explicit 6-node counterexample (Theorem 3.2), verbatim

"Theorem 3.2. Let
  N = [ 0 0 0 0 1 1
        0 0 0 0 0 1
        0 0 0 0 0 1
        1 1 0 0 1 0
        0 0 1 1 0 0
        0 0 0 1 1 0 ]                                   (19)
There is an eps0 > 0 such that, for every 0 < eps < eps0 and every theta > 0, the CTLN with full nonedge matrix N and
  delta = (29/25) eps
is legal and nondegenerate, its full support is permitted, and the corresponding fixed point is stable.
The full support is not a clique, and
  q = (29/25)(1 - eps) -> 29/25   (eps down to 0)."

(The matrix rows above are the pdftotext rows, re-flowed. HTML fetch printed the same six rows.)

The paper gives the graph ONLY as this nonedge matrix N. There is no figure, edge list, or adjacency matrix in the PDF text.

## 2. Convention used by the paper, and conversion to ours

### Verbatim definitions (Section 2)

Dynamics and weights:
"The CTLN associated with G has state x(t) in R^n_{>=0}, weight matrix W in R^{n x n}, and dynamics
  x_i-dot = -x_i + [ sum_{j=1}^n W_ij x_j + theta ]_+ ,  i in [n],  theta > 0,   (2)
where [u]_+ = max{u, 0} and
  W_ij = 0,        i = j,
         -1 + eps, j -> i,
         -1 - delta, j -/-> i,
at parameters  delta > 0,  0 < eps < delta/(1+delta).   (3)
We call (3) the legal range."
Context sentence before that: "Fix a positive integer n, and let G be a directed graph on [n], with no self-loops and at most one edge in each ordered direction between distinct vertices. Edges in both directions are allowed."

Nonedge matrix:
"Fix a nonempty candidate support sigma = {v_1, ..., v_m} subset of [n] and fix the displayed ordering. Its nonedge matrix is the m x m zero-one matrix N defined, for a, b in [m], by
  N_aa = 0,   N_ab = 1  <=>  v_b -/-> v_a   (a != b)."
Then:
"Set  t = (1 - eps)/eps,  beta = 1 + delta/eps,  q = delta(1 - eps)/eps.   (5)"
"Write (I_n - W)_sigma for the principal submatrix indexed by the displayed ordering of sigma, and set
  A_sigma := (I_n - W)_sigma,   K := eps^{-1} A_sigma = I_m + beta N + t J_m.   (6)"
Clique / target:
"A clique is a vertex set joined in both directions. A vertex outside a clique is a target if it receives an edge from every vertex in the clique. A clique is target-free if it has no target."
Legality:
"The legal condition is equivalent to q > 1, and also to t(beta - 2) > 1."
Intro: "the legal CTLN parameter range is q > 1".

### Convention verdict

- Paper's W_ij = -1 + eps iff j -> i. THIS IS IDENTICAL to our CLAUDE.md convention (W[i][j] = -1+eps iff A[j][i]==1, A[i][j]==1 iff i->j). No sign or transpose difference in W.
- Paper's q, t, beta, legal range: identical to ours.
- CAUTION, NAME COLLISION: the paper's symbol A_sigma (eq. 6) is the principal submatrix of I - W, NOT an adjacency matrix. The paper never uses an adjacency matrix A.
- CAUTION: the paper's N is a NONEDGE matrix indexed receiver-first: N[a][b] = 1 iff there is NO edge b -> a. Reading N as our adjacency matrix is wrong twice over (complement AND transpose).
- Conversion (I TRANSPOSED AND COMPLEMENTED): for a != b, A[b][a] = 1 - N[a][b], i.e. A = J - I - N^T. Paper vertices v_1..v_6 = our nodes 1..6 (1-indexed); our code is 0-indexed, subtract 1.
- I verified in exact Fractions, using the repo's own ctln.weight_matrix, that I - W == eps * (I + beta N + t J) entrywise with this A, at eps = 1/50, 1/1000 (delta = 29 eps/25). So the conversion is consistent with the repo's W convention.

### The graph in OUR convention (A[i][j] == 1 iff i -> j), 0-indexed

    A = [[0,1,1,0,1,1],
         [1,0,1,0,1,1],
         [1,1,0,1,0,1],
         [1,1,1,0,0,0],
         [0,1,1,0,0,0],
         [0,0,0,1,1,0]]

1-indexed (paper labels) directed edges i->j (19 of 30 ordered pairs):
 1->2 1->3 1->5 1->6 | 2->1 2->3 2->5 2->6 | 3->1 3->2 3->4 3->6 | 4->1 4->2 4->3 | 5->2 5->3 | 6->4 6->5
Non-edges read off from N_ab = 1 (b -/-> a), 11 of them:
 5-/->1 6-/->1 6-/->2 6-/->3 1-/->4 2-/->4 5-/->4 3-/->5 4-/->5 4-/->6 5-/->6
Structure: bidirected pairs {1,2},{1,3},{2,3},{2,5},{3,4} (so {1,2,3} is a bidirected triangle); the only pair with NO edge in either direction is {4,5}; the other 9 adjacent pairs are one-way: 1->5 1->6 2->6 3->6 4->1 4->2 5->3 6->4 6->5.

## 3. Parameters, and the real shape of the parameter set

### What the paper states
- delta = (29/25) eps, "for every 0 < eps < eps0 and every theta > 0". theta is arbitrary positive (the paper notes scaling x = theta y changes values, not supports or stability).
- beta = 54/25 and q = (29/25)(1 - eps) -> 29/25 as eps -> 0. The abstract: "Its values of q = delta(1 - eps)/eps tend to 29/25."
- eps0 is existential and NOT explicit. Proof ending: "Choose T larger than T_1, the stability threshold supplied by Lemma 3.1, and 25/4. For t > T, set eps = 1/(t+1), delta = (29/25) eps ... Finally, taking eps0 = 1/(T + 1) proves the uniform assertion." (T_1 comes from a nondegeneracy argument, Lemma 3.1's threshold from a Rouche argument; neither is computed.)
- So "q ~ 29/25" is a LIMIT as eps -> 0, not a single point. The pair (eps, q) matters, not q alone.

### What I computed myself (exact; NOT in the paper)
(a) Explicit threshold at beta = 54/25. Full support is permitted for every t (Lemma 2.1: permission of K does not depend on t; here det C = 179298949/9765625 > 0 and C^{-1}1 > 0). Stability of K = I + (54/25)N + tJ: characteristic polynomial coefficients det(sI+K) = s^6 + 6(t+1)s^5 + 3(1300t+2153)/625 s^4 + 4(154725t+162689)/15625 s^3 + 3(5288000t+5551469)/390625 s^2 + 6(36878825t+40171729)/9765625 s + (209141700t+179298949)/9765625.
   Hurwitz minors Delta_3, Delta_4, Delta_5 each have exactly one positive root (about 21.0, 46.4, 52.93); Delta_6 = a6 * Delta_5. Stable iff t > t*, where t* is the root of Delta_5 isolated in
   t* in [4031542479/76165654 , 935858668/17680649]  (about 52.9312),  i.e. eps* = 1/(t*+1) about 0.018542.
   Exact checks: at t = lower endpoint the Delta signs are [+,+,+,+,-,-]; at the upper endpoint all six are positive.
   Concrete exact decisions with the repo verifier (ctln.analyze_support, full support, theta=1, delta = 29 eps/25):
     eps=1/2, 1/4: NOT LEGAL (q = 0.58, 0.87) ; eps=1/8..1/20: legal, positive, UNSTABLE (Hurwitz minor 3 negative) ; eps=1/50: UNSTABLE (minor 5 negative) ;
     eps=1/100 (q=2871/2500=1.1484): STABLE, nonclique -> agrees=False (this is the counterexample) ; eps=1/1000 (q=28971/25000): STABLE ; eps=1/10^6: STABLE.
(b) eps -> 0 limit (Lemma 3.1 asymptotics), with b := beta = 1 + q. Stable full support as eps -> 0 iff the restriction D(b) of P(I + bN) to 1-perp is stable. I computed det(zI + D(b)) exactly and its Hurwitz minors:
     stable for 2 < b < b_c and unstable just above, where b_c = 2.16635502785... is the smallest root > 2 of the degree-6 polynomial
     797 b^6 - 6468 b^5 + 22544 b^4 - 44440 b^3 + 53280 b^2 - 38016 b + 13824   (sympy factors it as irreducible over Q; it is a common factor of Delta_4 and Delta_5).
     So this graph's eps->0 stable range is q in (1, q_c), q_c = b_c - 1 = 1.1663550..., an algebraic number of degree 6. Geneson's 29/25 = 1.16 sits just inside it.
     Permitted limit: (adj(C)1)_5 = -(b^3 - 3b^2 + 2b - 1); full support is permitted iff b < b_p = 2.3247179..., and b_p - 1 is the plastic number (root of x^3 = x + 1). All other components of adj(C)1 and det C are positive for all b > 2.
     Consequence: for q >= b_p - 1 = 1.3247 this graph is NOT even permitted at any eps (beta = 1 + q/(1-eps) > 1 + q). In particular it is invisible at q = 3/2 for every eps.
(c) Sampled table (ctln.analyze_support, full support; SAMPLED ONLY, not region-certified). eps = 1/k, delta chosen so q is exact. "." = positive but unstable, "pos-" = not a positive solution, S = stable nonclique:
   q \ k:     2    3    4    6    8   12   16   24   32   48   64  100  200  400 1000 1e4 1e6
   101/100  pos- pos- pos-  .    .    S    S    S    S    S    S    S    S    S    S    S    S
   21/20    pos- pos- pos-  .    .    .    S    S    S    S    S    S    S    S    S    S    S
   11/10    pos- pos- pos-  .    .    .    .    S    S    S    S    S    S    S    S    S    S
   29/25    pos- pos- pos- pos- pos-  .    .    .    .    .    .    .    .    S    S    S    S
   6/5      pos- pos- pos- pos- pos-  .    .    .    .    .    .    .    .    .    .    .    .
   5/4 ... 1.3: no S anywhere on the grid ; 3/2, 2, 3, 4: pos- everywhere on the grid.
(d) Nondegeneracy, exact at eps = 1/100 and 1/1000: all 255 determinants (63 principal minors of I-W and 192 uniform-input column-replacement determinants) are nonzero; the fixed point x* (theta=1) is strictly positive and proportional to the paper's vector v.

## 4. Which support, fixed point, and why it is not a target-free clique

- Support: the FULL vertex set {1,...,6}. Quote: "its full support is permitted, and the corresponding fixed point is stable. The full support is not a clique".
- Fixed point (formula): "x*_sigma = (theta/eps) K^{-1} 1_m, x*_{[n]\sigma} = 0." (13). For this graph the paper prints the proportional vector (eq. 20)
  v = (21597025, 34234375, 34234375, 46065775, 5850625, 67159525)^T, det C = 179298949/9765625, Cv = 179298949 * 1_6, C = I_6 + (54/25)N, "Thus C^{-1} 1_6 = v/179298949 > 0."
  and, with Lemma 2.1, K^{-1}1 = C^{-1}1 / (1 + t 1^T C^{-1} 1). The paper does not print x* itself. My exact x* at eps = 1/1000, theta = 1 is v * (1000/209111857249) componentwise (x_1 = 21597025000/209111857249, ...).
  Lemma 2.1 verbatim: "K is permitted <=> C is permitted. When these conditions hold, K^{-1} 1_d = C^{-1} 1_d / (1 + t 1_d^T C^{-1} 1_d)."
- Why it is not a target-free clique: the paper says only "Since N != 0, its support is not a clique." Concretely, 11 ordered pairs are non-edges; the pair {4,5} has no edge in either direction (N_45 = N_54 = 1). The full support has no outside vertices, so condition (14) is vacuous and "target" cannot arise; failure is purely non-cliqueness.
- Proof of stability (all exact): block reduction giving det B = 4,482,473,725 = 25 * 179,298,949 for B = 25C; restriction D = M/25 of P_6 C to 1-perp with M printed; det(xI_5 + M) = (x+25)(x^4 + x^3 + 4101x^2 + 3235x + 1394278); quartic Routh-Hurwitz (H1,H2,H3,H4) = (15625, 8457031250, 343562500000000, 479021635375000000000); Lemma 3.1 (Rouche) gives K = C + tJ_6 stable for all sufficiently large t; Jacobian at the fixed point is -eps K.
  EVERY ONE of these printed numbers reproduced exactly in my independent sympy computation (det B, det C, v, Cv, D = M/25, the factorization, the quartic, and H1..H4).

## 5. Positive result, verbatim

Abstract: "In the complementary direction, for any CTLN on n >= 3 vertices, we prove that in the parameter range
  q >= n - 2 - ((n-3)/2) eps,
no nonclique support can satisfy both the fixed-point positivity and linear stability conditions. Consequently, throughout this range, every nondegenerate CTLN has exactly its target-free cliques as supports of stable fixed points. In particular, this holds when eps <= delta/(delta + n - 2)."

"Theorem 4.1. Let delta > 0 and 0 < eps < delta/(1 + delta), and define t, beta, q by (5). Let m >= 3 and let N in {0,1}^{m x m} be nonzero with zero diagonal. If K = I_m + beta N + t J_m is permitted and stable, then
  q < m - 2 - ((m-3)/2) eps.   (24)
Consequently, at these parameters, every CTLN on a graph with n >= 3 vertices has the algebraic target-free clique classification whenever
  q >= n - 2 - ((n-3)/2) eps.   (25)
For nondegenerate CTLNs, the same condition gives the dynamical target-free clique classification."

"Corollary 4.2. Let n >= 3, delta > 0, and 0 < eps < delta/(1 + delta). Every CTLN at these parameters has the algebraic target-free clique classification, and every nondegenerate CTLN has the dynamical target-free clique classification, if
  eps <= 2 delta / ( delta + n - 2 + sqrt( delta^2 + 2 delta + (n-2)^2 ) ).   (33)
In particular, either conclusion holds if  eps <= delta/(delta + n - 2).   (34)"
Proof line: "Finally, (34) is equivalent to q >= n - 2, which implies (25)."

Secondhand statement CHECKED: "q >= n - 2 - (n-3) eps / 2" is exactly the paper's (25). CONFIRMED. The counterexample q ~ 29/25 is CONFIRMED as the eps -> 0 limit of q = (29/25)(1-eps).

Definitions the theorem relies on (verbatim): "A real square matrix is called stable if every eigenvalue has positive real part." and "A real d x d matrix A is permitted (for uniform input) if it is nonsingular and A^{-1} 1_d > 0 entrywise. We call sigma permitted when K is permitted." and "We say that sigma satisfies the algebraic stable-support conditions if it is permitted, K is stable, and all inequalities (14) hold." (14): "for every k not in sigma, theta + sum_{j in sigma} W_kj (x*_sigma)_j <= 0." Nondegeneracy: "A CTLN is nondegenerate if, for every nonempty subset rho of [n], the principal matrix (I_n - W)_rho is nonsingular and, for every i in rho, the determinant obtained by replacing its column indexed by i with 1_|rho| is nonzero [3]."

Structural remarks that matter for the project (my reading of the proof, not stated as such in the paper):
- Thm 4.1 never uses the outside inequalities (14). It bounds q using only the m x m matrix K of the support, so it applies to a support of size m inside any n. The bound is increasing in m ("The right side of (24) increases with m by 1 - eps/2 > 0"), so for support size m the safe range is q >= m - 2 - (m-3)eps/2. For m = 5: q >= 3 - eps. For m = 6: q >= 4 - 1.5 eps. For m = 3 the bound reads q < 1, which is impossible at legal parameters, so no nonclique size-3 support can be permitted and stable (paper also does m<=2 separately).
- Proof route is "not-permitted-and-stable" via coefficients c1, c2 of det(zI+K) only (inequalities (28)-(29), (31)-(32)); it drops a -beta^2 w term, so the bound is likely not tight.
- Old bound for comparison, verbatim: "the all-graph sufficient bound in [3] is q > n^2 - n - 1."

## 6. Other results, proof method, scope statements, future work

- Other items: Lemma 2.1 (permitted K iff permitted C; formula (15)); clique rule recomputed: "(I_n - W)_sigma = eps I_s + (1 - eps) J_s", input to outside vertex u missing r clique in-edges is theta[eps - r(eps + delta)]/(s - (s-1)eps), "positive for r = 0 and strictly negative for r >= 1"; an observation that a nonzero 0/1 N of order <= 2 cannot make K stable at legal parameters (det K = 1 + 2t - beta t L - beta^2 w); Lemma 3.1 (large-t stability from stability of (P_m C) on 1-perp, via Rouche).
- Method of counterexample: EXACT, not numerical. Integer determinants, polynomial factorization with integer coefficients, quartic Routh-Hurwitz, a mod-2 argument for nondegeneracy (B_rho = 25 I + 54 N_rho is = I mod 2), and an ASYMPTOTIC (existence) argument for large t. No explicit eps0 and no computer-algebra output is cited.
- Status claims about small cases (verbatim): "The conjecture that these are the only stable fixed-point supports of a nondegenerate CTLN [11] is known for symmetric graphs [4], oriented graphs [11], supports of size at most four, and several other graph and parameter classes [3]. It appears as Open Question 10 in [5]."
  The paper says NOTHING about n = 5, nor does it claim six is minimal. "Six-neuron" is in the title only. A 5-neuron full-support counterexample is NOT excluded by this paper. (It excludes only: support sizes <= 4 by citation, m=3 by Thm 4.1, and the high-q range.)
- Future work / open range: NO conclusion or future-work section. No statement of which q-range remains open and no mapping of the counterexample region. The only parameter statements are: counterexample at delta = 29 eps/25 for eps < eps0 (q just below 29/25), and the safe range (25). The gap between them is implicit.
- Related earlier work the paper mentions: "The analogous minimal-support conjecture for general competitive threshold-linear networks has recently been disproved [7], but that construction is not combinatorial."
- Acknowledgments, verbatim: "[AI tool names omitted here per project rule; see the paper] were used for proof exploration, proof criticism, exposition, and revision." The paper is a single-author arXiv v1 preprint; it has not been peer reviewed as far as I can tell. That is a reason the independent exact reproduction above matters.

## 7. Bibliography (verbatim reference list, relevant ones flagged)

[1] L. V. Ahlfors, Complex Analysis, 3rd ed., McGraw-Hill, 1979.
[2] C. Curto, J. Geneson, K. Morrison, Fixed points of competitive threshold-linear networks, Neural Computation 31 (2019), no. 1, 94-155, doi 10.1162/neco_a_01151.
[3] C. Curto, J. Geneson, K. Morrison, Stable fixed points of combinatorial threshold-linear networks, Advances in Applied Mathematics 154 (2024), 102652, doi 10.1016/j.aam.2023.102652.   <- source of the old q > n^2-n-1 bound, nondegeneracy, and "several other graph and parameter classes". HIGH PRIORITY to retrieve.
[4] C. Curto, K. Morrison, Pattern completion in symmetric threshold-linear networks, Neural Computation 28 (2016), no. 12, 2825-2852.
[5] C. Curto, K. Morrison, Graph rules for recurrent neural network dynamics, Notices AMS 70 (2023), no. 4, 536-551, doi 10.1090/noti2661; extended version, arXiv:2301.12638.   <- Open Question 10 lives here.
[6] F. R. Gantmacher, The Theory of Matrices, Vol. II, Chelsea, 1959.
[7] J. Geneson, Stable non-minimal fixed points of threshold-linear networks, preprint, arXiv:2511.05517 (2025).   <- 2025; disproves minimal-support conjecture for general competitive TLNs. Not yet retrieved by me.
[8] R. H. R. Hahnloser, H. S. Seung, J.-J. E. Slotine, Permitted and forbidden sets in symmetric threshold-linear networks, Neural Computation 15 (2003), no. 3, 621-638.
[9] R. A. Horn, C. R. Johnson, Matrix Analysis, 2nd ed., Cambridge University Press, 2013.
[10] J. Londono Alvarez, K. Morrison, C. Curto, Attractor-based models for sequences and pattern generation in neural circuits, Neural Computation 38 (2026), no. 3, 257-291, doi 10.1162/NECO.a.1492.   <- 2026
[11] K. Morrison, A. Degeratu, V. Itskov, C. Curto, Diversity of emergent dynamics in competitive threshold-linear networks, SIAM J. Appl. Dyn. Syst. 23 (2024), no. 1, 855-884, doi 10.1137/22M1541666.   <- source of the conjecture and of oriented-graph result.
[12] C. Parmelee, S. Moore, K. Morrison, C. Curto, Core motifs predict dynamic attractors in combinatorial threshold-linear networks, PLOS ONE 17 (2022), no. 3, e0264456.
Only [7] (2025) and [10] (2026) are recent. The paper cites no other 2025-2026 CTLN work.

## 8. What I could NOT verify

1. I did not read the TeX source or arXiv's own abstract-page abstract verbatim (the abstract-page fetch was paraphrased by the fetch tool). The abstract I quote is from the PDF text. Differences in typography between the PDF and the abstract page are possible but not expected.
2. The proof of Theorem 4.1 was read, not independently re-derived or machine-checked. I did not test it numerically on random graphs.
3. I verified every number printed in the proof of Theorem 3.2 (determinants, v, D, the factorization, Hurwitz numbers), the exact stability of the full support at eps = 1/100 and 1/1000 and 1/10^6 with the repo's own verifier, and nondegeneracy at eps = 1/100, 1/1000. I did NOT verify the Rouche-based Lemma 3.1 as a proof, nor the paper's specific T_1 / eps0 (not given).
4. My threshold t* ~ 52.93 (eps* ~ 0.01854), the limiting interval q in (1, 1.16636), and the permitted bound q < 1.3247 are my own exact computations, NOT claims of the paper. Table (c) is sampled and not a certified region.
5. The paper does not say whether any n <= 5 counterexample exists, and I did not check. I did not retrieve arXiv:2511.05517 or [3] (Curto-Geneson-Morrison 2024).
6. v2 or later versions: only v1 was listed at fetch time. I did not check for a later revision after that.
7. No scoop check was done in this task (only a WebSearch to confirm the ID/title). One search-engine snippet claims arXiv:2511.05517 states the conjecture was proven for "all sizes n <= 4"; that is an unverified secondhand paraphrase by the search layer, not something I read.

## 9. Convention notes for the project (not changed, flagged)

- CLAUDE.md section 2 and the ctln.py docstring write the closed form as I - W_sigma = -delta I + (1+delta) J - (eps+delta) A_sigma. With A[i][j]=1 iff i->j, this literal reading is WRONG by a transpose: the entry is (1+delta) - (eps+delta) A[j][i]. The code is right: ctln.identity_check_IminusW reads A[sj][si] ("A_sigma[i][j] uses the same transposed reading as W"), i.e. A_sigma there means the TRANSPOSED restriction. Suggest rewording the docs to "A_sigma^T" or define A_sigma[i][j] := A[sigma_j][sigma_i]. Paper-equivalent: I - W = eps K = eps (I + beta N + tJ) with N = J - I - A^T.
- Suggested ANOMALIES.md lines (I did not edit the file): (i) the Geneson graph's eps->0 critical q_c = 1.16635502785... has algebraic degree 6 (irreducible sextic above); (ii) its permitted threshold is the plastic number 1.3247 (root of x^3 = x + 1); (iii) the counterexample needs eps < ~0.0185, and a sweep at eps = 1/2 (CLAUDE.md section 5) can never see it, even at q in (1, 1.17), because the graph is not legal/stable there.

## 10. Reproduction files (scratchpad)

- /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/pdfwork/geneson.pdf and geneson.txt
- .../scratchpad/verify_geneson.py (N -> A conversion, paper certificates), verify_geneson2.py (verifier runs, Hurwitz minors in t), verify_geneson3.py (tight t*, eps->0 limit, sampled table), verify_geneson4.py (sign structure and factorizations), verify_geneson5.py (nondegeneracy). Run each with `python3 -I`.
