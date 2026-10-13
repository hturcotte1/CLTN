# Anomaly log

Append-only. Do not resolve entries by deleting them. Anything surprising goes
here: a case far slower than its neighbours, a critical value of unexpected
algebraic degree, a count that did not match, a reduction behaving oddly.

Format: `## YYYY-MM-DD  short title` then what happened, what you expected, and
what you did. Leave it open if unresolved.

---

## 2026-10-12  Boundary cases initially misreported as counterexamples
Setup sweep flagged 2 "violations" that were actually undecided boundary cases
(vanishing Hurwitz minor) being collapsed into the violation bucket. Fixed by
making `agrees` three-valued in `ctln.py`. Noted here because it is exactly the
failure mode that would waste days of an unattended run. Resolved.

## 2026-10-12  Irrational critical q values are common
`q_region` on a mixed n=4 graph returned a critical value of `sqrt(41)/2 - 2`;
an n=5 case returned `sqrt(2261)/14 - 2`. Implication: rational sampling of q
will systematically miss thin stable intervals. This is why Stage B exists.
Open — worth tracking whether algebraic degree correlates with anything.

## 2026-10-12  Undecided rate measured at 2.8% (n=6, eps=1/2, q=3/2)
2,000 random n=6 graphs, full support. Expect ~43,000 escalations over the full
1,540,944. Recheck this rate at other q values before sizing Stage B.

## 2026-10-13  Stage 0: Geneson's counterexample lives only at small eps, q near 1
arXiv:2607.21396 Thm 3.2 gives the graph as a receiver-first NONEDGE matrix N
(N[a][b]=1 iff b -/-> a); our A[b][a] = 1 - N[a][b] (complement AND transpose).
Parameters delta = 29 eps/25 for all sufficiently small eps (eps0 existential).
Verifier: positive FP for every eps tried; UNSTABLE at eps = 1/10, 1/20, 1/50,
STABLE at eps = 1/100, 1/1000. Literature agent's exact computation: stable iff
t = (1-eps)/eps > t* in [4031542479/76165654, 935858668/17680649] (~52.93),
i.e. eps < ~0.01854. In the eps -> 0 limit with b = 1 + q: stable iff
2 < b < b_c, b_c the smallest root > 2 of
797b^6 - 6468b^5 + 22544b^4 - 44440b^3 + 53280b^2 - 38016b + 13824 (irreducible,
degree 6, b_c ~ 2.16636); permitted iff b < b_p where b_p - 1 is the plastic
number (root of x^3 = x + 1). Implications: (a) the parameter plane is genuinely
2-D and a sweep at eps = 1/2 CANNOT see this example (it is not even a legal FP
there) -- the setup sweep's "1,495 pairs, zero counterexamples" is explained;
(b) critical values of degree 6 occur already at n = 6; (c) the sweep must be
organised around small eps (large t) and q in (1, ~1.2). Open.

## 2026-10-13  Rank-one structure: the eps -> 0 limit is a one-parameter problem
I - W = eps (C + t J) with C = I + beta N, t = (1-eps)/eps, beta = 1 + delta/eps
(paper's eq. 6; verified entrywise here on 200 random cases). Consequences,
all verified exactly: adj(C + tJ) 1 = adj(C) 1 (independent of t); det(C+tJ)
= det C + t 1^T adj(C) 1 (linear in t); charpoly(C+tJ) = chi_C - t psi_C with
psi_C(lam) = 1^T adj(lam I - C) 1 (linear in t). So positivity of the full-
support fixed point depends on beta alone (up to the sign of det), and the
t -> infinity (eps -> 0) stability condition is "roots of psi_C(beta; lam) in
the open right half-plane" -- a 1-parameter condition in beta = 1 + q. This is
the regime Geneson's example lives in. Plan: sweep the limit problem first.

## 2026-10-13  Documentation wording: closed form uses the transposed restriction
CLAUDE.md §2 and ctln.py write I - W_s = -d I + (1+d) J - (e+d) A_s. With
A[i][j]=1 iff i->j the entry is (1+d) - (e+d) A[j][i], so "A_s" there means the
TRANSPOSED restriction. The code (identity_check_IminusW, symbolic.py) reads
A[sj][si] and is correct; only the prose is ambiguous. Not changing §2 (fixed
by policy); noting here so nobody "fixes" the code to match the prose.

## 2026-10-13  symbolic.q_region used a numerical heuristic (nsimplify) and float probes
Cross-validating fastregion against symbolic.q_region on n=4 raised TypeError in
16 of 1,800 cases: `sp.nsimplify(root)` returned garbage closed forms such as
2**(52/963)*3**(830/963)*5**(40/963)*7**(503/963)/2 for an algebraic root, and
the probe point for each cell was a float midpoint pushed through
limit_denominator (a thin cell could be probed outside itself). Both violate
§3 (no floats in a decision path). Fixed: roots kept exact (CRootOf/radicals),
probe points suggested numerically but verified exactly (is_positive on the
differences) before use. fastregion.py never had this issue (VAS isolation
with rational endpoints). Resolved; the 1,784 cases that did run agreed.

## 2026-10-13  n=6 limit-mode sample: 96% of graphs die at positivity, all dominated ones do
3,000-graph sample (limit mode): 2,884 never_positivity, 116 reach the
stability stage, 0 stable, 0 escalations. Every one of the 116 is domination-
free; every dominated graph failed positivity -- consistent with Lemma D
(Curto Cor. 3.6). 13.3% of the 1,540,944 graphs are domination-free.

## 2026-10-13  n=5 full support: no counterexample in the limit or at eps = 1/100, 1/10, 1/4
run_region_sweep over all 9,608 digraphs on 5 nodes: limit mode and the three
slices each give exactly one stable_somewhere graph, the complete graph K5;
8,902 (8,916 at eps=1/4) die at positivity, 705 (691) reach the stability stage
and are never stable; 0 escalations; 960 analyze_support cross-checks per
slice all agree. Consistent with the two unrefereed GitHub claims that n <= 5
is clean. The two-parameter CAD on the 706 survivors is running; if all are
'never', n=5 full support is settled for ALL legal parameters (modulo the
lemmas in notes/lemmas.md), which by the reduction lemma (notes/lemma_F_
reduction.md) also settles every size-5 support inside n=6 graphs.

## 2026-10-13  Delta_5(beta, t) factors as (linear in t) x (quartic in t) for Geneson's graph
region2d on the Geneson graph: Delta_5 has two irreducible factors over Z with
t-degrees 1 and 4 (projection names lc_Delta_f0, lc_Delta_f1, disc_Delta_f1).
The 2-D stability region is {2 < beta < beta_c, t > max(1/(beta-2), T(beta))}
with beta_c the sextic root and T the quartic-factor curve; for beta below a
degree-14 algebraic number in (17/8, 15/7) the graph is stable at EVERY legal
t. Matches the structure claimed by github.com/ipitchford/ctln-six-is-minimal
(C6/C7). Whether the (1)x(4) factorisation is general is not yet known --
the 2D sweep records the t-degrees of the factors; check the distribution.

## 2026-10-13  n=5 full support settled for ALL legal parameters (two-parameter CAD)
run_region2d_sweep on the 706 positivity survivors among the 9,608 digraphs on
5 nodes: 705 'never', 1 'stable_somewhere' (K5), 0 escalations, 0 errors,
239 s. Every cell verdict was cross-checked against ctln.analyze_support.
So no 5-node non-complete graph has a positive, linearly stable full-support
fixed point at any legal (eps, delta) -- the claim of the two unrefereed
GitHub repos, now with an independent exact method. By the reduction lemma
(notes/lemma_F_reduction.md, proof pending) this also disposes of every
size-5 support inside every 6-node graph, so CLAUDE.md's Stage 4 is implied by
Stage 2/3 at n = 5 and n = 6. Conditional on notes/lemmas.md (Lemmas A-E).

## 2026-10-13  For Geneson's graph the counterexample region reaches q -> 1+: there is no q* > 1
Our two-parameter CAD (region2d.py) gives, for the Geneson graph, a first
beta-cell (2, beta_1) with beta_1 in (17/8, 15/7) a degree-14 algebraic number,
on which the full support is stable at EVERY legal t. Since q = (beta-1) t/(t+1)
-> 1+ as t -> 1/(beta-2)+, counterexamples exist arbitrarily close to q = 1.
Verified pointwise: eps = 1/12, delta = 11/120 (q = 121/120; the rational
witness of github.com/ipitchford/ctln-six-is-minimal) and beta = 2001/1000,
t = 2000 (q = 1.0009995) are both stable positive full-support fixed points.
The literature agent confirms that repo's seed is identical to Geneson's graph
(only the identity permutation maps one to the other), that its sextic is ours
under u = beta - 2, and that its Theorem 5.1 region {0<u<u*, t>max(1/u,T(u))}
has the shape of our cell decomposition. Consequence for Stage 5: the
interesting quantity is not a lower q* but the UPPER reach of counterexamples
(1.1664 for this graph vs Geneson's bound 4 - 1.5 eps at n = 6).

## 2026-10-13  Prover corrections (notes/lemmas.md, Lemmas A-E)
(a) Lemma D (dominated node => x_j = 0 at every fixed point) does NOT need
eps < 1, only eps > 0, delta >= 0 -- but it DOES need uniform theta: with
theta_j = 10, theta_k = 1 on the single edge j -> k there is a stable fixed
point with x_j = 10. Our CTLNs have uniform theta, so the prune stands.
(b) Lemma A.2: a stable K has det K > 0, so the fixed-point positivity
condition is exactly "all a_i > 0" (not merely a common sign) and the
off-support condition is "o_k < 0". fastregion's stage-1 filter is tightened
accordingly (sound: it only removes graphs that could never be stable).
(c) Reporting gap in fastregion.region: an escalated cell made the verdict
'escalate' and suppressed stable_intervals even when another cell was decided
stable. No escalation has occurred in any run so far (n=5 complete, n=6
sample, n=6 limit sweep in progress), so no result is affected; fixed so that
a decided stable cell is always a hit.
(d) Geneson's degree-6 critical value is the polynomial whose roots are
-2/(nu_i + nu_j) over pairs of roots of the quartic factor of the
characteristic polynomial of the compressed nonedge matrix -- at most
C(4,2) = 6, explaining the degree.
(e) Independent replay of sympy's root isolation added (Sturm counts and
sign changes at every isolating interval); enabled in both sweep drivers.
The running n=6 limit sweep started before this; a replay pass over all
1,540,944 graphs is scheduled after it.

## 2026-10-13  Reduction lemma proved (notes/lemma_F_reduction.md); two corrections
Theorem F: at every legal parameter point, the target-free clique conjecture
holds for all n-node graphs iff for every 2 <= k <= n no non-complete k-node
graph has a full support satisfying conditions (1) and (3). Correction to my
sketch: condition (2) IS needed in the clique branch (to exclude targets: the
graph 0 -> 1 with sigma = {0} satisfies (1),(3) but node 1 is a target); it is
never needed in the non-clique branch. Consequence: for n = 6 the open content
is exactly FS(5) and FS(6) (full-support questions), and FS(5) is now certified
over the whole legal plane (region2d), so there is NO size-5 counterexample
support inside any 6-node graph; Stage 4 is empty. Also proved: Lemma O (the
stable-support parameter set is open, so every counterexample has a rational
witness) and Lemma H (positivity and condition (2) depend on beta alone).

## 2026-10-13  The setup's "2.8% undecided" rate was an artifact of the sample point (1/2, 3/2)
Prover F: on the same 3,000 random 6-node graphs, (eps, delta) = (1/2, 3/2)
gives 59 undecided cases (in all of them the cubic coefficient of the
characteristic polynomial vanishes), while four other legal points, including
delta = 3/2 + 1/997, give 0 of 3,000. 56 of the 59 already fail positivity.
Furthermore a genuine counterexample point is never 'undecided' when evaluated
there (all Hurwitz minors strictly positive, Lemma O): counterexamples live in
open regions BETWEEN sample points, not in the undecided pile. CLAUDE.md §5's
sizing of Stage B from the 2.8% figure should not be used; the exact region
machinery supersedes the sample-and-escalate design.

## 2026-10-13  region2d raised on graphs whose boundary curves share a factor
Prover F found ArithmeticError("common factor") on three dominated 4-node
graphs and on all 2-node graphs (det K equals minus the legality polynomial
for a one-way edge). Fixed: the curve set is made pairwise coprime by splitting
off gcds before taking resultants; all 2- and 3-node graphs now give the
expected verdicts; n = 4 check running.

## 2026-10-13  Spectral form of the limit proved (notes/lemma_I_spectral.md, I.1-I.7)
Proved: limit-stable beta-set is exactly (2, beta_c), beta_c = 1/(1+alpha),
alpha = spectral abscissa of the adjacency compressed to 1^perp; positivity
at EVERY t is "w = (A^T + gamma I)^{-1} 1 > 0" with gamma = 1 - 1/beta (plus a
sign of det(A^T + gamma I) that condition (3) supplies); |E| >= 16 at n=6 and
>= 22 at n=7 are necessary for limit stability; the Bendixson sufficient
condition holds only for complete graphs (useless), and tournaments have every
compressed eigenvalue at real part exactly -1/2.
New Lemma I.7 (not in the literature as far as we know): every root of
gcd(chi_A, chi_{A+J}) yields an eigenvalue 1 - beta - beta*mu of K at EVERY
(beta, t); if one has Re mu >= -1/2, condition (3) fails at every legal
parameter. Graph patterns: two non-adjacent nodes with the same in-neighbours
(or the same out-neighbours); a module that every outside node treats
uniformly and that induces a directed 3-cycle. Kills 896 of the 1,433
domination-free n=5 graphs and 211 of 400 in a 3,000-graph n=6 sample
(exploration counts). NOT used as a production prune: the two-parameter CAD
is complete without it; kept as a structural lead for Stage 5.
Open: alpha(A - tau J) is provably NOT monotone in tau (3-node example), so
"the stable set in t at fixed beta is an interval" is not established; the
CAD handles finite t exactly, so no result depends on it.
At n=5, 738 limit-mode verdicts are 'undecided' in the prover's exploration
(167 with alpha = -1/2 exactly, incl. all 12 tournaments): these graphs are
limit-boundary cases; our fastregion limit sweep never sampled a boundary
(it samples inside cells), and the 2D CAD decided all of them 'never'.

## 2026-10-13  RUN STOPPED: second stable non-clique full support at n = 6 (new graph)
The limit-mode sweep (all 1,540,944 digraphs on 6 nodes, three shards) produced
two hits. Hit 1 (idx 1000255, key 1002304259) is Geneson's graph, as expected.
Hit 2 (idx 1387528, key 1066791683) is a NEW graph, not isomorphic to Geneson's:
    0110 00 / 000011 / 000111 / 110010 / 100101 / 111110   (rows = who i sends to)
18 edges, 4 bidirected pairs, 10 one-way edges, one mutually non-adjacent pair
{1,2}, trivial automorphism group, domination-free, no target (full support).
Limit certificate (eps -> 0): stable iff 1 < q < q_c, q_c the root of
q^5 + 2q^4 + q^3 - 2q^2 - 2q - 1 in (1, 13/12)  (DEGREE 5; Geneson's is degree 6).
Compressed-adjacency charpoly 6 lam^5 + 18 lam^4 + 30 lam^3 + 29 lam^2 + 16 lam + 4.
Pointwise exact witnesses (ctln.analyze_support, beta = 53/26): t = 100, 1000,
10^4, 10^5 (eps = 1/101 ... 1/100001, delta = (27/26) eps) all legal, positive,
all Hurwitz minors strictly positive, agrees = False. E.g. eps = 1/101,
delta = 27/2626, q = 1350/1313.
Per CLAUDE.md §3 the run is STOPPED here pending human hand-verification:
certificates/COUNTEREXAMPLE_CANDIDATE_n6_idx1387528.json. The sweep had
processed about 1.45M of 1.54M graphs (shard files intact, resumable with
run_region_sweep.py --resume). Nothing has been built on the candidate.

## 2026-10-13  Third limit-sweep hit (idx 1424629) is cospectral with Geneson's graph but not isomorphic
After the human restarted the run, the sweep flagged idx 1424629 (key
1069347459): 19 edges, 5 bidirected pairs, one mutually non-adjacent pair,
trivial automorphism group, domination-free, NOT isomorphic to Geneson's graph
and not its reversal (reversal keys of all three hits are not hits: reversing
transposes W, preserving the stability spectrum but not positivity). Yet its
limit boundary is the SAME sextic 797q^6 - 1686q^5 + ... and its compressed
adjacency characteristic polynomial is the same 6l^5+19l^4+32l^3+32l^2+17l+4.
Pointwise witnesses at beta = 53/26 (t = 100..10^5) are stable, positive,
non-clique. Certificate: certificates/COUNTEREXAMPLE_CANDIDATE_n6_idx1424629.json.
Lead for Stage 5: counterexamples may come in cospectral families (a
switching-type operation on the nonedge set that preserves A|perp's spectrum).
Hit 2 (idx 1387528, 18 edges, degree-5 critical value) is NOT cospectral with it.

## 2026-10-13  Two-parameter certification was slow on dense survivors; a wrong Lienard-Chipart shortcut was caught
Profiling: on an 18-edge survivor region2d took 96 s, 94 s of which were my
Sturm-sequence replay in Python Fractions on degree-100 projection polynomials
(sympy's own Sturm counter is equally slow). Replaced by a Descartes-bisection
counter on the squarefree part (integer arithmetic): same graph now ~2-4 s.
While rewriting I tried to drop Delta_5 in favour of the even Hurwitz minors
(Delta_2, Delta_4) via Lienard-Chipart. WRONG: the halved minor sets always
include the top minor (Delta_6 = a_6 Delta_5 for degree 6), and the Geneson
graph immediately showed the discrepancy (a stable cell disappeared because the
limit boundary beta_c is a root of lc_t(Delta_5), not of lc_t(Delta_4)).
Reverted to the minimal complete set {a_i, o_k, D, Delta_{m-1}, L} of Lemma B.
Lesson: every change to the boundary set must be checked against the Geneson
cells before any sweep. v2 also restricts the decomposition to the positivity
cells (everything outside is 'never' by Lemma A.2) and isolates jointly.

## 2026-10-13  Candidate idx 1387528 independently verified (workflow, 7 agents)
Verdict VERIFIED by three from-scratch exact re-derivations (Routh array,
Lyapunov/Sylvester, Fraction-only Hurwitz), two skeptics (conventions against
both papers' verbatim definitions; margins, exact zeros, Geneson nondegeneracy:
63 + 192 determinants nonzero) and a context agent (not isomorphic to Geneson's
graph or its reverse; not in the paper or the ipitchford repo). Stability
margin at P1 (eps = 1/101): max Re(lambda) < -25914987/68719476736, i.e. about
-3.8e-4 = 3.8% of eps, structural not coincidental. Transposing the adjacency
HIDES the example (x gets three negative entries), so no convention error can
manufacture it. Full report: notes/candidate_1387528_verification_workflow.json.

## 2026-10-13  Two non-isomorphic 6-node graphs with IDENTICAL full-support data at every (beta, t)
Geneson's graph G (idx 1000255) and hit 3 H (idx 1424629) are not isomorphic
(out-degree sequences [2,2,3,4,4,4] vs [2,2,3,3,4,5]; brute force over all
720 permutations, also against G^T), yet beta_data(G) == beta_data(H) up to a
permutation of the a_i: same det C(beta), same phi (charpoly of -C), same omega
(so I - W has the same characteristic polynomial for ALL beta, t), same
multiset of adjugate column sums. Consequently their positivity, stability
and whole 2-D regions coincide exactly, and both slices (eps = 1/100, 1/10)
returned identical q-windows with identical minimal polynomials. Both have 19
edges, 5 two-cycles, 12 three-cycles, the same in-degree sequence. This is a
non-trivial "CTLN full-support equivalence" between non-isomorphic graphs; the
transformation relating them is not yet identified (not reversal). Lead for
Stage 5: count how many such equivalence classes the 34,071 survivors fall
into, and characterise the operation.
Slice facts: at eps = 1/10 only G and H have a stable window ((1, root in
(81/80, 36/35))); hit 2 (idx 1387528) is stable at eps = 1/100 for
1 < q < root in (693/650, 429/400) but has no window at eps = 1/10.

## 2026-10-13  Full-support equivalence classes among the 34,071 n=6 survivors
Grouping survivors by (det C, phi, omega, sorted adjugate sums) gives 21,800
classes: 17,147 singletons, 2,994 pairs, ..., and classes of size up to 259.
Geneson's class is exactly {idx 1000255, idx 1424629}; hit 2 (idx 1387528)
is a singleton. All members of a class have identical 2-D regions, which gives
a consistency check on the 2-D sweep (verdicts must agree within a class) and
a 36% saving had the sweep run on representatives. The operation generating
these classes is unidentified (certificates/n6_survivor_equivalence_classes.json).

## 2026-10-13  Candidate idx 1424629 verified; the equivalence operation is a "twin move"
Second adversarial workflow (7 agents): VERIFIED at P1 (eps = 1/101, delta =
27/2626) and P2 (eps = 1/1001), certified max Re(lambda) <= -498/872525 and
<= -23/390905, no zero minor, nondegenerate over all 63 supports. The context
agent found the relation to Geneson's graph G: after relabelling, H = G minus
the edge 3->1 plus the edge 5->1, where nodes 3 and 5 are IN-TWINS (identical
in-neighbourhoods {1,2,4}, mutually bidirected). Moving out-edges between
in-twins preserves chi_A, chi_{A+J}, the full-support characteristic
polynomial AND the fixed point identically in (eps, delta, theta) (machine-
checked on this pair; the general lemma is with prover J). Consequences:
(a) x_3 = x_5 exactly and -eps is an exact eigenvalue of -I + W (e_3 - e_5 is
a left eigenvector), which is why the slices showed a factor (101 lam + 1);
(b) H has no target-free clique at all, and its full support is its only
stable fixed point at P1, P2 (scan_supports); (c) H carries no new region
information relative to G. Open: whether twin moves generate all 21,800
equivalence classes among the survivors.
Side notes from the verifiers: ctln.hurwitz_stable builds the transpose of
the textbook Hurwitz matrix (harmless: leading principal minors are
transpose-invariant; do not "fix"); a sympy pitfall (Symbol with assumptions
can make subs/factor a silent no-op) -- symbolic.py uses the same Symbol
objects throughout, so it is not affected, and it was cross-validated against
fastregion on 8,094 cases.

## 2026-10-13  Lemma J (notes/lemma_J_equivalence.md): the full-support equivalence is characterised
Theorem J.2 (proved): D(A) = D(B) iff, after relabelling, A and B have the same
characteristic polynomial and the same in-walk counts (A^T)^k 1 for k = 0..n-1.
This reproduces our 21,800 classes exactly; the walk sequence cannot be
truncated (length <= 4 merges 54 classes). Rank-one moves N' = N + y x^T
preserve D iff x is orthogonal to all N^k 1 and N^k y (J.1); bidirected
in-twin switching is the special case x = e_u - e_v. G2 = idx 1424629 is
Geneson's G1 with the edge 1->4 replaced by 2->4 after relabelling (1,2
bidirected in-twins). J.4 (proved): bidirected in-twin switching preserves
EVERY fixed point of EVERY support at EVERY (eps, delta, theta) -- fixed-point
status, vector, Jacobian characteristic polynomial, target-free cliques. So
counterexamples should be reported modulo bidirected twin switching; hit 3
is not independent of Geneson's graph. Other D-equivalences are NOT CTLN-
natural (exact counterexamples with different FP sets inside one class).
Rank-one moves do not generate all classes (95 residual classes need rank 2-3;
the 259-member class is plain cospectrality with all in-degrees 3).
J.4' (proved): two in-twins inside a support give -I + W_sigma the exact
eigenvalue -eps (bidirected; explains the (101 lam + 1) factor) or +delta
(non-adjacent) -- a support containing a NON-ADJACENT in-twin pair is never
stable: cheap sound prune for Stage 4 / n = 7. Not yet checked against the
literature.

## 2026-10-13  certificates/n6_survivors.idx mixes two positivity filters
The limit sweep's first ~318,000 indices per shard ran before commit b68c5af
(common-sign filter); the rest used the tightened all-a_i > 0 filter. Both are
sound, but the 34,071-graph survivor set is a superset of the 28,350 strict
survivors and is not closed under D-equivalence (17,307 strict classes vs
21,800). The 2-D certification runs on the superset, so nothing is lost; the
replay pass with the final code (task 11) will produce the reproducible set.
