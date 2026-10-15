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

## 2026-10-13  Literature status of the twin lemmas; a sign typo in CGM24 Thm 2.5
J.4' is a special case of Curto-Geneson-Morrison's simply-added Lemma 3.13 /
Cor 2.7 (which is a stronger prune: any simply-added uniform-in-degree tau
that is not a stable motif kills sigma). J.4 is a short corollary of their
machinery; its full statement (all supports, Jacobian charpoly, stability)
was not found in print. J.2 (charpoly + in-walk counts) was not found.
CGM24 Theorem 2.5 prints the clique eigenvalue as |sigma|(1-eps) - eps; the
trace forces |sigma|(1-eps) + eps (correct in CGM19a). Our test 3 uses the
correct value. Scoop re-check: no change. See notes/literature_twins_check.md.

## 2026-10-13  n = 6 full-support map COMPLETE
Two-parameter certification of all 34,071 positivity survivors: 34,067 never,
4 stable somewhere (K6 + idx 1000255, 1424629, 1387528), 0 escalations,
0 errors; consistent with the limit sweep and all five eps-slices. With the
n = 5 result and Theorem F this settles the conjecture at n = 6: it fails
exactly for the two twin-equivalence classes above, only for full supports,
only for q below ~1.166 (resp. ~1.079), and for every q > 1 close to 1.
Hit 2's region is a single cell: stable iff 2 < beta < beta_1' at EVERY legal
t (so its t -> infinity boundary and all-t boundary coincide), unlike
Geneson's class where a quartic curve T(beta) separates them.

## 2026-10-13  Replay pass with the final code: 28,350 strict survivors, same four stable graphs
run_region_sweep limit mode over all 1,540,944 graphs with the final code
(tightened positivity filter, Sturm replay): 1,512,594 never_positivity,
28,350 survivors (a subset of the earlier 34,071; 5,721 removed), 4 stable
(K6 + the 3 counterexamples), 0 escalations, 22,666 s. The 2-D certification
covered the superset, so every strict survivor is certified.

## 2026-10-13  Partial adversarial audit (3 of 6 agents completed before a usage limit; resumed)
CAD-soundness audit found no verdict error: projection set is the Collins/
McCallum set; delineability verified on 9,300 random points in 27 graphs;
1,308 random points outside the positivity cells all non-stable; Descartes
replay verified on 800 constructed and real polynomials. Findings acted on:
(1) 298 records written by the first-generation engine before the container
restart were kept by --resume; the auditor re-ran them with the final engine
(0 verdict changes); they are now re-certified in place. (2) The first
isolating interval above a hull end is now refined explicitly instead of
relying on the legality factor's degenerate interval to force separation.
(3) The aggregator now checks that stable records contain no undecided cell
(the four stable records have none). (4) Documentation: beta_critical in the
records lists all roots of the joint isolation on the positivity hull,
including the degenerate root 2 and projection roots outside the positivity
cell; the two-parameter form of Lemma B (constancy of the verdict on
connected components of {D != 0, Delta != 0}) plus openness (Lemma O) is what
the CAD uses and should be written out in notes/lemmas.md. (5) The replay is
independent of the VAS isolator but uses sympy's squarefree part.

## 2026-10-13  Pre-restart 2-D records re-certified with the final engine
The 298 records (79 + 132 + 87) that the first-generation region2d wrote
before the container restart and that --resume had kept were recomputed with
the final engine in place: 0 verdict changes. certificates/region2d_n6_shard*
are now uniformly produced by the committed engine; n6_summary regenerated
(34,067 never / 4 stable / 0 escalations / 0 errors; no stable record has an
undecided cell).

## 2026-10-13  Lemma B' (two-parameter sign-invariance) written; sympy resultant sign quirk
notes/lemma_Bprime_2d.md: Lemma B' (verdict constant on connected subsets of
the legal region avoiding the zero sets of a_i, o_k, D, Delta_{m-1}),
Corollary B'' (one sample per sign-invariant cell; undecided samples; only
full-dimensional cells needed, via openness), Lemma B''' (elementary
delineability for the projection set {lc, disc, pairwise resultants} with
vertical polynomials), Lemma N (no boundary polynomial is identically zero:
all equal 1 or are positive at (beta, t) = (0, 0)). Findings: sympy 1.14
resultant(f, g) can return MINUS the Sylvester determinant when deg f < deg g
and deg f * deg g is odd (resultant(t-1, t^3-8, t) = 7); harmless here (only
zero sets are used) but never rely on its sign. The stable-cell tuples in the
map are slices at the sample beta, not boxes (presentation fixed).

## 2026-10-14  Full audit verdict and the follow-up checks
Resumed audit workflow (6 agents): "major gaps, none blocking, no verdict
error". Independent from-scratch symbolic re-derivation of the idx 1387528
region reproduced {2 < beta < beta_1', every legal t} with the same minimal
polynomial beta^5 - 3beta^4 + 3beta^3 - 3beta^2 + 2beta - 1 (isolating
interval (152/73, 25/12), beta_1' ~ 2.08302; an earlier report draft wrote
(152/73, 27/13), wrong, fixed). Fixes made: every HITS file is now tracked
(the broad .gitignore pattern had excluded them), gzipped bulk records and a
sha256 manifest committed; the aggregator uses the final-code replay run as
the limit run (1,512,594 / 28,350) and reports the 34,071 as the certified
superset; the report states the filter semantics correctly (no legal beta
with all a_i > 0 does not exclude positive UNSTABLE fixed points -- the empty
graph has one). Follow-up checks, no model in the loop: (i) the projection
algebra of region2d (Delta_{m-1}, 1,900+ resultants, discriminants, covering)
re-verified by independent Bareiss determinants on 321 cases
(notes/projection_check_321_cases.log); (ii) adversarial spot-check of the
151 'never' graphs with the widest positivity cells or >= 17 edges: 22,650
exact random legal points concentrated near beta -> 2+, t -> t_min+ and large
t, 0 stable points (certificates/n6_never_spotcheck.json). Remaining trust
base: sympy 1.14 exact polynomial algebra, Routh-Hurwitz, Orlando; the
never-side records are regenerable by the committed engine rather than
self-contained.

## 2026-10-14  n = 5 and n <= 4 artifacts regenerated with the committed engine
n = 5 limit sweep (final filter): 9,218 never_positivity, 390 strict
survivors, K5 the only limit-stable graph; 2-D certification of the 390:
389 never + K5, 0 escalations. All 237 graphs on 2-4 nodes: only the complete
graphs are ever stable. Identical verdicts to the superseded-code runs (706
superset: 705 never + K5). The k <= 4 literature theorem is now corroborated,
not depended on, by the committed engine.

## 2026-10-14  Human sign-off on the two new n = 6 counterexamples
The project owner hand-verified certificates/COUNTEREXAMPLE_CANDIDATE_n6_
idx1387528.json and _idx1424629.json and signed them off; the CLAUDE.md §3
gate is closed and the n = 6 theorem stands as stated in notes/REPORT_n6.md.
Stage 6 (n = 7) starts now.

## 2026-10-14  Stage 6 (n = 7) tooling: nauty + exact C filter, validated at n = 6
nauty 2.8.8 (from the pynauty sdist) built; geng | directg reproduces
A000273 for n = 3..6 (16, 218, 9608, 1540944) and showg confirms the digraph6
bit (i,j) = 1 iff i -> j convention on Geneson's graph (all 19 arcs). filter7.c
applies two PROVEN necessary conditions for limit stability of the full
support: domination-free (Lemma D) and alpha(A|perp) < -1/2 (Lemma I.3/I.5),
the latter as an exact __int128 Hurwitz test of 2^(m-1) m chi_perp(mu - 1/2);
a vanishing minor is kept as 'B'. At n = 6: domination-free 204,378 (exactly
our Python count; reversal-sensitive), 3,594 survivors (28 'S', 3,566 'B'),
the four limit-stable graphs all 'S'; verdicts agree with an independent
Python implementation on all survivors plus 528 random domination-free
graphs (0 mismatches once the same coefficient rule is used; ctln.hurwitz_
stable is over-conservative when a zero and a negative coefficient coexist).
Since limit stability at some legal beta <=> alpha < -1/2 strictly, 'B'
graphs are certified never limit-stable; only 'S' graphs go to the exact
limit stage. Whole n = 6 pipeline: 2.6 s; n = 7 (880,471,142 graphs) runs in
four geng res/mod shards.
SCOPE of Stage 6: the eps -> 0 (limit) classification at n = 7 is complete and
certified by this route; the finite-eps two-parameter map is NOT (the
positivity survivors would number in the millions and the CAD costs ~1 s
each). Finite-eps counterexamples invisible in the limit cannot be excluded
at n = 7 by this stage.

## 2026-10-14  n = 7 enumeration: 882,033,440 digraphs (OEIS A000273(7)); an earlier figure was misquoted
geng 7 r/4 | directg | filter7 over four shards: total 882,033,440, which IS
A000273(7) (the value 880,471,142 quoted earlier in this log and in the
report draft came from memory and was wrong; CLAUDE.md lists the sequence
only through n = 6). graphs.A000273 extended with the n = 7 value.
Domination-free: 117,232,922 (13.3%, same fraction as n = 6). Spectral
survivors: 638,042, of which 1,889 strict (alpha < -1/2) and 636,153 boundary
(alpha = -1/2 exactly or a zero coefficient/minor; never limit-stable).

## 2026-10-14  RUN STOPPED: 123 limit counterexamples at n = 7
Exact limit stage on the 1,889 strict survivors (fastregion limit mode with
replay): 1,726 never_positivity, 39 never (positivity possible, never limit-
stable), K7 stable, and 123 NON-COMPLETE graphs limit-stable, 0 escalations,
37.5 s. Every hit carries four exact pointwise witnesses (ctln.analyze_support
at beta inside the stable cell, t = 100, 1000, 10^4, 10^5), all stable
positive non-clique full supports at legal parameters. Critical values
include the PLASTIC NUMBER q^3 - q - 1 = 0 (q_c ~ 1.3247), higher than any
n = 6 reach (1.1664). Certificates: certificates/COUNTEREXAMPLE_CANDIDATES_n7_
limit.jsonl (+ _summary.json with D-equivalence classes, edge counts, twin
pairs). Per CLAUDE.md §3 the run is STOPPED pending human hand-verification;
no 2-D regions, structural analysis or theorems have been built on them.

## 2026-10-14  "Finite-eps map at n = 7 out of reach" was a cost estimate, not a fact; measured now
The earlier entry put the n = 7 two-parameter map out of reach on an estimate
(positivity survivors "in the millions", CAD "~1 s each"). Measured today:
* region2d at m = 7: 10.1 s/graph mean over 15 seven-node graphs (positive
  somewhere, never limit-stable), 9.9-13.2 s, all 'never', 0 escalations;
  the n = 6 map cost 1.9 s/graph (65,000 CPU-s / 34,071).
* New exact C prefilter n7/filter7p.c (domination-free AND some beta > 2 with
  all a_i(beta) > 0, the Lemma A.2 condition, decided by exact integer
  Cramer evaluation + a Descartes drop test on a_i(2+u); undecided lines go
  to n7/resolve_U.py = fastregion.region positivity_only, exact). Validated
  on ALL 1,540,944 six-node digraphs: keeps 28,739 (28,312 P + 427 U), the
  resolver keeps 28,350 (38 of the U lines positive, 389 never), and that set
  equals the n = 6 strict positivity survivor set EXACTLY (labelg canonical
  forms, diff zero in both directions); domination-free count 204,378 matches
  the report. 2.8 microseconds per graph including generation. The n = 7
  run over all 882,033,440 digraphs is in progress (four geng shards).
* n = 6 evidence on monotonicity in t: among the 4 stable n = 6 graphs no
  stable cell is bounded above in t (t_hi = None everywhere), consistent with
  "stable at a legal finite t => limit-stable at the same beta" (Lemma M,
  unproved; the general-matrix analogue is false, and at m = 2 the one-edge
  graph is stable exactly for illegal t < 1/(beta-2)). If Lemma M holds the
  n = 7 finite-eps map reduces to the 638,042 spectral survivors; if not the
  full positivity survivor set must be certified in two parameters.

## 2026-10-14  A cheap exact necessary condition kills 91% of the n = 6 positivity survivors
n7/cheap_prune.py decides exactly whether some legal (beta, t) has all a_i > 0
AND all coefficients of det(mu I + K) = phi + t*omega positive (Stodola; the
constraints are linear in t, so the feasible t-set at each beta is an
interval; root isolation of a_i, phi_k, omega_k, the pairwise numerators and
the legality numerators gives sign-invariant beta-cells, one exact sample
each). On the 28,350 strict n = 6 survivors: 25,812 cheap_never, 2,538
survive (9.0%), 38 ms/graph, in 1,085 s. All 4 stable n = 6 graphs, K5 and
all 124 limit-stable n = 7 graphs survive; 30 of the 39 n = 7 'never with
cert' graphs are killed. n = 5: 343 of 390 killed, K5 survives; 17,197 random
exact legal points on cheap_never graphs were all infeasible (self-check).
Consequence: the two-parameter CAD (10 s/graph at m = 7) is needed on roughly
a tenth of the positivity survivors, not all of them. Unexpected: the first
4,000 survivors in geng order (sparse) are ALL cheap_never; the survivors are
dense graphs, as the counterexamples are.

## 2026-10-14  n = 7 exact positivity prefilter complete: 9,632,133 of 882,033,440 digraphs
geng 7 r/4 | directg | filter7p, four shards, 49 minutes wall on 4 cores:
total 882,033,440 (= A000273(7), reproduced again), domination-free
117,232,922 (as before), dropped by the exact Descartes test 107,600,789,
kept 9,632,133 = 8,845,537 'P' (positivity certified at a rational beta > 2)
+ 786,596 'U' (undecided in C; exact Python resolution running; on shard 3,
3,417 of 35,006 U lines are positive, 9.8%). 2,442,088 kept graphs carry two
non-adjacent in-twins (Remark J.4', never stable). Kept fraction 1.09% vs
1.84% at n = 6.
Cheap necessary-condition pass (n7/cheap_prune.py) on a uniform random sample
of 20,000 shard-3 survivors: 2,479 survive (12.4%), 17,521 cheap_never,
64 ms/graph, 55,042 random exact legal points on cheap_never graphs all
infeasible. None of the 2,479 survivors carries non-adjacent in-twins (same
at n = 6: 0 of 2,538), so the J.4' prune is subsumed by the coefficient test.
Projected CAD workload for the complete finite-eps n = 7 map: about 1.1
million graphs at 10 s each, roughly 3,100 CPU-hours; the cheap pass itself
about 160 CPU-hours. Previously "out of reach"; it is a rented-compute job.

## 2026-10-14  RUN STOPPED again: a finite-eps-ONLY counterexample at n = 7, invisible in the limit
Workflow wf_945cf15e-c68 (one prover, two independent refuters, six exact referees) set out to
prove "Lemma M": stability at a legal finite t implies limit stability at the same beta
(which would have reduced the finite-eps n = 7 map to the spectral survivors). It is FALSE,
with exact certificates at every level:
* pure positive stability of K, m = 3: A = 011/001/100, beta = 22/5: stable at t = 11/25
  (legal, q = 187/180), unstable at t = 1 and t = 100 (also 011/101/100 at beta = 4 and
  011/100/010 at beta = 50/11). Positivity fails there, so no CTLN consequence.
* WITH positivity, n = 7: graph &FNFL~[udo? (labelg canonical &FS@{`]f|~?), 28 edges,
  domination-free, not a clique, spectral abscissa alpha(A|perp) = -1/2 EXACTLY
  (7 chi_perp = (l+1)(l^2+l+1)(7l^3+14l^2+14l+5); filter7 flag 'B'), hence NEVER limit-stable
  and NOT among the 123 limit candidates. ctln.analyze_support: eps = 1/49, delta = 33/1568
  (beta = 65/32, t = 48, q = 99/98): legal, x > 0, all Hurwitz minors positive -- a stable
  positive non-clique full-support fixed point; also eps = 1/36, delta = 17/588 (beta = 100/49,
  t = 35). At the same beta and t = 100, 1000: unstable (minor 6 negative, decided).
  region2d: stable_somewhere; at beta = 53/26 the stable t-range is (26, 53), bounded above.
  Found independently by both refuters (different labellings, recognised as isomorphic only by
  labelg), each confirmed by an independent exact referee. Certificate:
  certificates/COUNTEREXAMPLE_CANDIDATE_n7_finite_FNFL.json (re-runnable:
  n7/write_finite_certificates.py). Per CLAUDE.md section 3 the run is STOPPED: no sweep
  stage was launched after this; the candidate awaits human hand-verification together
  with the 123 limit candidates.
* "Bulge": for limit candidates &FNDL}[u|o? and &FK@dyUfXp? the finite-eps stable region
  extends to beta ABOVE beta_c and is bounded above in t there (e.g. &FNDL}[u|o?: beta =
  52/25 > beta_c in (56/27, 27/13), stable at t = 20, unstable at t = 60 and 10^5). Lemma AS
  (prover, PROVED): nu_A(c) = -alpha + R*/c + O(c^-2); the bulge occurs iff R* > 0; a float
  table predicts 39 of the 123 candidates bulge. So the t -> infinity classification
  understates the stable regions of some candidates (n = 6: none bulge; no stable cell
  bounded in t). certificates/n7_limit_hits_finite_eps_bulge_witnesses.json.
Consequence for the map: "finite-stable => limit-stable" fails, so the finite-eps n = 7 map
cannot be read off the limit stage, and the finite-eps counterexample set is strictly larger
than the limit set (at least 124 non-clique graphs at n = 7).

## 2026-10-14  Lemma R (ray reduction) makes the finite-eps n = 7 map cheap after all
Prover (PROVED, re-derived here): K(lambda beta, lambda t) = (1-lambda) I + lambda K(beta, t)
(Remark F), so stability propagates down every ray to the legality curve t(beta-2) = 1;
f(lambda) = lambda t (lambda beta - 2) - 1 changes sign on (2/beta, 1). Hence the full support
is stable at some legal (eps, delta) iff K(2 + 1/t, t) is stable for some t > 0 -- an exact
ONE-dimensional problem: t K = L(t) = tI + (2t+1)N + t^2 J, integer polynomial matrix,
critical set = positive roots of det L and Delta_{m-1}(t). n7/curve_test.py (from the prover's
curve1d.py): 35-45 ms per 7-node graph; prover's validation: all 9,608 five-node digraphs
(95 stable, each with an exact legal witness), all 124 limit-stable n = 7 graphs stable (as
Lemma R requires), 16,000 random n = 7 positivity survivors all 'never'. Strictly stronger than
the coefficient-sign prune (which it subsumes) at the same cost. Also PROVED and re-derived:
e_2(K) = B(1+2t) + Nn(1-beta)(1+beta+2t) + O(1 - t(beta-2)) over bidirected / non-adjacent /
one-way pairs, so legal stability forces B > Nn, i.e. edges >= m(m-1)/2 + 1 = 22 at n = 7
(4,045,243 of the 9,632,133 prefilter lines). Revised cost of the complete finite-eps n = 7
map: curve test on <= 8.9 M survivors (~100 CPU-hours, ~50 with the edge bound) + region2d on
the curve-stable graphs (hundreds, 10 s each). It fits on this 4-core machine in about a day;
previously estimated at 3,300 CPU-hours (CAD on 1.1 M graphs). NOT launched: the run is stopped.
Scratch-directory collision: two workflow agents wrote into the same scratch directory and one
overwrote the other's scan_pos.py; the prover moved its work and reported it. No repository file
was touched; future workflow prompts assign disjoint directories per agent label.

## 2026-10-14  Curve test validated: on the n = 6 positivity survivors it isolates EXACTLY the 4 stable graphs
n7/validate_curve_test.py: n = 5, all 9,608 digraphs: 9,513 never, 95 stable_somewhere,
0 escalate, 0 disagreements between curve_test and curve_test_fast, every stable graph
confirmed by an exact legal witness (ctln.analyze_support at beta = 2 + 1/t + 10^-6);
identical to the prover's independent count. n = 6, all 28,350 strict positivity
survivors: 28,346 never, 4 stable_somewhere (K6, Geneson idx 1000255, idx 1387528,
idx 1424629), 0 escalate, 192.6 s (6.8 ms/graph, the coefficient fast path decides most).
So at n = 6 "positive somewhere" and "stable somewhere in the legal region" together
already pin down exactly the stable-fixed-point graphs; the 2-D map was needed only for
the regions themselves. The prover's figure of 28 legal-stable n = 6 graphs counts ALL
digraphs (most not positive anywhere). Projected n = 7 cost of the curve stage: 8.9 M
survivors at ~10-15 ms = 25-40 CPU-hours. Not launched (run stopped).

## 2026-10-14  Audit workflow wf_56c77840-9ba: prefilter, Lemma R, curve test, edge bound, FNFL certificate -- no blocking issue
Six agents (notes/audit_n7_finite_workflow.json, artifacts in notes/audit_n7_finite/).
Prefilter + resolver: SOUND, both lenses. 979,800 exact sympy-vs-C a_eval values on 6,000
graphs, 0 mismatches; all 786,596 n = 7 U lines recomputed by an independent sympy-free
Fraction/Descartes decider, byte-identical (63,586 positive); n = 6 counts reproduced;
domination test matches Curto 2510.05098 Def 1.1 read from the paper; 105,000-graph
relabelling test, 0 decision changes. Lemma R + n7/curve_test.py + edge bound: SOUND.
Remark F re-derived symbolically; critical-set sufficiency confirmed by an independent
all-minors Routh-Hurwitz arrangement over ZZ[t] on 3,129 graphs (21,805 cells, 0
mismatches); 0 disagreements curve_test vs curve_test_fast on 3,000 random n = 7 graphs;
300 never-graphs x 20 random exact legal points, 0 stable; e_2 - (B-Nn)(1+2t) <= 0 at every
legal point, so edges >= 22 at m = 7. FNFL certificate: HOLDS under both lenses (own sympy
rebuild of W, exact RHP root counts 0/0/2; chi_perp factorisation two ways, alpha = -1/2
exactly; labelg non-membership among the 123 also by brute force over 5040 permutations;
region2d reproduced key-by-key). Non-blocking follow-ups applied: overflow guard on the
prefilter samples (proven safe iff max(2q+p', q) <= 193 for n <= 8; the current 14 samples
are), resolver keeps unexpected verdicts, assert -> raise in curve_test, adaptive witness
eta in the validator (&F]VhW[rSp? needs eta = 10^-7 at t = 1917: an exact witness at a fixed
offset is not guaranteed, only for small enough offset). Logged, not changed: CLAUDE.md
section 2's closed form I - W_sigma = -delta I + (1+delta) J - (eps+delta) A_sigma holds with
A^T (W reads A[j][i]); ctln.identity_check_IminusW already uses the transpose and every
computation is consistent, so this is a wording note for the human, per section 2. Also:
geng res/mod classes at n = 7 are very uneven (17 of 300 empty, some ~2e8 digraphs); the
total count check covers truncation. All prunes are proved for the FULL support only; the
survivor file must not be reused to prune size-5/6 supports (those reduce to induced
subgraphs by Theorem F anyway).

## 2026-10-15  Research workflow wf_a743d04d-73f: infinite families, q_n* lower bounds, n = 8 filter, scoop check
(notes/research_n7_parallel_workflow.json; artifacts notes/research_wf3/; theory notes
notes/research_wf3/theory/notes_theory.md. Every n >= 7 instance is a CANDIDATE pending
hand-verification; the theorems are exact algebra on the verified n = 6 bases.)
* THEOREM F (prover; REFEREE CONFIRMED): Geneson's graph plus a bidirected clique class of k
  nodes (closed-twin blow-up, n = 6 + k) is limit-stable for EVERY k >= 0: alpha < -1/2 for all
  k (Hurwitz minors 16(k+1), 64(k+1)(k+16), 8192(k+1)^2(3k+16)), positive on (2, 1 + plastic).
  So the target-free clique conjecture fails for every n >= 6, with sup q = 1.1664 (n = 6),
  1.2747 (n = 7), the plastic number 1.3247 (n = 8..14), then decreasing to 1. F_1 is one
  of the 123 n = 7 candidates.
* Lemma W / X (REFEREE CONFIRMED): closed-twin blow-ups preserve the full-support fixed-point
  data by explicit charpoly identities; 9 (n = 6) + 617 (n = 7) clique-class families proved
  for all k >= 1 by exact root isolation in k (not yet refereed).
* THEOREM E (not yet refereed): a finite-eps-ONLY family E_a (n = 7 + a), alpha = -1/2 EXACTLY
  for every a (Phi = z(z^2+z+1)((a+7) z(z^2+z+1) + 2)), stable positive at legal points with
  eps -> 0, beta -> 2+ (R* = 1/4). &FNFL~[udo? itself is the n = 6 hit idx 1387528 with one
  node doubled; the factor l^2 + l + 1 is a MOVING eigenvalue pair sitting on Re = -1/2.
* Growing q (not yet refereed; instances re-verified by a referee): families Q (n = 7 + k,
  q -> 1.8442), Q2 (n = 7 + 4k, q -> 3), Y (n = 8 + 5k, q_n* >= 3.8 for n >= 1823); exact
  lower bounds on the full-support q_n*: n = 8: 1709/1001, n = 9: 17691/10010, n = 10:
  18418/10001, n = 11: 17220/9091, n = 47: 25000/10001. Conjecture: sup_n q_n* = infinity
  (Geneson's safe bound n - 2 is the only upper bound known). Conjecture SL (every limit
  polynomial has a root on Re = -1/2) REFUTED exactly (55 of 135).
* Exact q* tool n7/qstar.py (not yet refereed): q_6* = root of 797q^6 - 1686q^5 + 2159q^4 -
  3004q^3 + 2499q^2 - 2158q + 1521 in (561/481, 7/6), Geneson's limit sextic, so at n = 6 the
  conjecture holds for all legal parameters with q >= q_6* and fails for every q in (1, q_6*).
  For &FNFL~[udo? the supremum of q over its stable region is a DEGREE-25 algebraic number in
  (685003/657781, 2170535/2084278) ~ 1.0414, attained along the Hurwitz curve Delta_6 = 0 at
  t ~ 25, NOT at the positivity boundary (beta_hi - 1 ~ 1.0830 is not reached). The maximum
  of beta_c - 1 over the 123 limit candidates is the plastic number (6 graphs); 6 of 123 have
  limit intervals starting strictly above q = 1. FNFL is 'undecided' (vanishing Hurwitz
  minor) at beta = 3 away from any boundary.
* n = 8 spectral filter n7/filter8.c (every integer operation overflow-guarded; overflow ->
  kept with flag 'O'): byte-identical to filter7 on all 1,540,944 n = 6 digraphs and on the
  638,042 n = 7 survivors; forced-overflow build loses nothing; 0 overflows on 155 M n = 8
  digraphs (max |charpoly coefficient| 5376). Rate ~0.9 M digraphs/s per core, directg-bound.
  Full n = 8 limit map: 500-800 CPU-hours, expected 1e5-5e5 strict survivors -- RENTED
  COMPUTE needed. Incidental n = 8 limit candidate &GLjlNZ}NzH\o (41 arcs) found in a random
  dense sample (notes/research_wf3/filter8/CANDIDATE_n8_limit_GLjlNZ.json).
* Scoop check: NOT scooped. Partial overlap: anonymous Zenodo/GitHub preprint "Six Is
  Minimal" (ipitchford/ctln-six-is-minimal, zenodo.23089682, 2026-10-01): n <= 5 exclusion,
  Geneson's graph as the 6-witness, a family on 5r + 6 vertices; no n = 6 classification, no
  n = 7, no finite-eps-only example, no threshold; its counts (9,846 classes; 1,052,741
  labelled) match neither A000273(5) = 9,608 nor 2^20. A claims diff is running.
* Hygiene: the n = 73 witness in family_Y_certificate.json uses a non-standard digraph6
  header (chr(n+63)); decode from the adjacency rows. Referee workflow wf_452fbfe3-f07 is
  checking Theorems E, Q, Q2, Y, the catalogue, claim 16 and the q* tool.

## 2026-10-15  n = 7 curve stage complete: 165 graphs are stable somewhere in the legal plane
n7/run_n7_curve.py over all 8,909,123 positivity survivors, 4 shards, 4.1 h wall (16.4 CPU-h,
6.6 ms/line): 5,135,775 dropped by the edge bound (edges <= 21), 3,773,183 'never' by the exact
Lemma-R curve test, 165 'stable_somewhere', 0 escalations, 0 errors. Every digraph on 7 nodes
is now accounted for by a proved necessary condition except these 165 (which include K7 and the
124 known candidates if the pipeline is consistent -- checked by the aggregator). Far below the
earlier projection of ~1.1 M two-parameter certifications. Two-parameter stage running on the 165.

## 2026-10-15  COMPLETE finite-eps map at n = 7: exactly 124 non-clique stable full supports
Two-parameter stage (region2d, crosscheck + replay) on the 165 curve-stable graphs: 125
stable_somewhere (K7 + 124 non-complete), 40 never (stable somewhere but never positive at
the same point), 0 escalations, 0 errors, 3,467 CPU-s. Aggregation (n7/aggregate_n7.py,
accounting_ok = True): 882,033,440 digraphs = 764,800,518 dominated (Lemma D) + 107,600,789
never positive (exact Descartes, Lemma A.2) + 723,010 never positive (exact resolver) +
5,135,775 with <= 21 edges (e_2 bound) + 3,773,183 never stable on the legality curve (Lemma R)
+ 40 stable-but-never-positive-together (CAD) + 125 stable positive (CAD). The 124 non-clique
graphs are EXACTLY the 123 limit candidates plus &FNFL~[udo?: no other finite-eps-only graph
exists at n = 7. 40 of the 124 have a stable cell bounded above in t (the "bulge": the limit
classification understates their regions); 91 D-equivalence classes (71 singletons, 16 pairs,
two of size 4, one of 5, one of 8); 64 carry bidirected in-twins; edge counts 23..31.
Certificates: certificates/n7_finite_summary.json, certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl
(every hit with its full two-parameter certificate, exact limit status and class key),
certificates/n7_region2d/ (all 165 records). All 124 remain CANDIDATES pending hand-verification.

## 2026-10-15  Referee workflow wf_452fbfe3-f07: every refereed theory claim holds; Conjecture P refuted; preprint count note corrected
(notes/referee_n7_theory_workflow.json; artifacts notes/referee_wf4/.) Four independent referees,
all "holds, high confidence": Theorem E + Lemma B (finite-eps-only family, R* = 1/4 for all a,
every proof step re-derived, instances a = 0..5 recomputed; write-up gaps: pick a rational gamma
in Lemma B, mention the roots at z = 1 in step 2, a missing 3c/4 coefficient in the notes, an
unrecorded AMAX = 5 run); Theorems Q, Q2, Y (thresholds k >= 25/79/286, 3/10/33, 13/24/46/117/363
proved sharp by exact root isolation in k; Y covers n = 8 + 5k only; the "many-class ceiling" is
a class-size-limit statement, not a bound for finite classes); the 626-entry catalogue (all
re-proved independently, for all REAL k >= 1; entry 616, base &FNFL~[udo?, fails at k = 0 exactly
as alpha = -1/2 predicts) and claim 16 (21 of 123 candidates contain an n = 6 hit; none contains
K6); the q* tool (q_6* confirmed as Geneson's sextic to 4e-31; the FNFL supremum, a degree-25
number ~1.0414, is the CORNER where Delta_6 = 0 meets the positivity line beta = beta_pos at
t ~ 25.009, not the t -> infinity end; qstar_fnfl.json's aggregate block is stale and reports only
the weaker sample bound 1.0195; --selftest is documented but not implemented).
Shared dependency of every "for all k" statement: Lemma I.3/I.4 (limit-to-finite-t transfer,
notes/lemma_I_spectral.md), refereed in the n = 6 phase and corroborated at every sampled finite t
here, not re-proved by these referees.
CONJECTURE P REFUTED exactly (side finding of the catalogue referee): 33 of the 54 "S = all"
catalogue entries, e.g. base &FMAlW\d\y? with chi_D = z^7 - z^5 - 7z^4 - 10z^3 - 9z^2 - 4z - 1,
have every eigenvalue of D = J - I - A with Re > -1/2 AND u = (D^T + rI)^-1 1 > 0 at r = 49/100;
so positivity does not force an eigenvalue at or beyond Re = -1/2, and the joins G v K_k are
limit-stable for every k >= 1 (notes/research_wf3/theory/notes_theory.md item 13 is wrong as
stated OPEN; the catalogue proves the opposite).
Correction (append-only) to the entry of 2026-10-15 on the "Six Is Minimal" preprint: its counts
are NOT inconsistent -- 9,846 = A000273(1) + ... + A000273(5) and 1,052,741 = sum of 2^(m(m-1)) for
m = 1..5; my note compared them with the m = 5 values alone. The claims diff (literature agent;
their repository cloned read-only): their conventions equal ours up to a transpose of the nonedge
matrix; their n <= 5 exclusion equals our n = 5 map; their Geneson-class region gives
1 + u* = 1.16635502785 = our q_6* (same sextic, same degree-14 beta polynomial): an INDEPENDENT
confirmation of q_6*; their 5r + 6 family is a Lemma-W clique-expansion of the same seed with
q -> 1. They have no n = 6 classification, no n >= 7 statement, no finite-eps-only example, no
Theorem F, no growing-q families. They have an explicit rational parameter schedule for all r and
global NONDEGENERACY certificates (parity trick D K = I mod 2); Geneson's conjecture is stated for
nondegenerate CTLNs, so n7/nondegeneracy_check.py now certifies nondegeneracy exactly (all 2^n - 1
supports: nonsingular, no zero coordinate, no zero drive) at a witness point of every hit.
Their open question, which Theorem F does not answer: at FIXED (eps, delta) the clique-class family
F_k is stable only up to some k (exact: at eps = 1/21, delta = 11/210 stable for k <= 38, unstable
at k = 41, 45, 60), so unbounded support size at a fixed parameter remains open for both groups.
Positioning: cite as an unrefereed preprint; "independently obtained / consistent with"; no
priority language.

## 2026-10-15  Nondegeneracy certified exactly at every witness point
n7/nondegeneracy_check.py: for Geneson's graph at the paper's parameters, the two other n = 6
hits at their witnesses, and all 124 n = 7 candidates at a rational point inside their first
stable cell: every one of the 2^n - 1 supports has det(I - W_sigma) != 0, a fixed point with no
zero coordinate and no off-support drive exactly 0, and the full support is stable there. So all
witnesses are nondegenerate CTLNs in the sense the conjecture requires.
certificates/nondegeneracy_witnesses.json.

## 2026-10-15  Lean 4 + Mathlib available after all; 26 theorems machine-checked (no sorry)
elan installs; the Mathlib olean cache downloads (9,031 archives in 95 s); toolchain
leanprover/lean4:v4.35.0, Mathlib tag v4.35.0. notes/lean/CtlnCert/ (sources, lakefile, manifest,
build and axiom-audit logs; the 7.3 GB .lake directory is not committed): Remark F's affine
identity K(l beta, l t) = (1-l) I + l K(beta, t) over any commutative ring; the characteristic
polynomial scaling eval x charpoly(c I + l M) = l^m eval((x-c)/l) charpoly(M); eigenvalues of
K along a ray are (1-l) + l nu; positive stability propagates down rays (0 < l <= 1); the IVT
step and the explicit root lamStar = (t + sqrt(t^2 + beta t))/(beta t) of l t (l beta - 2) = 1 in
(2/beta, 1), unique positive root; the CLAUDE.md section 2 closed form I - W = -delta I +
(1+delta) J - (eps+delta) A^T (machine-checking the transpose reading); I - W = eps K(J - I - A^T,
beta, t); q > 1 <-> t(beta-2) > 1; Hurwitz stability of -I + W <-> positive stability of K; and
Lemma R's forward direction in eps/delta language (stable at a legal point => stable at a point
with q = 1 exactly). Every theorem depends only on propext, Classical.choice, Quot.sound; 0
sorryAx. NOT formalised: the converse of Lemma R (openness of stability needs eigenvalue
continuity, absent from Mathlib) and everything spectral/exact-arithmetic downstream.

## 2026-10-15  Lean batch 2: the pipeline's prune lemmas are machine-checked (57 theorems total, no sorry)
notes/lean/CtlnCert/{LemmaA2,EdgeBound,Domination}.lean. Lemma A.2: positive stability => det K > 0
(proved by a continuity argument: x -> det(xI + K) is monic with no root in [0, inf), so positive at
0), Cramer in adjugate form, and the CTLN statement "under Hurwitz stability, x > 0 iff adj(K) 1 > 0"
(positivity is t-free). The e_2 edge bound: e_2(K) = B(1+2t) + Nn(1-beta)(1+beta+2t) + O(1-t(beta-2))
over any commutative ring for 0/1 zero-diagonal A, pair counting, and -- the spectral step included,
via Mathlib's charpoly coefficients as sums of principal minors and Gauss-Lucas -- positive stability
at a legal point implies Nn < B and m(m-1) + 2 <= 2 #arcs, in eps/delta language. Lemma D: a
dominated node forces x_j = 0 at any fixed point with theta_j <= theta_k (so no positive full-support
fixed point), including the induced-subgraph form. Axiom audit: 57 theorems, only propext,
Classical.choice, Quot.sound; 0 sorryAx. Unexpected: nothing had to be left as a hypothesis.

## 2026-10-15  Lean batch 3: Stodola and Remark J.4' machine-checked (87 theorems total, no sorry)
notes/lean/CtlnCert/{Stodola,Twins}.lean. Stodola: a monic real polynomial whose complex roots
all have negative real part has every coefficient positive (strong induction on the degree,
splitting off a real linear or a conjugate-pair quadratic factor; products of positive-coefficient
polynomials); corollaries: every coefficient of charpoly(-K) is positive when K is positive stable
(the k = 0 case is det K > 0, agreeing with Lemma A.2), and at the CTLN level for -I + W. This is
the statement behind the curve test's fast path (coeff_never); the two remaining links -- that
L(t) = tK(2 + 1/t, t) has K's stability for t > 0 (positive scaling, eig_affine) and that the
Python builds exactly this L(t) -- are stated in the README, not as Lean theorems. Twins: for
in-twins u, v (same in-neighbours off {u, v}) the vector e_u - e_v is a left eigenvector of W with
eigenvalue 1 + delta (non-adjacent) or 1 - eps (bidirected); hence -I + W has the eigenvalue delta
> 0 and is never Hurwitz stable when two non-adjacent in-twins lie in the support, for every eps
and every delta > 0, including the induced-subgraph form. Axiom audit: 87 theorems, standard
axioms only, 0 sorryAx.

## 2026-10-15  The bulge criterion, exactly: bounded-in-t stable cells <=> R* > 0 AND positive at beta_c (124/124)
Exact R* (Lemma AS) for all 124 n = 7 hits (notes/bulge_exact/, certificates/n7_bulge_Rstar_exact.json):
51 positive, 73 negative, none zero; every critical root simple and not shared with chi. R* > 0 alone
matches the 40 graphs with a stable cell bounded above in t on 113 of 124; the 11 exceptions all have
R* > 0 but FAIL positivity at beta_c (their limit interval ends at a positivity root a_k below the
Hurwitz boundary), so their bulge beyond beta_c is pure stability, never a stable fixed point. With
the positivity condition the rule is exact: 40 true/true, 84 false/false. All 40 bounded cells start
exactly at beta_c (gcd test against the CAD critical points); no bounded cell occurs anywhere else.
FNFL: R* = 1/4 exactly, beta_c = 2. The SIX PLASTIC-NUMBER graphs (largest limit q_c) do NOT bulge:
their limit interval ends at 1 + rho because positivity ends there, with the stability boundary
further right (beta_c ~ 2.33 to 2.5). Validation: the right end of all 93 Delta-terminated limit
intervals equals the exact beta_c (0 mismatches); the float first-order table agreed on sign(R*) for
all 123; the CAD t_hi values match beta_c^3 R*/(beta - beta_c) to first order (e.g. 13.4 vs 16).

## 2026-10-15  q_7* is the plastic number, exactly, and it is not attained
From the 124 two-parameter certificates: every stable cell unbounded in t has sup q = beta_hi - 1 with
beta_hi a limit endpoint (max 1 + rho, six graphs); every one of the 40 bounded cells has the outer
isolating bound beta_hi <= 16/7, so q < 9/7 there; hence q_7* = rho = root of q^3 - q - 1 ~ 1.3247,
approached only as t -> infinity. Four of the six plastic-number graphs have limit intervals starting
at q = 1, so the union of the candidates' intervals covers (1, rho) with no gap. Conditional on hand-
verification: for n <= 7 the conjecture holds for every legal parameter with q >= rho and fails for
every q in (1, rho). The exact qstar --sharpen run (independent machinery) is still computing; its
output is the confirmation certificate (certificates/qstar/qstar_n7_finite.json when done).

## 2026-10-15  Human hand-verification of the 124 n = 7 candidates recorded; independent machine verification launched
The user reports having hand-verified all 124 candidates (certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl)
without full confidence and asked for independent machine verification. Workflow wf_e7b72ed9-c86:
four verifiers with separately written exact implementations (sympy Rationals + own Hurwitz minors;
pure Fractions with own elimination/charpoly/minors; exact right-half-plane root counting without
Hurwitz minors; structural and limit claims incl. pairwise non-isomorphism, domination-freeness and the
exact limit test), none using the project's code. Sign-off is recorded only after this workflow
reports every graph confirmed by at least three implementations. Rented compute for n = 8 is deferred
by the user.

## 2026-10-15  ALL 124 n = 7 COUNTEREXAMPLES VERIFIED: human + four independent exact implementations
Workflow wf_e7b72ed9-c86 (notes/verify_124_workflow.json, artifacts notes/verify_124/): four verifiers
written from the CLAUDE.md definitions without the project's code. sympy/Hurwitz: 124/124 at 249
prescribed points and all 492 stored witnesses (flags and x identical to the files). Pure Fractions:
124/124 at 422 points (first legal stored witness + an interior point of EVERY stable cell, 298 cells),
stored unstable points re-decided unstable with strictly negative minors; with the transpose
convention reversed 0/124 survive. Root counting without Hurwitz minors (exact Routh, Sturm/argument
principle and exact isolation, all three required to agree, calibrated on 3,000 random polynomials):
124/124 at 742 points, 0 undecided. Structural: 124 pairwise non-isomorphic (own canonical form and
nauty), 124 domination-free, none complete (11-19 missing arcs), the 123 limit candidates identified
up to isomorphism, limit stability of all 123 confirmed exactly at the cell sample (given Lemma I),
FNFL's spectral abscissa -1/2 exact. Certificate statuses updated; the user's hand-verification of all
124 (2026-10-15) is recorded alongside. Data-hygiene finding: 15 of the 123 limit candidates have an
ILLEGAL first stored pointwise witness (t = 100, q <= 1, already flagged legal = false; &F]VhW[rSp?
also at t = 1000) -- a hand check that used the first witness would have examined a non-CTLN point;
certificates/n7_candidates_legal_witnesses.json now lists one legal witness per graph plus an
interior cell point. The field agrees_with_conjecture is meaningless at those illegal points. A
verifier's interim run flagged 46 graphs because repeated eigenvalues (twin structure, off the
imaginary axis) were treated as undecided by one method; fixed, final: 0 undecided.
