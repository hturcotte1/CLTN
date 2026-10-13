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
