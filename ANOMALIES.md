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
