# Issue 2 of the external report (2026-11-01): thm:AS invoked at nu_inf = 1/2 in lem:Bmech

Status: the report's point is correct as a citation gap; the mathematics of lem:Bmech and thm:famE is unaffected.
Exact check: `check_AS_endpoint.py` in this directory (`python3 check_AS_endpoint.py`, exit 0 iff all 30 checks
pass; output `check_AS_endpoint_out.json`, transcript `check_AS_endpoint_stdout.txt`). All arithmetic is exact
(sympy over QQ, the field QQ[z]/(z^2 - z + 1), fractions.Fraction, exact Hurwitz minors via `ctln.py`, exact
real-root counting). Conventions: the project conventions file, section 2; paper Section 5 ((E): K = beta(cJ - A^T - gamma I),
c = 1 + t/beta, gamma = 1 - 1/beta, beta_c = 1/(1 - nu_inf)).

## 1. Where the range nu_inf in (1/2, 1) is used

The paper's thm:AS bundles two statements of notes/lemma_M_refuted.md section 4:

* **Lemma AS** (no hypothesis on Re z0): at a simple root z0 of psi1 with chi1(z0) != 0, F(z, e) = psi1(z) - e chi1(z)
  has F(z0, 0) = 0 and F_z(z0, 0) = psi1'(z0) != 0, so the implicit function theorem gives an analytic branch
  z(e) = z0 + e chi1(z0)/psi1'(z0) + O(e^2), and chi1(z0)/psi1'(z0) = 1/s^'(z0) = chi(z0)/psi'(z0) because psi1(z0) = 0
  (also psi'(z0) = h(z0) psi1'(z0), chi(z0) = h(z0) chi1(z0)).
* **Theorem AS**, proof step 1: for large c the eigenvalues of cJ - A are the fixed roots of h = gcd(chi, psi), the
  moving roots of chi1 - c psi1 (all but one converge to the roots of psi1; one real one tends to +inf, since deg chi1 =
  deg psi1 + 1). Under the hypothesis "every root of psi on Re z = nu_inf is simple and not a root of chi", no root of h
  lies on that line and every other root of psi has real part > nu_inf, so nu_A(c) = min Re z(1/c) over the critical
  branches = nu_inf + R*/c + O(c^-2) (a minimum of finitely many analytic functions of e = 1/c).

Neither step uses nu_inf > 1/2 (nor nu_inf < 1, except to make beta_c finite). The boundary formula
beta = beta_c + beta_c^3 R*/t + O(t^-2) is obtained by solving nu_A(c) = gamma with gamma = 1 - 1/beta, c = 1 + t/beta,
and also uses only the expansion; at nu_inf = 1/2 it reads beta = 2 + 8R*/t + O(t^-2).

The range is used in exactly one place: proof step 2 of (a) (lemma_M_refuted.md section 4) needs a single c2 that is
legal for every gamma in [nu_inf, (1 + nu_inf)/2], which is possible because gamma^2/(2 gamma - 1) is bounded there, i.e.
because nu_inf > 1/2. In words: for nu_inf > 1/2 the legality constraint is invisible near (beta_c, inf); at
nu_inf = 1/2 it is not.

## 2. What is true at nu_inf = 1/2

Legal and stable at (c, gamma) iff nu_P2(c) < gamma < nu_A(c), with nu_P2(c) = c - sqrt(c^2 - c) (lemma_M_refuted.md (3)).
Exactly, nu_P2(c) = 1/(1 + sqrt(1 - 1/c)) = 1/2 + 1/(8c) + 1/(16c^2) + 5/(128c^3) + ..., and (checked by exact root
counting after clearing the square root)

    1/2 + 1/(8c)  <  c - sqrt(c^2 - c)  <  1/2 + 1/(8c) + 1/(8c^2)        (lower for all c > 1, upper for all c >= 2).

(Lower bound by hand: with s = sqrt(1-x), x = 1/c, it is s < (4-x)/(4+x), i.e. (1-x)(4+x)^2 < (4-x)^2, i.e.
-8x^2 - x^3 < 0.) So at nu_inf = 1/2 both the stability boundary and the legality boundary are 1/2 + O(1/c), and the
comparison is between R* and 1/8:

* **(a) as stated is false at nu_inf = 1/2 in general**: if 0 < R* < 1/8 then for large c,
  nu_A(c) <= 1/2 + R*/c + K/c^2 < 1/2 + 1/(8c) < nu_P2(c), so no legal point with c large (beta bounded, t large) is
  stable, although R* > 0. Hence thm:AS cannot simply be restated with nu_inf in [1/2, 1) for all of its parts.
* **(b)** at nu_inf = 1/2 is vacuous (R* < 0 < 1/8 gives an empty stable set near (2, inf)).
* **R* > 1/8**: nu_A(c) - nu_P2(c) >= (R* - 1/8)/(2c) > 0 for c >= c1, so a rational gamma in between exists for every
  large c; beta = 1/(1 - gamma) -> 2+. This is lem:Bmech.

In the (beta, t) plane: legality is beta > 2 + 1/t exactly, stability (to first order) beta < 2 + 8R*/t. For E_0,
R* = 1/4 gives beta < 2 + 2/t: at beta = 53/26 the first-order prediction of the upper t-end is 8R/(beta - 2) = 52
against legality t > 26, consistent with the recorded window 26 < t < t_1, t_1 in (53, 54) (paper Section 5).

## 3. Decision (minimal correct fix)

The expansion holds at nu_inf = 1/2, so (per the brief) thm:AS is restated with nu_inf in [1/2, 1) for the expansion and
the boundary curve; (a) and (b) are kept for nu_inf > 1/2 and a new item (c) records the endpoint (threshold 1/8
instead of 0, citing lem:Bmech for the positive half). The proof sketch of thm:AS gains one sentence saying where the
range is used and one for (c). The proof of lem:Bmech is changed to cite exactly what it uses: the expansion of
thm:AS at nu_inf = 1/2 (not (a), (b)), with the hypothesis translation (i)+(ii) => nu_inf = 1/2 and the critical roots
are exactly z0, z0bar, simple, not roots of chi; and R* = R because 1/s^'(z0) = chi(z0)/psi'(z0) and conjugate roots
have the same real part. No circularity: (c)'s positive half is a pointer to lem:Bmech, whose proof uses only the
expansion.

Source notes: notes/research_wf3/theory/notes_theory.md section 4.2 (Lemma B) already cites the local **Lemma AS**
(which has no range hypothesis), not Theorem AS, so it has no gap; a one-clause clarification is proposed. In
notes/lemma_M_refuted.md section 4 the restriction nu_inf in (1/2, 1) on Theorem AS is correct as stated there (it is
needed for (a)); a remark recording the endpoint case is proposed so that the paper's restated thm:AS has a source.

## 4. Exact check on E_0 = &FNFL~[udo? (script output)

* digraph6 decodes to the paper's rows 0011110 0011100 1101111 1110111 0011011 0100101 1100000; chi and psi are
  invariant under A -> A^T; psi(z) = 7 chi_perp(-z) with the paper's factorisation, i.e.
  psi(z) = (z - 1)(z^2 - z + 1)(7z^3 - 14z^2 + 14z - 5) (exact identity checked).
* z^2 - z + 1 (irreducible over QQ; roots z0, z0bar = 1/2 -+ i sqrt(3)/2, real part exactly 1/2) divides psi exactly
  once: simple. gcd(chi, z^2 - z + 1) = 1: chi(z0) != 0.
* Cofactor C(z) of psi: C(1/2 - x) has coefficients 7, 7, 7, 13/4, 5/16 and Hurwitz minors 7, 105/4, 70, 175/8, all
  strictly positive: every other root of psi has Re z > 1/2. Hence nu_inf = -alpha = 1/2 exactly.
* In QQ[z]/(z^2 - z + 1): chi(z0)/psi'(z0) = 1/4 - (sqrt(3)/12) i; cross-checked with radicals and as 1/s^'(z0).
  R = 1/4 > 1/8 (Fraction comparison).
* Legality boundary: series of 1/(1 + sqrt(1 - x)) to O(x^5) as above; identity c - sqrt(c^2 - c) = 1/(1+sqrt(1-1/c));
  the two-sided bound above by exact root counting (lower: numerator x^2(x + 8) times denominator (x + 4)^2, i.e. x^3 + 16x^2 + 80x + 128 after removing x^2,
  no root in [0, 1]; upper polynomial of degree 7, no root in [0, 1/2]).
* Finite-c bracket via the repository verifier `ctln.analyze_support` (exact, three-valued): at c = 10, 100, 1000, 10^4
  and gamma = 1/2 + 3/(16c) the point is legal and the full support is positive with all Hurwitz minors positive (an SFP;
  e.g. c = 10: eps = 77/1517, delta = 83/1517); at gamma = 1/2 + 5/(16c) the point is legal and the full support is
  decided unstable (Hurwitz minor 6 negative; at c = 10 this point is also outside the positivity cell). So
  3/16 < c(nu_A(c) - 1/2) < 5/16 at each tested c, consistent with R = 1/4. No minor vanished.
