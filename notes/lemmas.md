# Lemmas: rank-one structure, boundary sets, the eps -> 0 limit, domination

Prover notes, 2026-10-13. The conventions are exactly CLAUDE.md §2; nothing here changes them.

Each statement below is one of three things:

- **PROVED here.** The proof is in this file.
- **CITED.** A classical theorem, with a precise reference.
- **UNPROVED.** Labelled as such.

sympy was used only for sanity checks (listed in the Appendix). No identity below rests on a
test. The one computational claim (the Geneson sextic, §4) is exact, and its re-runnable code
is inlined in the Appendix.

| item | status | hypotheses actually used |
|---|---|---|
| A (i)-(v), A.1, A.2 | PROVED | A: eps > 0 only; A.1, A.2: (H) |
| B (a)-(e), B-inf | PROVED, citing Routh-Hurwitz and Rouché; the Orlando step is proved self-contained (B.1 b3) | (H) |
| C (a)-(d) | PROVED (Rouché) | beta > 1; legality only for "legal" |
| C' (compression), Lemma I (proof of notes/limit_spectral.md Claims 1-3) | PROVED; Claim 3 holds with NO nondegeneracy hypothesis | beta > 1 |
| D | PROVED. **eps < 1 is NOT needed**; theta_j <= theta_k IS needed (counterexample given) | eps > 0, delta >= 0 |
| 0 (cliques), E, E.1-E.4 | PROVED | (H) |
| F (ray monotonicity of stability) | PROVED | none |
| Review of fastregion.py (§8) | no soundness bug found; 1 reporting gap; 4 hardening/efficiency items | - |

---

## 0. Setting

Graph: `A[i][j] = 1` iff `i -> j`, nodes `0..n-1`.
`W[i][j] = 0` if `i = j`; `-1+eps` if `A[j][i] = 1`; `-1-delta` otherwise.

**(H)** is the standing hypothesis `0 < eps < 1, delta > 0, theta > 0`, with the same `theta` at every node.
Where a lemma needs less, it says so.

Reparametrisation: `t := (1-eps)/eps`, `beta := 1 + delta/eps`. The map `(eps,delta) -> (t,beta)` is a bijection
`(0,1) x (0,inf) -> (0,inf) x (1,inf)` with inverse `eps = 1/(1+t)`, `delta = (beta-1)/(1+t)`. Also
`q := delta(1-eps)/eps = (beta-1) t/(t+1)`, and

    legal  <=>  q > 1  <=>  (beta-1) t > t+1  <=>  t(beta-2) > 1  <=>  beta > 2 + 1/t.

Support `sigma = (sigma_1 < ... < sigma_m)`, `m >= 1`. Indices `a, b` run over `1..m`.

- `N` (m x m), receiver-first nonedge matrix: `N_ab = 1` iff `a != b` and `A[sigma_b][sigma_a] = 0` (sigma_b does not send to sigma_a).
- Off-support row, for `k` not in `sigma`: `N_kb = 1` iff `A[sigma_b][k] = 0`.
- `C := I + beta N`, `K := C + t J`, `J = 1 1^T`.
- `e := 1^T N 1` = number of ordered non-adjacent pairs inside sigma. `A_sigma` = induced adjacency, `(A_sigma)_ab = A[sigma_a][sigma_b]`, so `N = J - I - A_sigma^T`.
- Input to node `l` at a state `x in R^n`: `y_l(x) := theta + sum_i W[l][i] x_i`. For `x in R^m` on sigma, `x^ in R^n` is its extension by zeros.

**SFP(sigma)** ("sigma supports a stable fixed point") means there is `x in R^m` with:

1. `(I - W_sigma) x = theta 1` and `x > 0` entrywise;
2. for every `k` not in sigma, `y_k(x^) < 0`;
3. every eigenvalue of `-I + W_sigma` has negative real part.

If (3) holds then `I - W_sigma` is nonsingular, so `x` is unique. If (3) fails, SFP is false whatever (1) and (2) say.
So "the solution of (1)" is never ambiguous where it matters.

**Hurwitz.** A real polynomial is *Hurwitz* if all its roots have negative real part. For
`f(z) = c_0 z^n + c_1 z^(n-1) + ... + c_n`, the k-th Hurwitz minor is
`Delta_k(f) := det[c_(2i-j)]_(i,j=1..k)`, with `c_l := 0` for `l < 0` or `l > n`, and `Delta_0 := 1`. This is exactly
the matrix built by `ctln.hurwitz_stable` and `fastregion.hurwitz_det_poly`. Some books use the transpose; leading principal
minors are unchanged.

---

## 1. Lemma A (rank-one structure)

Two tools, valid over any commutative ring `R`. We use `R = Z[beta, t, mu]`.

**(T1) Cramer, column form.** For `M in R^(m x m)` and `b in R^m`: `(adj(M) b)_i = det M[i<-b]`, where `M[i<-b]` is `M`
with column `i` replaced by `b`.
*Proof.* `adj(M)_ij = (-1)^(i+j) det M(j|i)`, where `M(j|i)` deletes row j and column i. So
`(adj(M) b)_i = sum_j (-1)^(i+j) b_j det M(j|i)`, which is the Laplace expansion of `det M[i<-b]` along column i. ∎

**(T2) Rank-one update.** For `M in R^(m x m)` and `tau in R`: `det(M + tau 1 1^T) = det M + tau 1^T adj(M) 1`.
*Proof.* Column j of `M + tau 1 1^T` is `M_*j + tau 1`. Expanding det multilinearly in the columns gives a sum over
subsets `S` of `[m]` of `det(column j = tau 1 for j in S, = M_*j otherwise)`. Terms with `|S| >= 2` have two equal
columns and vanish. So `det = det M + tau sum_i det M[i<-1] = det M + tau sum_i (adj(M) 1)_i` by (T1). ∎
No invertibility is assumed. This adjugate form of the matrix determinant lemma is why `det C = 0` needs no special
treatment anywhere below.

**Lemma A.** Assume `eps > 0`. Then:

- **(i)** `(I - W)_sigma = eps K`.
- **(ii)** `adj(K) 1 = adj(C) 1 =: a = (a_1..a_m)`, for every value of `t` (not only `t > 0`). Each `a_i` is in `Z[beta]` with degree `<= m-1`.
- **(iii)** `det K = d + t s`, where `d := det C` (in `Z[beta]`, degree `<= m`) and `s := 1^T a`.
- **(iv)** `det(mu I + K) = phi(mu) + t omega(mu)`, where `phi(mu) := det(mu I + C)` and `omega(mu) := 1^T adj(mu I + C) 1`.
  `phi` is monic of degree m in mu, with `phi(0) = d`. `omega` has mu-degree **exactly** `m-1`, leading coefficient `m`, and `omega(0) = s`.
  All coefficients are in `Z[beta]`.
- **(v)** If `D := d + t s = det K != 0`, then the unique solution of (1) is `x = (theta/eps) a / D`. For every `k` not in sigma,

      y_k(x^) = theta * o_k / D,     o_k := d - beta * sum_b N_kb a_b.

*Proof.*

- (i) Using `eps*beta = eps + delta` and `eps*t = 1 - eps`, entry (a,b) of `eps K` is `eps[a=b] + (eps+delta) N_ab + (1-eps)`.
  - `a = b`: `N_aa = 0`, so the entry is `1`, which equals `1 - W[sigma_a][sigma_a]`.
  - `a != b`, `N_ab = 0` (sigma_b -> sigma_a): the entry is `1 - eps = 1 - W[sigma_a][sigma_b]`.
  - `a != b`, `N_ab = 1`: the entry is `1 + delta = 1 - W[sigma_a][sigma_b]`. ∎
- (ii) In `K[i<-1]`, column j (j != i) is `C_*j + t 1` and column i is `1`. Subtracting `t` times column i from every
  other column preserves the determinant and gives `C[i<-1]`. So `det K[i<-1] = det C[i<-1]` for every `t`, and (T1)
  gives the claim.
  Degree: in `C[i<-1]` the entries of column i are constants, and all other entries lie in `{0, 1, beta}`. Each Leibniz
  term takes exactly one entry from column i, so it has degree `<= m-1`. ∎
- (iii) Apply (T2) with `M = C`, `tau = t`. ∎
- (iv) Apply (T2) with `M = mu I + C` over `Z[beta, mu]`.
  - `phi`: only the identity permutation reaches `mu^m`.
  - `omega = sum_i det (mu I + C)[i<-1]`. In `(mu I + C)[i<-1]`, `mu` occurs only in the `m-1` diagonal positions `(l,l)`
    with `l != i`, each as `mu + 1` (because `C_ll = 1`). A Leibniz term for a permutation `pi` has mu-degree
    `<= #{l != i : pi(l) = l}`. The identity gives `(mu+1)^(m-1) * 1`, since the (i,i) entry of the new column is 1.
    Any `pi != id` fixes at most `m-2` points. So every summand has leading term `mu^(m-1)`, and `omega` has leading term
    `m mu^(m-1)`.
  - `omega(0) = 1^T adj(C) 1 = s` and `phi(0) = det C = d`. ∎
- (v) By (i), (1) reads `K x = (theta/eps) 1`. If `det K != 0`, then `x = (theta/eps) adj(K) 1 / det K = (theta/eps) a/D` by (ii) and (iii).
  For `k` not in sigma, `W[k][sigma_b]` is `-1+eps = -eps t` if `sigma_b -> k` (`N_kb = 0`), and `-1-delta = -eps(t+beta)`
  if not (`N_kb = 1`). That is, `W[k][sigma_b] = -eps (t + beta N_kb)`. Hence

      y_k = theta - eps sum_b (t + beta N_kb) x_b = theta - theta (t s + beta sum_b N_kb a_b)/D = theta (D - t s - beta sum_b N_kb a_b)/D = theta o_k / D. ∎

  Nothing divides by `det C`. So the case `det C = 0` needs no continuity argument: (ii), (iii) and (v) hold at every
  `(beta, t)` with `det K != 0`. Geneson's Lemma 2.1 form, `K^-1 1 = C^-1 1 / (1 + t 1^T C^-1 1)`, needs `det C != 0`;
  this one does not.

**Corollary A.1 (pointwise criterion).** Under (H), at `(beta, t)`:

    SFP(sigma)  <=>  D != 0,  a_i/D > 0 for all i,  o_k/D < 0 for all k outside sigma,  and  p(mu) := det(mu I + K) is Hurwitz.

*Proof.* `-I + W_sigma = -eps K` with `eps > 0`. So (3) holds iff every eigenvalue `lambda` of `K` has `Re lambda > 0`, iff
every root `-lambda` of `p` has negative real part, iff `p` is Hurwitz. If `p` is Hurwitz then `D = p(0) = det K != 0` and
(v) applies: (1) holds iff `(theta/eps) a_i/D > 0` for all i, iff `a_i/D > 0`. (2) holds iff `theta o_k/D < 0`, iff `o_k/D < 0`.
If `p` is not Hurwitz, both sides are false. ∎

**Lemma A.2 (sign of D; permission is t-free).** If `p` is Hurwitz then `D > 0`. Hence, under (H),

    SFP(sigma)  <=>  a_i(beta) > 0 for all i,  o_k(beta) < 0 for all k outside sigma,  and  p Hurwitz.

*Proof.* `D = det K` is the product of the eigenvalues of `K`. The real ones are `> 0`, and the non-real ones come in
conjugate pairs with product `|lambda|^2 > 0`. Then `a_i/D > 0` iff `a_i > 0`, and `o_k/D < 0` iff `o_k < 0`. ∎

Consequence: `a` and `o_k` do not involve `t`. Let

    B+ := { beta : a_i(beta) > 0 for all i,  o_k(beta) < 0 for all k }.

In the `(beta, t)` quarter-plane, the SFP set is `(B+ x (0,inf)) ∩ {p Hurwitz} ∩ legal`. So permission and the
off-support condition cut out vertical strips in beta, and only stability depends on `t`.

---

## 2. Lemma B (boundary set at fixed t)

### B.1 Facts about Hurwitz polynomials

Let `f(z) = c_0 z^n + ... + c_n` be real, with `c_0 > 0` and `n >= 1`.

- **(b0) CITED: the Routh-Hurwitz theorem.** `f` is Hurwitz iff `Delta_1(f) > 0, ..., Delta_n(f) > 0`.
  Sources: A. Hurwitz, *Ueber die Bedingungen, unter welchen eine Gleichung nur Wurzeln mit negativen reellen Theilen
  besitzt*, Math. Ann. 46 (1895) 273-284; F. R. Gantmacher, *The Theory of Matrices*, Vol. II, Chelsea 1959, Ch. XV.
  The Gantmacher chapter is certain; I quote its section number from memory (§6), so verify it before citing externally.
  We use both directions. "<=" is what makes `ctln.hurwitz_stable` returning True sound. "=>" is what makes a zero
  or negative minor a proof that `f` is not Hurwitz.
- **(b1) PROVED.** If `f` is Hurwitz, all `c_l > 0`.
  *Proof.* `f = c_0 prod (z + r) prod (z^2 + 2uz + u^2 + v^2)` with `r, u > 0` (real roots `-r`, pairs `-u ± iv`).
  Every factor has positive coefficients. ∎
- **(b2) PROVED.** `Delta_n(f) = c_n Delta_(n-1)(f)`.
  *Proof.* Row n of the Hurwitz matrix is `(c_(2n-1), ..., c_n) = (0, ..., 0, c_n)`; expand along it. ∎
- **(b3) PROVED, self-contained. This is the consequence of Orlando's formula that we actually use.** If `n >= 2` and
  `f(iw) = 0` for some real `w != 0`, then `Delta_(n-1)(f) = 0`.
  *Proof.*
  1. Put `P(x) := sum_(l=0..n) c_l x^l = x^n f(1/x)`, so `P(0) = c_0 != 0`.
  2. `f` is real, so `f(±iw) = 0`, hence `P(∓i/w) = 0`. So the real quadratic `1 + w^2 x^2` divides `P` in `R[x]`:
     `P = (1 + w^2 x^2) R` with `R` real, `deg R <= n-2`, and `R(0) = c_0 != 0`.
  3. Let `V(x) := x R(-x) = sum_(j=1..n-1) v_j x^j`. Then `V != 0`.
  4. `P(x)V(x) = (1 + w^2 x^2) * x * R(x)R(-x)` is odd, because `R(x)R(-x)` is even. So every even-degree coefficient of `PV` is 0.
  5. For `i = 1..n-1`, `[x^(2i)] (P V) = sum_(j=1..n-1) c_(2i-j) v_j = (H_(n-1) v)_i`, where
     `H_(n-1) = [c_(2i-j)]_(i,j=1..n-1)` is the matrix whose determinant is `Delta_(n-1)`.
  6. So `H_(n-1) v = 0` with `v != 0`, and `Delta_(n-1)(f) = 0`. ∎
- **(b4) CITED for context: Orlando's formula.** If `f` has roots `z_1..z_n`, then
  `Delta_(n-1)(f) = (-1)^(n(n-1)/2) c_0^(n-1) prod_(1<=i<l<=n) (z_i + z_l)`.
  Sources: L. Orlando, *Sul problema di Hurwitz relativo alle parti reali delle radici di un'equazione algebrica*,
  Math. Ann. 71 (1911) 233-245 (page numbers from memory); Gantmacher Vol. II Ch. XV (§7 from memory).
  (b3) is its special case `z_i + z_l = 0` with `z_i = iw`. Only C'(i) uses (b4), and only to explain spurious critical values.
  Sanity-checked: the symbolic identity holds for `n = 2..5`, and 60 exact random instances hold for `n = 6, 7`.

### B.2 Continuity of roots

**CITED: Rouché's theorem** (e.g. Ahlfors, *Complex Analysis*, 3rd ed., McGraw-Hill 1979). Let `f_beta` be real polynomials
of fixed degree `n`, with coefficients continuous in `beta` and leading coefficient never 0. Let `z_0` be a root of
`f_(beta_0)` of multiplicity `r`, and suppose the closed disk `B(z_0, rho)` contains no other root of `f_(beta_0)`. Then for
`beta` close to `beta_0`, `f_beta` has exactly `r` roots (with multiplicity) in `B(z_0, rho)`.
The proof: on the circle `|z - z_0| = rho` we have `|f_(beta_0)| >= c > 0`, while `sup |f_beta - f_(beta_0)| -> 0`.

### Lemma B

Fix `t in (0, inf)`. Rationality of `t` matters only for computation. Let `L_t := (2 + 1/t, inf)` be the legal beta-range.
Let `p_t(mu; beta) := det(mu I + K) = phi + t omega`, and let `P_t` be the finite family of polynomials in beta:

    a_1, ..., a_m;   D_t := d + t s;   o_k (k outside sigma);   Delta_(m-1)(p_t)   (for m = 1 use Delta_0 := 1).

- **(a)** `p_t` is monic of degree m in mu, its coefficients are polynomials in beta, and `p_t(0) = D_t`.
- **(b)** (CITED, (b0)) If `p_t(.; beta)` is Hurwitz, then `Delta_1, ..., Delta_m > 0` at beta.
- **(c)** If `p_t(.; beta_0)` has a root on the imaginary axis, then `D_t(beta_0) = 0` or `Delta_(m-1)(beta_0) = 0`.
- **(d)** If some member of `P_t` is the zero polynomial, SFP(sigma) is false at every beta in `L_t`.
  Otherwise let `Z_t := { real roots in L_t of members of P_t }`, a finite set. Then:
  - SFP(sigma) is **constant on each connected component** (open interval) of `L_t \ Z_t`;
  - SFP(sigma) is **false at every point of `Z_t`**.

  Hence `{ beta in L_t : SFP }` is a finite union of components, i.e. of open intervals with endpoints in
  `Z_t ∪ {2 + 1/t, +inf}`.
- **(e)** For any `beta_0` in a component `I`, the exact verdict at `beta_0` is the verdict on all of `I`. In particular,
  suppose that at `beta_0` some coefficient of `p_t` is `<= 0`, or some `Delta_k(p_t)(beta_0) <= 0` for any k. This
  includes a **vanishing minor other than `Delta_(m-1)`** at a point where `D_t != 0` and `Delta_(m-1) != 0`. Then `p_t`
  is not Hurwitz at `beta_0`, and SFP is false on all of `I`.

*Proof.*

- (a) Leibniz expansion, and `p_t(0) = det K = D_t` by A(iii).
- (b) Is (b0).
- (c) `p_t` is real, so a root on `iR` is either 0 or one of a pair `±iw` with `w != 0`. A zero root means `D_t = p_t(0) = 0`.
  A pair gives `Delta_(m-1) = 0` by (b3) when `m >= 2`. For `m = 1` the only possible root on `iR` is 0.
- (d) *Zero polynomials.*
  - `a_i ≡ 0`: wherever `D_t != 0`, `x_i = 0` by A(v), so (1) fails; where `D_t = 0`, (3) fails.
  - `D_t ≡ 0`: (3) fails.
  - `o_k ≡ 0`: wherever `D_t != 0`, `y_k = 0`, so the strict (2) fails.
  - `Delta_(m-1) ≡ 0`: by (b), `p_t` is never Hurwitz.

  In every case SFP is false (Corollary A.1).

  *Constancy.* Let `I` be a component of `L_t \ Z_t`. Each of `a_i`, `D_t`, `o_k` is continuous and has no zero on the
  interval `I`, so its sign is constant there. So the first three clauses of Corollary A.1 have constant truth value on `I`.
  By (c), `p_t(.; beta)` has no root on `iR` for beta in `I`. Let `n_+(beta)` be the number of roots with positive real part.
  1. Fix `beta_0` in `I`. Choose `rho > 0` smaller than half the distance between distinct roots of `p_t(.; beta_0)`, and
     smaller than the distance of each root to `iR`.
  2. By B.2, for beta near `beta_0` each disk `B(z_l, rho)` contains exactly `mult(z_l)` roots of `p_t(.; beta)`. These
     account for all `m` roots, and each disk lies in one open half-plane.
  3. So `n_+(beta) = n_+(beta_0)` near `beta_0`. Thus `n_+` is locally constant on the connected set `I`, hence constant.
  4. No root is on `iR`, so "`p_t` Hurwitz" holds iff `n_+ = 0`, which is constant on `I`.

  With Corollary A.1, SFP is constant on `I`.

  *Points of `Z_t`.* Let `beta_0 in Z_t`.
  - `D_t(beta_0) = 0`: (3) fails.
  - `a_i(beta_0) = 0` with `D_t != 0`: `x_i = 0`, so (1) fails.
  - `o_k(beta_0) = 0` with `D_t != 0`: `y_k = 0`, violating the strict (2).
  - `Delta_(m-1)(beta_0) = 0` with the others nonzero: not Hurwitz by (b).

  So SFP is false on `Z_t`. ∎
- (e) The first sentence is (d). For the second: by (b1) and (b), the stated condition implies `p_t` is not Hurwitz at
  `beta_0`. So SFP(`beta_0`) is false, and therefore false on `I` by (d). ∎

Remarks.

- (e) is a mathematical fact about open cells; it does not contradict CLAUDE.md §3. At an arbitrary sampled parameter a
  vanishing minor says nothing about the neighbourhood, which is why a pointwise sweep must call it *undecided*. Inside a
  cell of the decomposition the neighbourhood question has already been settled by (d). Whether code should use this is a
  policy decision; fastregion.py currently re-samples instead, which is conservative. §3 is unchanged by this note.
- Points of `Z_t` are never SFP under conditions (1)-(3) as stated, which are all strict. At such a point a fixed point, if
  one exists, is degenerate or non-hyperbolic. That is the precise meaning of "boundary" in the region computation.
- The same lemma holds with the roles swapped: fix beta and vary `t` over `(1/(beta-2), inf)`. Then `a_i` and `o_k` are
  constants, `D = d + t s` is linear in `t`, and the critical set is the real roots of `D` and of `Delta_(m-1)(p_t)` viewed
  as a polynomial in `t`. The proof is word for word the same. This is how the threshold `t* ~ 52.93` at `beta = 54/25`
  for Geneson's graph arises.

**Lemma B-inf (the t = inf problem).** For `beta > 2`, define

    SFP_inf(beta)  :<=>  s(beta) != 0,  a_i/s > 0 for all i,  o_k/s < 0 for all k,  and  omega(.; beta) Hurwitz.

Let `P_inf := {a_i, s, o_k, Delta_(m-2)(omega)}`. For `m = 2` use `Delta_0 := 1`. For `m = 1`, `omega = 1` has no roots
and counts as Hurwitz.

Then (c), (d) and (e) hold verbatim on `(2, inf)`, with `omega`, `s` and `Delta_(m-2)(omega)` in place of `p_t`, `D_t` and
`Delta_(m-1)(p_t)`.
*Proof.* `omega` has degree `m-1` in mu, constant leading coefficient `m > 0` (A(iv)), and constant term `s`. The proof of
(c)-(e) used only these properties. ∎ What SFP_inf means for the dynamics is Lemma C.

---

## 3. Lemma C (eps -> 0, i.e. t -> inf, at fixed beta)

**Lemma C.** Fix `beta > 1` and sigma with `m >= 1`. Let `z_1..z_r` be the distinct roots of `omega(.; beta)`, with
multiplicities `m_1..m_r` (summing to `m-1`).

- **(a) Spectrum.** For every `rho > 0` there is `T` such that for all `t > T`:
  - for each `l`, exactly `m_l` eigenvalues of `K(beta, t)` (with multiplicity) lie within `rho` of `-z_l`;
  - the remaining eigenvalue `lambda_m(t)` is real, with `lambda_m(t) = m t + O(1)`. By C'(c), more precisely
    `lambda_m(t) = m t + 1 + beta e/m + o(1)`.
- **(b)** Suppose `omega(.; beta)` is Hurwitz, all `a_i(beta)` are nonzero with a common sign, and `o_k(beta) s(beta) < 0`
  for every `k` outside sigma. Then there is `T(beta)` such that for every `t > T(beta)`, SFP(sigma) holds at
  `eps = 1/(1+t)`, `delta = (beta-1)/(1+t)`, and these parameters are legal if `beta > 2`. By C'(h), the common sign is
  then necessarily `+`, and `s > 0`.
- **(c)** If `omega(.; beta)` has a root with positive real part, then there is `T` such that `K(beta, t)` is not stable
  for all `t > T`, so SFP fails.
- **(d) Positivity in the limit.**
  - If `s(beta) != 0`, then for `t > |d(beta)/s(beta)|` we have `sign D = sign s`, and:
    - (1) holds iff `a_i(beta)/s(beta) > 0` for all i, iff all `a_i(beta)` are nonzero with one common sign;
    - (2) holds iff `o_k(beta) s(beta) < 0` for all `k` outside sigma.
  - If `s(beta) = 0`, then (1) fails for every `t` with `D != 0`. If moreover `d(beta) = 0`, then `D ≡ 0` and (3) fails.

  So at every beta: "(1) holds for all large `t`" iff "all `a_i(beta)` are nonzero with a common sign".

*Proof.*

- (a) For `t > 0`, the roots of `p_t = phi + t omega` are those of `g_t := omega + phi/t`.
  1. Choose `rho' <= rho` so the closed disks `B_l := B(z_l, rho')` are pairwise disjoint. On their boundary circles,
     `|omega| >= c > 0` and `|phi| <= M`.
  2. For `t > M/c`, `|phi/t| < |omega|` on the circles, so Rouché gives exactly `m_l` roots of `g_t` in `B_l`. That
     accounts for `m-1` roots.
  3. `g_t` has degree `m` (leading coefficient `1/t`), so exactly one root `zeta(t)` lies outside the union of the disks.
  4. The roots of the real polynomial `omega` are closed under conjugation and the radii are equal, so the union of disks
     is conjugation-invariant. Hence `conj(zeta(t))` is also a root outside the union; since there is only one such
     root, `zeta(t)` is real.
  5. The coefficient of `mu^(m-1)` in `p_t` is `tr K = m + m t`, because `C` has unit diagonal. So the roots sum to
     `-m(1+t)`, and `zeta(t) = -m(1+t) - (sum of the m-1 roots in the disks) = -m t + O(1)`.
  6. The eigenvalues of `K` are the negatives of the roots of `p_t`. Letting `rho' -> 0` gives convergence with multiplicity. ∎
- (b) Legality: `t(beta-2) > 1` once `t > 1/(beta-2)`.
  Positivity and off-support: `s = sum a_i` has the common sign of the `a_i`, so `s != 0` and (d) applies.
  Stability: `omega` is Hurwitz, so choose `rho' < min_l |Re z_l|`; then every `B_l` lies in `Re < 0`. For `t > M/c`, the
  `m-1` roots in the disks have `Re < 0`, and `Re zeta(t) <= -m(1+t) + (m-1)(max_l |z_l| + rho') < 0` for large `t`.
  So `p_t` is Hurwitz, and Corollary A.1 gives SFP. Take `T(beta)` to be the largest of these thresholds. ∎
- (c) Take a disk around a root with `Re z > 0` lying inside `Re > 0`. For large `t`, `p_t` has a root in it, i.e. `K`
  has an eigenvalue with negative real part. ∎
- (d) For `t > |d/s|`, `D = d + t s` has the sign of `s`; then apply A(v). If all `a_i` share a nonzero sign, then `s != 0`
  with that sign, so `a_i/s > 0`. Conversely, `a_i/s > 0` for all i forces a common sign. If `s = 0`, then `D = d` for
  all `t` and `sum_i a_i = 0`, so the `a_i` cannot all have the sign of `d`. ∎

Caveats:

- **Uniformity.** `T(beta)` depends on beta and in general blows up as beta approaches a root of a member of `P_inf`. So
  Lemma C certifies, for each beta in an SFP_inf cell, a stable fixed point for all `eps < 1/(1+T(beta))`. It gives **no
  explicit eps0**. An explicit finite-eps certificate needs a finite-t region computation (Lemma B), as the literature
  agent did for Geneson's graph (`t* ~ 52.93` at `beta = 54/25`).
- **One-sided negatives.** "SFP_inf false" at beta only says SFP fails for all large `t`. It says nothing about moderate `t`.
- **Boundary values.** At a beta where `omega` has a root on `iR` (`s = 0` or `Delta_(m-2)(omega) = 0`, by (b3)), Lemma C decides nothing.

---

## 4. Lemma C' and Lemma I: the t = inf problem in closed form

This section proves what `notes/limit_spectral.md` calls Lemma I (its Claims 1-3). Claim 3 holds in a stronger form than
stated there: no nondegeneracy hypothesis is needed.

**Lemma C' (compression).** Let `m >= 2`. Let `U in R^(m x (m-1))` have orthonormal columns spanning `1^perp`, and put
`G := U^T N U`. For any `X`, write `X|perp := U^T X U`.

- **(a)** `omega(mu; beta) = m det((mu + 1) I_(m-1) + beta G)`, identically in mu and beta.
  Equivalently, `(1/m) 1^T adj(mu I + C) 1 = det(mu I + C|perp)`; this is the "classical identity" quoted in limit_spectral.md.
- **(b)** `s(beta) = m det(I + beta G) = m prod_(nu in spec G) (1 + beta nu)`.
- **(c)** `tr G = -e/m`, so the large eigenvalue in C(a) is `m t + 1 + beta e/m + o(1)`.
- **(d)** For `beta > 0`: `omega(.; beta)` is Hurwitz iff `1 + beta Re nu > 0` for every `nu in spec G`, iff
  `beta < beta_c := -1 / min Re spec(G)`. If `min Re spec G >= 0`, set `beta_c := +inf`.
- **(e)** On the t = inf legal range `beta > 2`, the stability clause of SFP_inf holds exactly on the interval
  `(2, beta_c)`. This interval is nonempty iff `omega(.; 2)` is Hurwitz, iff `min Re spec(G) > -1/2`.
- **(f)** Necessary condition: if `omega(.; beta)` is Hurwitz for some `beta > 0`, then `e < m(m-1)/beta`. In particular,
  t = inf stability at a legal `beta > 2` needs `e < m(m-1)/2`, i.e. more than half of the ordered pairs in sigma are
  edges. In the limit it also needs `q = beta - 1 < m(m-1)/e - 1`.
- **(g) Rational form.** Let `V in Q^(m x (m-1))` have columns forming a basis of `1^perp` (e.g. the columns `e_i - e_m`).
  Then `G_V := (V^T V)^-1 V^T N V` is similar to `G`, so `spec G`, `beta_c` and the test in (e) are exact over `Q`.
- **(h)** If `omega(.; beta)` is Hurwitz, then `s(beta) > 0`. Hence SFP_inf(beta) holds iff `a(beta) > 0`,
  `o_k(beta) < 0` for all k, and `beta < beta_c`.
- **(i)** (Uses (b4), CITED.) `Delta_(m-2)(omega) = ± m^(m-2) prod_(i<j) (2 + beta(nu_i + nu_j))`. So its real roots are
  `beta = -2/(nu_i + nu_j)` over index pairs with a real sum. Conjugate pairs give the genuine crossings `-1/Re nu`.
  Pairs of real eigenvalues give **spurious** critical values, where no root crosses.

*Proof.*

- (a)
  1. `Q := [u | U]` with `u := 1/sqrt(m)` is orthogonal, and `Q^T 1 = sqrt(m) e_1`.
  2. So `Q^T (mu I + C + t 1 1^T) Q = X + t m e_1 e_1^T`, where `X := mu I + Q^T C Q`.
  3. The determinant is affine in the (1,1) entry, with slope equal to the complementary minor `det(mu I_(m-1) + U^T C U)`.
     Hence `det(mu I + K) = det(mu I + C) + t m det(mu I_(m-1) + U^T C U)` for all `t`.
  4. Comparing the coefficient of `t` with A(iv) gives `omega = m det(mu I_(m-1) + U^T C U)`, and
     `U^T C U = U^T U + beta U^T N U = I + beta G`. ∎
- (b) Put `mu = 0` in (a); A(iv) gives `omega(0) = s`. ∎
- (c) `tr(U^T N U) = tr(N U U^T) = tr(N (I - J/m)) = tr N - (1^T N 1)/m = -e/m`.
  In C(a), the `m-1` bounded eigenvalues converge to `-(roots of omega) = 1 + beta nu` by (a), with sum
  `(m-1) + beta tr G = m - 1 - beta e/m`. The total is `tr K = m + m t`. ∎
- (d) By (a), the roots of `omega` are `mu = -1 - beta nu`, `nu in spec G`, with multiplicity. `Re mu < 0` iff
  `1 + beta Re nu > 0`. For `beta > 0` this holds for all `nu` iff `Re nu > -1/beta` for all `nu`. That is equivalent to
  `beta < -1/min Re spec G` when this minimum is negative, and it always holds otherwise. ∎
- (e) Immediate from (d): the condition `Re nu > -1/beta` gets harder as beta grows. At `beta = 2` it reads
  `min Re spec G > -1/2`. ∎
- (f) Hurwitz at beta gives `Re nu > -1/beta` for all `m-1` eigenvalues. Summing,
  `-e/m = tr G = sum Re nu > -(m-1)/beta`. ∎
- (g) `V = U S` with `S := U^T V` invertible. Then `V^T V = S^T S` and `V^T N V = S^T G S`, so `G_V = S^-1 G S`. ∎
- (h) `s = omega(0)` is the constant coefficient of a Hurwitz polynomial with positive leading coefficient, so `s > 0`
  by (b1). Then `a_i/s > 0` iff `a_i > 0`, and `o_k/s < 0` iff `o_k < 0`. ∎
- (i) Apply (b4) to `omega` (leading coefficient `m`, degree `m-1`, roots `mu_i = -1 - beta nu_i`). Note
  `mu_i + mu_j = -(2 + beta(nu_i + nu_j))`. ∎

**Lemma I (limit form of all three conditions; proves limit_spectral.md Claims 1-3).** Let `beta > 1` and
`gamma := 1 - 1/beta`, so `gamma in (0,1)`. For a legal `beta > 2`, `gamma = q/(q+1) in (1/2, 1)` with `q = beta - 1` the
limiting q. Write `A^T` for `A_sigma^T`, so `N = J - I - A^T` and `C = -beta (A^T + gamma I) + beta J`.

- **(I.1) = Claim 1.** As `t -> inf`, one eigenvalue of `K` is `m t + O(1)` and the other `m-1` converge to the
  eigenvalues of `C|perp`. *Proof:* C(a) and C'(a). ∎
- **(I.2) = Claim 2.** `C|perp = (1-beta) I - beta A^T|perp`. Further, `omega(.; beta)` is Hurwitz iff every eigenvalue
  of `A|perp` has `Re < -gamma`, iff `alpha := max Re spec(A|perp) < -gamma`. The limit-stable range in beta is
  `(2, beta_c)` with `beta_c = 1/(1 + alpha)`. Always `alpha >= -1`, and `beta_c = +inf` when `alpha = -1`. The range is
  nonempty iff `alpha < -1/2`.
  *Proof.*
  1. `U^T J U = 0`, so `C|perp = (1-beta) I - beta (A^T)|perp` and `G = -I - (A^T)|perp`.
  2. `(A^T)|perp = (A|perp)^T` has the same spectrum as `A|perp`. So `spec G = {-1 - lambda : lambda in spec A|perp}`.
  3. C'(d): `1 + beta Re nu > 0` iff `1 - beta - beta Re lambda > 0` iff `Re lambda < -gamma`.
  4. `alpha >= -1`: the mean of the `m-1` eigenvalues is `tr(A|perp)/(m-1) = -|E|/(m(m-1)) >= -1`, where `|E|` is the
     number of directed edges inside sigma.
  5. `min Re spec G = -1 - alpha`, so `beta_c = 1/(1+alpha)`. ∎
- **(I.3) = Claim 3, strengthened.** For every `beta > 1`, with no nondegeneracy assumption:

      a(beta) is one-signed (all a_i nonzero, common sign)   <=>   A^T + gamma I is nonsingular and w := (A^T + gamma I)^-1 1 > 0.

  In that case `a = ((beta s - d)/beta) w`, and `beta s - d` has the sign of `a`. Note that `w` solves
  `gamma w_i + sum_(j in sigma, j -> i) w_j = 1` for every `i in sigma`, and that a one-signed solution `w` is
  automatically positive.
  *Proof.*
  1. Since `C adj(C) = d I`, we have `C a = d 1`. Substituting `C = -beta(A^T + gamma I) + beta J` gives

         (A^T + gamma I) a = ((beta s - d)/beta) 1.        (*)

  2. Apply (T2) with `M = C`, `tau = -beta`: `det(C - beta J) = d - beta s`. Hence
     `det(A^T + gamma I) = (-1/beta)^m (d - beta s)`.
  3. (=>) Suppose `a` is one-signed with sign `sigma_a`. Row i of the left side of (*) is
     `gamma a_i + sum_(j -> i) a_j`. This is nonzero with sign `sigma_a`, because `gamma > 0` and `A >= 0`. So
     `beta s - d != 0` with sign `sigma_a`. By step 2, `A^T + gamma I` is nonsingular, and `w = beta a/(beta s - d) > 0`.
  4. (<=) A(ii) holds for every value of `t`, including `t = -beta`, so `adj(C - beta J) 1 = adj(C) 1 = a`. Since
     `C - beta J = -beta (A^T + gamma I)` is nonsingular,
     `a = det(C - beta J) (C - beta J)^-1 1 = (d - beta s)(-1/beta) w`, a nonzero multiple of `w > 0`.
  5. Automatic positivity: if all `w_j < 0`, then row i of `(A^T + gamma I) w = 1` has a negative left side. ∎

  Claim 3 in limit_spectral.md assumes "C invertible and `1 - beta 1^T u != 0`". Neither assumption is needed: the
  determinant identity in step 2 shows that one-signedness of `a` already forces `A^T + gamma I` to be nonsingular,
  and `det C` never enters.
- **(I.4) Off-support in the limit (new; it completes the triple).** Suppose (I.3) holds, i.e. `w > 0`. For `k` outside
  sigma put `in_k(w) := sum_(b : sigma_b -> k) w_b`. Then

      o_k s < 0   <=>   in_k(w) < 1.

  *Proof.*
  1. Let `c := (beta s - d)/beta`, which is nonzero, so `a = c w`. Put `W1 := 1^T w > 0`. Then `s = c W1`, so
     `sign s = sign c`.
  2. `d = beta s - beta c = beta c (W1 - 1)`.
  3. `o_k = d - beta c sum_b N_kb w_b = beta c (W1 - 1 - sum_b N_kb w_b) = beta c (in_k(w) - 1)`, using
     `W1 - sum_b N_kb w_b = sum_b (1 - N_kb) w_b = in_k(w)`.
  4. So `o_k s = beta c^2 W1 (in_k(w) - 1)`. ∎

**Summary (limit graph rule).** For `beta > 2` with `gamma = 1 - 1/beta`:

    SFP_inf(beta)  <=>  (I1) A_sigma^T + gamma I nonsingular and w := (A_sigma^T + gamma I)^-1 1 > 0
                        (I2) sum_(b : sigma_b -> k) w_b < 1 for every k outside sigma
                        (I3) every eigenvalue of A_sigma|perp has real part < -gamma

There is no exceptional set of beta. In words:

- inside sigma, the self term plus the in-flow is exactly 1: `gamma w_i + in_i(w) = 1`;
- every outside node receives strictly less than 1;
- the compressed adjacency has spectral abscissa below `-gamma`.

For a clique, `w = 1/(m-1+gamma)`, and (I2) says exactly "no target" (cf. Lemma 0), as it should.
*Proof.* Combine C'(h), (I.2), (I.3) and (I.4). (I1) gives `s != 0`, because `a = c w` with `c != 0`. ∎

**Geneson's graph (full support, m = 6, e = 11, |E| = 19).** This is an exact computation; the code is in the Appendix.

- `charpoly(G_V) = z (6z^4 + 11z^3 + 16z^2 + 10z + 3)/6`.
- The zero eigenvalue comes from paper nodes 2 and 3 having identical `N`-rows: `(e_2 - e_3)^T N = 0` with
  `e_2 - e_3 in 1^perp`.
- The polynomial whose roots are `-2/(nu_i + nu_j)` (i < j, `nu` roots of the quartic) is computed by resultant. It
  equals **exactly** (ratio 1) `797 b^6 - 6468 b^5 + 22544 b^4 - 44440 b^3 + 53280 b^2 - 38016 b + 13824`. This is the
  irreducible sextic found independently by the literature agent and by fastregion.
- So `beta_c = -1/Re nu_1`, with `nu_1 ~ -0.461605 ± 1.130501 i` (decimals for orientation only), and `alpha = -0.538395`.
  This explains the "degree 6 already at n = 6" anomaly: `Re nu` for `nu` a root of a quartic has degree `<= C(4,2) = 6`.
- The second sextic root above 2, `~2.19751 = -1/Re nu_2` with `nu_2 ~ -0.455062 ± 0.358105 i`, is the second pair
  crossing. fastregion also reports it, isolated in `(68/31, 11/5)`.
- C'(f): `e = 11 < 15`, as required.

---

## 5. Lemma D (graphical domination kills every fixed point through the dominated node)

**Definition** (Curto, arXiv:2510.05098, Def. 1.1, as quoted verbatim in graphs.py by the literature agent). Node `k`
*graphically dominates* node `j` (`k != j`) iff:

- (i) for every `i` not in `{j,k}`, `i -> j` implies `i -> k`; and
- (ii) `j -> k` and `k -/-> j`.

`graphs.graphically_dominates(A, k, j)` implements exactly (i) and (ii): `A[j][k] = 1`, `A[k][j] = 0`, and no `i` outside
`{j,k}` with `A[i][j] = 1` and `A[i][k] = 0`. I checked the code.

**Hypotheses for Lemma D.**

- *Used:* `eps > 0`, `delta >= 0`, and the same external input at `j` and `k` (`theta_j <= theta_k` suffices).
- *Not used:* `eps < 1`, legality `q > 1`, `theta > 0`, nondegeneracy, and anything about edges not mentioned in (i) and (ii).

**D.1 (comparison inequality).** If `k` dominates `j`, then for every `x in R^n` with `x >= 0`:

    y_j(x) - y_k(x) <= (1 - eps) x_j - (1 + delta) x_k,
    equivalently  (-x_j + y_j(x)) - (-x_k + y_k(x)) <= -eps x_j - delta x_k.

*Proof.* `y_j - y_k = (theta_j - theta_k) + sum_i (W[j][i] - W[k][i]) x_i`, and `theta_j - theta_k <= 0`. Term by term:

- `i = j`: `W[j][j] = 0`, and `W[k][j]` reads `A[j][k] = 1` by (ii), so it is `-1 + eps`. The coefficient is `1 - eps`.
- `i = k`: `W[j][k]` reads `A[k][j] = 0` by (ii), so it is `-1 - delta`, and `W[k][k] = 0`. The coefficient is `-1 - delta`.
- `i` not in `{j,k}`, `i -> j`: then `i -> k` by (i), so both entries are `-1+eps` and the difference is 0.
- `i` not in `{j,k}`, `i -/-> j`: then `W[j][i] = -1-delta <= W[k][i]`, since `W[k][i]` is `-1-delta` or `-1+eps` and
  `eps + delta >= 0`. The difference is `<= 0`, and `x_i >= 0`. ∎

**Lemma D.** If some `k` dominates `j`, then no fixed point `x*` of `dx/dt = -x + [W x + theta]_+` has `x*_j > 0`.
*Proof.* A fixed point satisfies `x* = [W x* + theta]_+`, so `x* >= 0` and `x*_l = [y_l(x*)]_+` for all `l`. Suppose
`x*_j > 0`; then `y_j(x*) = x*_j`. By D.1,

    y_k(x*) >= y_j - (1-eps) x*_j + (1+delta) x*_k = eps x*_j + (1+delta) x*_k > 0,

so `x*_k = y_k(x*)`. Then `x*_k >= eps x*_j + (1+delta) x*_k`, i.e. `0 >= eps x*_j + delta x*_k >= eps x*_j > 0`.
Contradiction. ∎

**Which hypotheses matter.**

- **eps < 1 is NOT needed** (contrary to the task sketch). The proof never uses the sign of `-1 + eps`. It uses
  `eps > 0` (final strict inequality), `delta >= 0` (final step), and `eps + delta >= 0` (D.1). Sanity check: 14,328
  supports containing a dominated node, with eps drawn from `(0,1)`, `{1}` and `(1,3]`; no fixed point was found.
- **The input condition IS needed.** Take `A` with the single edge `j -> k`, so `k` dominates `j`. With
  `theta_j = 10, theta_k = 1, eps = 1/4, delta = 1/2`, the point `x_j = 10, x_k = 0` is a stable fixed point, since
  `y_k = 1 - (3/4) 10 < 0`. Checked exactly.
- **Both halves of (ii) are needed.**
  - If `k -> j`, the `x_k` coefficient in D.1 becomes `-1 + eps`, and the argument only yields `x_k >= x_j`.
  - If `j -/-> k`, the `x_j` coefficient becomes `1 + delta`, and the argument only yields `y_k >= -delta x_j + ...`.

**Corollary D.2 (supports through j).** If `j` is dominated in `G` and `j in sigma`, then (1) and even the weak form of (2)
(`<= 0`) cannot both hold: together they would make `x^` a fixed point with `x^_j > 0`. So sigma is not a fixed-point
support, stable or not.

**Corollary D.3 (full-support prune).** Suppose `G` has a dominated node. Then:

- for every `eps > 0, delta >= 0` (in particular every legal parameter) there is no `x > 0` with `(I - W) x = theta 1`, so
  `[n]` is not an SFP support;
- `G` is not complete bidirected, because `k -/-> j`.

So the full-support instance of the conjecture holds trivially for `G`, and **the full-support sweep may be restricted to
domination-free graphs** without losing any counterexample. The certificate for a skipped graph is the witness pair
`(k, j)`, checkable in `O(n)`.

**Corollary D.4 (induced-subgraph prune, any support size).** If the induced subgraph `G[sigma]` has a node `j`
dominated by some `k` within `G[sigma]` (condition (i) quantified over `i in sigma` only), then (1) fails for sigma in `G`.
*Proof.* `W_sigma` is the weight matrix of the CTLN of `G[sigma]`, and (1) for sigma says exactly that `G[sigma]` has a
full-support fixed point. Apply D.3 to `G[sigma]`. ∎ Domination inside `G[sigma]` is weaker than in `G` (fewer `i` to
check), so this prunes strictly more than D.2.

---

## 6. Lemma 0, Lemma E and the Stage 4 prune

**Lemma 0 (cliques).** Under (H), let sigma be a clique. Then (1) and (3) hold, and (2) holds iff sigma has no target.
Hence a clique is an SFP support iff it is target-free. So **at given parameters, the conjecture is equivalent to
"every SFP support is a clique"**.
*Proof.*

1. `N = 0`, so `C = I`, `a = 1`, `d = 1`, `s = m`, and `D = 1 + t m > 0`.
2. `p = (mu+1)^(m-1) (mu + 1 + t m)` is Hurwitz, and `x = (theta/eps) 1/(1 + t m) > 0`.
3. For `k` outside sigma, `o_k = 1 - beta r_k` with `r_k := #{b : sigma_b -/-> k}`, so
   `y_k = theta (1 - beta r_k)/(1 + t m)`.
4. If `k` is a target, then `r_k = 0` and `y_k > 0`. Otherwise `r_k >= 1` and `1 - beta r_k <= 1 - beta = -delta/eps < 0`. ∎

This matches Geneson's clique formula `theta[eps - r(eps+delta)]/(m - (m-1)eps)`, after multiplying the numerator and
denominator by `eps`.

**Lemma E.** Assume (H); in fact `eps > 0`, `delta > 0` and uniform theta suffice. Let `k` dominate `j` in `G`, and let
`G' := G - j`, with nodes keeping their names. Then for every `sigma ⊆ V(G) \ {j}`:

    SFP(sigma) in G   <=>   SFP(sigma) in G'.

The same holds with (2) weakened to `<= 0`.

*Proof.*

1. `W_sigma` is the same in `G` and `G'`: its entries read `A[b][a]` for `a, b in sigma`, and `G'` is an induced
   subgraph. So (1) and (3) are the same statements, with the same solution `x`.
2. For `l` outside sigma with `l != j`, `y_l(x^) = theta + sum_(b in sigma) W[l][b] x_b` is the same number in `G` and
   `G'`. The term for `j` contributes nothing since `x^_j = 0`, and `W[l][b]` reads `A[b][l]`, which is unchanged.
   So "(2) in G" holds iff "(2) in G'" holds and `y_j(x^) < 0`.
3. The direction G => G' is now immediate.
4. For G' => G, let `x` satisfy (1), (2) and (3) in `G'`, and let `x^` be its lift to `V(G)` (so `x^_j = 0`). Apply D.1 in
   `G` at `x^ >= 0`:
   - **Case k in sigma.** `x^_k = x_k > 0` and `y_k(x^) = x_k` (row k of (1)). D.1 with `x^_j = 0` gives
     `y_j <= y_k - (1+delta) x_k = -delta x_k < 0`.
   - **Case k not in sigma.** `k != j`, so `k` is a node of `G'` outside sigma, and (2) in `G'` gives `y_k(x^) < 0`.
     By step 2 this is the same value in `G`. D.1 with `x^_j = x^_k = 0` gives `y_j <= y_k < 0`.

   In both cases `y_j(x^) < 0`, so (2) holds in `G`.
5. Weak version: case 1 still gives `< 0`, and case 2 gives `y_j <= y_k <= 0`. ∎

Strictness: case 1 needs `delta > 0`; with `delta = 0` it gives only `y_j <= 0`. Case 2 inherits strictness from `k`.
The proof uses only D.1, not Theorem 1 of arXiv:2510.05098.

**Corollary E.1.** With D.2: the SFP supports of `G` are exactly the SFP supports of `G'`, and none of them contains `j`.

**Corollary E.2 (counterexamples transfer).** Whether `sigma ⊆ V \ {j}` is a clique is the same in `G` and `G'`.
So, by Lemma 0 and E.1, at every parameter satisfying (H): sigma is a counterexample (a non-clique SFP support) for `G`
iff it is one for `G'`. Target-freeness itself can differ, since `j` may be a target in `G` but cannot be one in `G'`.
Lemma 0 shows that target-freeness only matters for cliques, and for cliques it coincides with SFP in both graphs.

**Corollary E.3 (automatic off-support condition).** If `l` outside sigma is dominated in `G` by some `k in sigma`, then
(1) for sigma implies `y_l(x^) < 0`. So condition (2) at `l` need not be checked. (This is case 1 of the proof; it uses
only (1).)

**Corollary E.4 (Stage 4 at n = 6, |sigma| = 5, sigma = V \ {l}).**

- (a) If `G[sigma]` has a dominated node (within `G[sigma]`), then sigma fails (1). Skip it (D.4).
- (b) Otherwise, if `l` is dominated in `G` (necessarily by some `k in sigma`), then SFP(sigma in G) holds iff the 5-node
  graph `G[sigma] = G - l` has an SFP full support (Lemma E with `j = l`). These cases are the n = 5 full-support
  problem at the same parameters, with `G[sigma]` domination-free by (a).
- (c) The genuinely new Stage 4 cases are those with `G[sigma]` domination-free and `l` not dominated in `G`.

This says which cases need computation; I have not measured how many there are.

---

## 7. Remark F (ray monotonicity of stability)

For real `beta, t`, `K(beta, t) = I + beta N + t J`. If `K(beta, t)` is stable (all eigenvalues have `Re > 0`), then so
is `K(lambda beta, lambda t)` for every `lambda in (0, 1]`.
*Proof.* `K(lambda beta, lambda t) = (1 - lambda) I + lambda K(beta, t)`, whose eigenvalues `(1 - lambda) + lambda kappa`
have positive real part. ∎

So the stability region in the `(beta, t)` plane is star-shaped about `(0, 0)`. Along the ray `t = r beta`,
`K = I + beta (N + r J)` is stable iff `beta < -1 / min Re spec(N + r J)` (all `beta > 0` if that minimum is `>= 0`).
Lemma C' is the `r -> inf` end of this. The remark concerns (3) only: (1) and (2) are not monotone along rays in general.
It is recorded because it is exact and cheap; nothing above uses it.

---

## 8. Review of fastregion.py

I reviewed `fastregion.py` as committed in 0e01ccd (sha256 prefix a011e23c5fa14138; unchanged in the working tree at
review time), with sympy 1.14.0 and Python ground types (gmpy2 is not installed). I read the code, not just the docstring.

**R1. The docstring mathematics is now proved.**

- "I - W_sigma = eps K", the rank-one identities, the off-support formula, and the conditions list: A(i)-(v) and Corollary A.1.
- "Boundary set: roots of a_i, D, o_k and Delta_(deg-1)": Lemma B(d), and Lemma B-inf for `t = INF`.
- The `t = INF` "one-parameter problem": Lemma C, with the uniformity caveat in §3 and in R9.
- The `assert omega[0] == [] and omega[1] == [m]` in `beta_data` is A(iv).

**R2. Stage 2 is complete.**

- `polys2 = a_i, o_k, D, Delta` is exactly `P_t` (finite `t`; `Delta = Delta_(m-1)` of `den * p`, whose roots match those
  of `p` because `Delta_k(c p) = c^k Delta_k(p)` with `c = den > 0`), or `P_inf` (`t = INF`; `Delta = Delta_(m-2)(omega)`,
  `D = s`).
- The early "identically zero -> never" branches are correct by B(d).
- `decide_point` implements Corollary A.1 for finite `t`, and the definition of SFP_inf (with `s` in place of `D`) for `t = INF`.
- The interpolation degree bounds are valid: `d <= m`; `a_i <= m-1`; the mu-coefficients of `det(mu I + C)` and of
  `det(mu I + C + J)` are sums of k x k principal minors, so their degree is `<= k <= m`; a Hurwitz minor has degree
  `<=` the sum of its row-maximum degrees. So every interpolation is exact, and the extra check point is redundant but harmless.
- `idet` (Bareiss with row swaps) and `icharpoly` (Faddeev-LeVerrier, exact division asserted) are exact.

**R3. Stage 1 is valid.**

- Claim: if SFP holds at some legal beta, then at that beta every `a_i` has the sign of `D` (nonzero) and every `o_k` has
  the opposite sign. So the t-free predicate `S(beta) := "all a_i nonzero with one sign sgn, and sgn*o_k < 0 for all k"`
  holds there.
- In particular that beta is not a root of any `a_i` or `o_k`. So it lies in an open cell of `L_t` minus those roots, on
  which `S` is constant (IVT), and `S` holds at that cell's sample.
- Hence "S false at every Stage 1 sample" implies "SFP false on all of `L_t`". This is valid; the sign pattern of `a` is
  indeed constant between consecutive roots of the `a_i`.
- **Sound sharpening available:** by Lemma A.2 (and C'(h) for INF), SFP forces `D > 0` (resp. `s > 0`). So Stage 1 can
  require `all a_i > 0 and all o_k < 0` instead of "common sign". This removes the all-negative cells. Not applied; it is
  a suggestion.

**R4. The first cell and beta_min: no bug, but it depends on sympy internals.**

The relevant internals of sympy 1.14.0 are:

- `_discard_if_outside_interval` keeps an isolating interval `(u, v)` only if `u >= inf`, discards it if `v < inf`, and
  otherwise refines it.
- Refinement is nested.
- A non-degenerate VAS interval never has a root at an endpoint. Every split point is tested exactly, and found roots are
  divided out.
- A rational root is eventually hit exactly, as a degenerate interval `(r, r)`. This is observed, and consistent
  with the finite continued-fraction expansion of a rational; I did not prove it from the full source. If it ever
  failed, the isolation call would loop rather than return a wrong interval.

Consequences, each tested in the scratch file `check_firstcell.py`:

- (i) A returned interval **never straddles** beta_min (`u < beta_min < v` cannot occur).
- (ii) A root exactly at beta_min comes back as `(beta_min, beta_min)` and is filtered by the `iv[0][0] == iv[0][1] == b_min`
  line. This is correct, because beta_min is not legal.
- (iii) The only edge case is `u = beta_min < root < v`. `_cells` refines it until `u > beta_min`; this terminates
  because the root is strictly greater than beta_min.

The previous version's condition `a <= b_min < b` was equivalent to the current one, given (i). Also:

- `strict=True` makes the intervals strictly disjoint (`_disjoint_p` with strict inequalities), so `_cells`' assertion
  `prev < a_` holds.
- Comparisons between sympy's `PythonMPQ` and `fractions.Fraction` are exact: `Fraction` is in
  `PythonMPQ._compatible_types`, and comparison is by integer cross-multiplication.

Risk: all of this rests on undocumented internals. If a future sympy returned straddling intervals for roots below
beta_min, the `while a_ <= b_min` loop would hang; that failure is loud, not silent. If gmpy2 is ever installed, the
ground type becomes `mpq` and the Fraction comparisons should be re-tested.

Recommendation: add an exact replay check of each critical interval. For a non-degenerate interval, check
`f(u) f(v) < 0` for the irreducible (hence square-free) `f`, and `u >= beta_min`. Also check that every irreducible factor
of every member of `P_t` with a root above beta_min appears in the critical list. That is the "certificate round-trips"
standard of CLAUDE.md §3.

**R5. Reporting gap: escalation hides decided stable cells.** In `region()`, if any cell escalates then
`verdict = 'escalate'` and `stable_intervals` is not written, even when other cells are decided stable. The cells remain
in `cert["cells"]`, but a summary that counts `verdict == 'stable_somewhere'` will miss such a case until the escalation
is resolved. Given CLAUDE.md §5 ("the counterexamples are in the undecided pile"), record decided stable cells even when
escalating, or make every downstream count read `cells`.

**R6. Efficiency: unnecessary escalations.** A cell is escalated when `stable is None` at all three samples, even if
`positive` or `off_ok` is already False there. In that case the cell is "not SFP" whatever the stability (Corollary A.1
plus B(d)). This costs work but is not a soundness issue.

**R7. Vanishing minor at an interior sample.** By B(e) the cell is provably not stable. Re-sampling and escalating is
conservative and correct. Using B(e) instead is a policy decision; CLAUDE.md §3 is not changed by this note.

**R8. Endpoint semantics.**

- Critical intervals are open (the root is strictly inside) unless they are degenerate `(r, r)` (an exact rational root).
- The first critical interval stored in `crit` is the *unrefined* one, whose left endpoint may equal beta_min; the
  refined one is used only for the cells. Both are valid isolating intervals.
- Stable regions are open intervals between consecutive critical values (or beta_min). The critical values themselves are
  never SFP (B(d)).

The certificate schema should state all three points.

**R9. Meaning of t = INF.**

- A positive verdict is pointwise in beta with non-uniform `T(beta)` (Lemma C). "stable_somewhere at t = INF" means: for
  every beta in the cell there is an `eps0(beta) > 0` below which SFP holds. It gives no explicit eps0.
- A negative verdict ("never at t = INF") only means that SFP fails *for all large t*. It does **not** exclude
  counterexamples at moderate eps. The limit sweep is therefore a filter for where to look, not a certificate of absence
  at finite eps.
- The q-label at t = INF is the limit `beta - 1`. At finite `t`, `q = (beta-1) t/(t+1) < beta - 1`.

**R10. A cheaper exact test is available at t = INF** (Lemma C'/I).

- The stability clause holds exactly on `(2, beta_c)`. It is nonempty iff `omega(.; 2)` is Hurwitz: one exact Hurwitz
  test of a fixed integer polynomial per (graph, support).
- `e < m(m-1)/2` (more than half the ordered pairs are edges) is a necessary pre-filter.
- Part of `Delta_(m-2)(omega)` is spurious (C'(i)): it adds extra cells, which is harmless.
- Positivity and off-support become (I1) and (I2): a positive solution `w` of `gamma w_i + in_i(w) = 1`, and `in_k(w) < 1`
  outside sigma.

**R11. No floating point.** There is none in fastregion.py or in the decision paths it calls. `ctln.hurwitz_stable` is
sound in all three outputs:

- True: by (b0) "<=".
- False: by (b0) "=>" and (b1).
- None: a zero coefficient or a zero minor. Mathematically this is "not Hurwitz at that point" (by (b0) and (b1)); the
  project labels it undecided, by policy.

**R12. `crosscheck`** compares against `ctln.analyze_support`, which uses the same strict (2) and the same three-valued
stability. The two pipelines are consistent, and my Geneson runs at `t = 40, 60, 100` passed the crosscheck.

---

## 9. Not proved here / open

- An explicit eps0 for t = INF cells (non-uniform `T(beta)`). Also unknown: whether, at fixed beta, the stable set in `t`
  is a single interval `(t*, inf)`. It is for Geneson's graph at `beta = 54/25`; in general this is UNPROVED.
- A finite-t analogue of the trace bound C'(f). Geneson's Thm 4.1 is the known finite-t bound, and nothing here improves it.
- The Gantmacher section numbers (§6, §7) and the Orlando page numbers are quoted from memory; check them before
  external citation. The theorems themselves are standard. The one consequence of Orlando's formula that we use, (b3),
  is proved here.
- The Lean formalization of the certificate lemma has not been started. Corollary A.1, Lemma A.2, (b1), (b3), D.1,
  Lemma D and Lemma E are all elementary and are good first targets. B and C need Routh-Hurwitz and Rouché, which may
  not be in Mathlib.
- The counts behind Corollary E.4 (how many Stage 4 cases survive) have not been measured.

---

## Appendix: sanity checks run (tests, not proofs)

Scratch scripts are in `/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/prover/`. That
directory is session-specific: copy the scripts into the repo if they are to be kept.

| script | what it checks | result |
|---|---|---|
| `check_lemmaA.py` | A(i)-(v), C'(a),(b),(c), symbolic in beta, t, mu | 40 random (graph, support) cases, n <= 7, m <= 6: all pass |
| `check_orlando.py` | (b4) symbolic for n = 2..5 and exact for n = 6, 7; the (b3) kernel vector; Geneson `G` spectrum | all pass; `-1/min Re spec G = 2.1663550278491...`, and the sextic vanishes there to 1e-26 |
| `check_geneson_sextic.py` | exact: the pair-sum polynomial of the quartic, with `x = -2/beta`, equals the sextic | ratio exactly 1; the sextic is irreducible |
| `check_domination.py` | Lemma D (incl. eps >= 1), Lemma E (strict and weak), the non-uniform theta example | 300 graphs with a dominated node; 14,328 D-checks and 26,856 E-checks pass; the example reproduces |
| `check_lemmaC.py` | omega-Hurwitz vs p_t-Hurwitz at t = 1e7; `lambda_max - m t -> 1 + beta e/m`; Remark F | 149/149 agree (1 boundary case skipped); `4.96000238` vs `4.96`; 27/27 |
| `check_lemmaI.py` | (*), the det identity, (I.3), (I.4), (I.2) | 400 cases: 111 one-signed matched `w > 0` exactly; 230 off-support and 400 stability comparisons agree |
| `check_firstcell.py` | sympy isolation at beta_min: exact root at beta_min, roots 1e-9 above and below, irrational near-misses | behaves as stated in R4 |
| `check_region_geneson.py` | `fastregion.region` on Geneson's graph at t = INF, 40, 60, 100 (with crosscheck) | INF: stable `(2, beta_c)` with the sextic root; finite t: one stable cell starting at beta_min |

The exact Geneson computation of §4 is re-runnable from here, standalone with sympy only:

```python
import sympy as sp
z, x, y, b = sp.symbols('z x y b')
N = sp.Matrix([[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]])  # Geneson (19)
m = 6
V = sp.Matrix(m, m-1, lambda i, j: 1 if i == j else (-1 if i == m-1 else 0))      # basis of 1-perp
G = (V.T*V).inv() * V.T * N * V                                                  # similar to U^T N U  (C'(g))
print(sp.factor(G.charpoly(z).as_expr()))                                        # z*(6z^4+11z^3+16z^2+10z+3)/6
g = 6*y**4 + 11*y**3 + 16*y**2 + 10*y + 3
R = sp.Poly(sp.resultant(g, g.subs(y, x - y), y), x)                             # roots nu_i+nu_j, all ordered pairs
Q, r = sp.div(R, sp.Poly(16*g.subs(y, x/2), x)); assert r.is_zero                 # drop i = j
c, fl = sp.factor_list(Q.as_expr()); assert len(fl) == 1 and fl[0][1] == 2        # Q = c * S^2
S = fl[0][0]
sext = 797*b**6 - 6468*b**5 + 22544*b**4 - 44440*b**3 + 53280*b**2 - 38016*b + 13824
print(sp.cancel(sp.expand(S.subs(x, -2/b)*b**6) / sext))                         # 1
```
