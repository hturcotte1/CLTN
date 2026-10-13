# Lemma I: the eps -> 0 limit as a spectral condition on the compressed adjacency

Prover notes, 2026-10-13. A rigorous, self-contained treatment of the derivation in
`notes/limit_spectral.md` (Claims 1-3 and its "Consequences"), plus graph-theoretic corollaries.
Conventions are exactly CLAUDE.md §2; nothing here changes them. This file builds on
`notes/lemmas.md` (Lemmas A-E, C', I) and `notes/lemma_F_reduction.md` (Lemma H), and is consistent
with both; where a statement here sharpens one there, the relation is spelled out (§9).

Each item is labelled:

- **PROVED here** -- the proof is in this file and uses only the tools in §1 (all proved here) plus
  linear algebra (determinants, multilinearity, `det(XY) = det X det Y`, the spectral theorem for real
  symmetric matrices, unique factorisation in `C[lambda]`, Bolzano-Weierstrass).
- **CITED** -- a classical theorem with a reference; used only where stated.
- **UNPROVED / OPEN** -- labelled as such.

sympy/numpy were used only for sanity checks (Appendix). No statement below rests on a test.

| item | status | hypotheses |
|---|---|---|
| T1, T2 (adjugate rank-one update), R (continuity of roots) | PROVED | none |
| I.1 compression identity | PROVED | any real (indeed complex) square X, unit u |
| I.2 limit eigenvalues (fixed C, and C_k -> C) | PROVED | m >= 2, C real |
| I.3 limit stability; uniform version on compact beta-sets; finite-t form | PROVED | beta > 1 (beta > 2 for legality) |
| I.4 positivity: exact for every t > 0, hence also in the limit | PROVED | beta > 1, t > 0 |
| I.5 density bound | PROVED | m >= 2 |
| I.6 Bendixson / symmetric part, both directions, graph corollaries | PROVED (Bendixson proved) | m >= 2 |
| I.7 exact eigenvalues from eigenvectors of A orthogonal to 1 (twins, modules) | PROVED | beta > 1, every t |
| uniform statement for beta -> infinity jointly with t | OPEN (§4, remark I.3.5) | - |

---

## 0. Setting

Graph `G` on `[n]`, `A[i][j] = 1` iff `i -> j`. Support `sigma = (sigma_1 < ... < sigma_m)`, `m >= 1`.
Throughout, **`A` denotes the induced adjacency `A_sigma`**: the `m x m` 0/1 matrix with
`A_ab = 1` iff `sigma_a -> sigma_b`; its diagonal is 0. (For the full support `sigma = [n]` this is the
adjacency matrix itself.) `|E|` is the number of directed edges inside sigma, `|E| = 1^T A 1`, and
`e := m(m-1) - |E|` is the number of ordered non-adjacent pairs.

- `J = 1 1^T`; `N := J - I - A^T` (receiver-first nonedge matrix: `N_ab = 1` iff `a != b` and
  `sigma_b -/-> sigma_a`); so `e = 1^T N 1`.
- `t := (1-eps)/eps`, `beta := 1 + delta/eps`, `gamma := 1 - 1/beta`.
- `C := I + beta N`, `K := C + t J`.
- `chi_X(lambda) := det(lambda I - X)` (monic, degree = size of X).
- "`X` is **positive stable**" means every eigenvalue of `X` has real part `> 0`.

**(S1)** `I - W_sigma = eps K`. *Proof.* `eps beta = eps + delta`, `eps t = 1 - eps`. Entry `(a,b)` of `eps K` is
`eps[a=b] + (eps+delta) N_ab + (1 - eps)`: for `a = b` it is `1`; for `a != b` with `sigma_b -> sigma_a`
(`N_ab = 0`) it is `1 - eps = 1 - W[sigma_a][sigma_b]`; for `a != b` with `sigma_b -/-> sigma_a` (`N_ab = 1`)
it is `1 + delta = 1 - W[sigma_a][sigma_b]`. (This is Lemma A(i) of lemmas.md.) ∎

**(S2)** For `eps > 0`: condition (3) holds iff `K` is positive stable. *Proof.* `-I + W_sigma = -eps K` by (S1),
so its eigenvalues are `-eps kappa`, `kappa in spec K`, and `Re(-eps kappa) < 0` iff `Re kappa > 0`. ∎

**(S3)** The map `(eps, delta) -> (t, beta)` is a bijection `(0,1) x (0,inf) -> (0,inf) x (1,inf)`, inverse
`eps = 1/(1+t)`, `delta = (beta-1)/(1+t)`. `q := delta(1-eps)/eps = (beta-1) t/(t+1)`, and
legal (`q > 1`) iff `t (beta - 2) > 1` iff `beta > 2 + 1/t`. In particular every legal point has `beta > 2`,
hence `gamma in (1/2, 1)`. **"The eps -> 0 limit at fixed beta"** means `t -> +inf` with `beta` fixed,
i.e. the ray `delta = (beta - 1) eps`, `eps -> 0+`; along it `q -> beta - 1` and `gamma = q_lim/(q_lim + 1)`
with `q_lim := beta - 1`.

**(S4) The identity behind everything.** For `beta != 0`,

    C = (1 - beta) I + beta J - beta A^T = -beta (A^T + gamma I) + beta J,
    K = (1 - beta) I - beta (A^T - tau J),       tau := 1 + t/beta.

*Proof.* Expand `I + beta(J - I - A^T)`; `-beta gamma = 1 - beta`; `beta J + t J = beta tau J`. ∎

**Compression.** For `m >= 2` fix `U in R^(m x (m-1))` with orthonormal columns spanning `1^perp`
(so `U^T U = I_(m-1)`, `U U^T = P := I - J/m`, `U^T 1 = 0`). For any `m x m` matrix `X` put
`X|perp := U^T X U`. If `U'` is another such matrix then `U' = U R` with `R := U^T U'` orthogonal
(`U' = P U' = U U^T U'` and `R^T R = U'^T P U' = I`), so `U'^T X U' = R^T (X|perp) R`. Hence the
characteristic polynomial, the spectrum, the spectral abscissa and the spectrum of the symmetric part of
`X|perp` do not depend on the choice of `U`, and the notation is unambiguous for all uses below.
Also: `J|perp = U^T 1 1^T U = 0`, `I|perp = I_(m-1)`, and `(X^T)|perp = (X|perp)^T`.

`alpha(M) := max Re spec(M)` (spectral abscissa); `sym(M) := (M + M^T)/2`.

---

## 1. Tools (all PROVED)

Valid over any commutative ring `R` unless stated (we use `R = R, C, Z[beta, t, lambda]`).

**(T1) Cramer, column form.** For `M in R^(m x m)` and `b in R^m`: `(adj(M) b)_i = det M[i<-b]`, where
`M[i<-b]` is `M` with column `i` replaced by `b`.
*Proof.* `adj(M)_ij = (-1)^(i+j) det M(j|i)` (`M(j|i)` deletes row `j` and column `i`), so
`(adj(M) b)_i = sum_j (-1)^(i+j) b_j det M(j|i)`, the Laplace expansion of `det M[i<-b]` along column `i`. ∎

**(T2) Adjugate rank-one update.** For `M in R^(m x m)`, `u in R^m`, `tau in R`:

    det(M + tau u u^T) = det M + tau u^T adj(M) u.

No invertibility is assumed. If `M` is invertible this reads `det M (1 + tau u^T M^-1 u)`.
*Proof.* Column `j` of `M + tau u u^T` is `M_*j + tau u_j u`. Expand `det` multilinearly in the columns: a term
in which two or more columns are replaced by multiples of `u` has two proportional columns and vanishes. The
terms with exactly one replaced column `j` sum to `tau sum_j u_j det M[j<-u] = tau sum_j u_j (adj(M) u)_j`
by (T1). If `M` is invertible, `adj M = det M M^-1`. ∎
(For `u = 1` this is (T2) of lemmas.md.)

**(R) Continuity of the roots of monic polynomials.** Let `p_k, p` be monic polynomials of degree `n >= 1` over
`C` with the coefficients of `p_k` converging to those of `p`. Write `p = prod_(j=1..n) (lambda - z_j)` (roots
listed with multiplicity). Then there are listings `w^(k)_1, ..., w^(k)_n` of the roots of `p_k` (with
multiplicity) such that `max_j |w^(k)_j - z_j| -> 0`.
Consequently: if `z` is a root of `p` of multiplicity `r` and `0 < rho` is less than the distance from `z` to
every other root of `p`, then for all large `k` the open disk `B(z, rho)` contains exactly `r` roots of `p_k`
(with multiplicity).

*Proof.*
1. *Cauchy bound.* A root `lambda` of a monic `lambda^n + c_1 lambda^(n-1) + ... + c_n` satisfies
   `|lambda| <= 1 + M`, `M := max |c_i|`. Indeed if `|lambda| > 1 + M` then
   `|c_1 lambda^(n-1) + ... + c_n| <= M (|lambda|^n - 1)/(|lambda| - 1) < |lambda|^n`.
   The coefficients of `p_k` converge, so they are bounded, and all roots of all `p_k` lie in one disk `|lambda| <= B`.
2. Let `d_k := min over listings w of the roots of p_k of max_j |w_j - z_j|` (a minimum over finitely many
   permutations of one fixed listing). Suppose `d_k` does not tend to 0: then `d_k >= eta > 0` along a subsequence.
   Fix any listing `v^(k)` of the roots of `p_k`; it lies in the compact set `{|lambda| <= B}^n`, so a further
   subsequence has `v^(k) -> v in C^n`.
3. The coefficients of `p_k` are, up to sign, the elementary symmetric polynomials of `v^(k)`; these are continuous,
   so along the subsequence they converge to those of `v`. They also converge to the coefficients of `p`. Hence
   `prod_j (lambda - v_j) = p = prod_j (lambda - z_j)`, and by unique factorisation in `C[lambda]` there is a
   permutation `pi` with `v_(pi(j)) = z_j` for all `j`.
4. Then `max_j |v^(k)_(pi(j)) - z_j| -> 0` along the subsequence, so `d_k -> 0` there, contradicting `d_k >= eta`.
5. Consequence. Let `z` have multiplicity `r`, let `D_z` be its distance to the other roots of `p` (`+inf` if
   there are none) and `0 < rho < D_z`. For all large `k`, `d_k < min(rho, D_z - rho)`; take an optimal listing.
   The `r` roots of `p_k` paired with the copies of `z` lie within `d_k < rho` of `z`. Every other root of `p_k` is
   paired with some `z_j != z` and lies within `d_k` of it, so its distance to `z` is
   `>= |z_j - z| - d_k >= D_z - d_k > rho`. So `B(z, rho)` contains exactly `r` roots of `p_k`. ∎

The sequential form implies the usual `for every rho there is T` form for one-parameter families
`p_t`, `t -> inf` (argue by contradiction along a sequence `t_k -> inf`); this is how it is used below.

---

## 2. Lemma I.1 (compression identity) -- PROVED

**Lemma I.1.** Let `m >= 2`, let `X` be any complex `m x m` matrix, let `u in R^m` be a unit vector, and let
`U in R^(m x (m-1))` have orthonormal columns spanning `u^perp`. Then

    det(U^T X U) = u^T adj(X) u.

Consequently, for every `lambda in C`,

    det(lambda I_(m-1) - U^T X U) = u^T adj(lambda I - X) u,
    and for u = 1/sqrt(m):   chi_(X|perp)(lambda) = (1/m) 1^T adj(lambda I - X) 1.          (I.1a)

*Proof.*
1. `Q := [U | u]` is an `m x m` real matrix with orthonormal columns, hence orthogonal, `det(Q)^2 = 1`, and
   `Q^T u = e_m`.
2. For every `tau in C`: `Q^T (X + tau u u^T) Q = Q^T X Q + tau e_m e_m^T`.
3. Left side: `det(Q^T (X + tau u u^T) Q) = det(Q)^2 det(X + tau u u^T) = det X + tau u^T adj(X) u` by (T2).
4. Right side: `Q^T X Q + tau e_m e_m^T` differs from `Q^T X Q` only in entry `(m,m)`. Expanding along row `m`,
   the determinant is affine in that entry with slope equal to its cofactor, which is the determinant of the
   leading `(m-1) x (m-1)` block of `Q^T X Q`, namely `det(U^T X U)`. So the right side is
   `det(Q^T X Q) + tau det(U^T X U) = det X + tau det(U^T X U)`.
5. Both sides are affine in `tau` and agree for all `tau`; compare the coefficients of `tau`. ∎

For (I.1a) apply the lemma to `lambda I - X`, using `U^T (lambda I - X) U = lambda I_(m-1) - U^T X U`. ∎

*Alternative proof (remark).* `adj(Q^T X Q) = adj(Q) adj(X) adj(Q^T) = det(Q)^2 Q^T adj(X) Q`; its `(m,m)` entry is
the `(m,m)` cofactor `det(U^T X U)` on one side and `u^T adj(X) u` on the other. (This needs `adj(YZ) = adj(Z) adj(Y)`,
which the proof above avoids.)

**Corollary I.1.1 (characteristic polynomial of a rank-one shift) -- PROVED.** For every `m x m` complex `X` and
`tau in C`:

    chi_(X + tau J)(lambda) = chi_X(lambda) - tau m chi_(X|perp)(lambda).

In particular `m chi_(X|perp) = chi_X - chi_(X+J)`, so for an integer matrix `X` the polynomial `m chi_(X|perp)`
has integer coefficients and is computable without square roots or a choice of basis.
*Proof.* (T2) with `M = lambda I - X`, `u = 1`, and `-tau` in place of `tau`; then (I.1a). ∎

*Dictionary with lemmas.md.* There `omega(mu) := 1^T adj(mu I + C) 1`. Since `adj(-Y) = (-1)^(m-1) adj(Y)`,
`omega(mu) = (-1)^(m-1) m chi_(C|perp)(-mu) = m det(mu I + C|perp)`, which is C'(a) there. The roots of `omega` are the
negatives of the eigenvalues of `C|perp`.

---

## 3. Lemma I.2 (limit eigenvalues) -- PROVED

The statement is proved for a convergent *sequence* of matrices `C_k -> C_inf` with `t_k -> +inf`. The fixed-`C`
statement is the constant sequence; the sequence form is what gives uniformity on compact `beta`-sets in I.3(f).

**Lemma I.2.** Let `m >= 2`. Let `C_k` (`k = 1, 2, ...`) be real `m x m` matrices with `C_k -> C_inf` entrywise, and
`t_k` real with `t_k -> +inf`. Put `K_k := C_k + t_k J` and `kappa_inf := (1/m) 1^T C_inf 1`. Let `z_1, ..., z_r` be the
distinct eigenvalues of `C_inf|perp`, with algebraic multiplicities `m_1, ..., m_r` (`sum m_l = m - 1`), and let
`rho_0` be half the least distance between two distinct `z_l` (`rho_0 := +inf` if `r = 1`).

- **(a) Secular equation.** For every real `t` and every `C`:

      det(lambda I - (C + t J)) = det(lambda I - C) - t 1^T adj(lambda I - C) 1 = chi_C(lambda) - t m chi_(C|perp)(lambda).

- **(b) The unbounded eigenvalue.** For all large `k`, `K_k` has exactly one eigenvalue `lambda*_k` (counted with
  multiplicity) in the region `|lambda/t_k - m| < m/2`. It is real and simple, `lambda*_k / t_k -> m`, and

      lambda*_k - m t_k -> kappa_inf.

- **(c) The bounded eigenvalues.** For every `rho in (0, rho_0)` and all large `k`, each open disk `B(z_l, rho)`
  contains exactly `m_l` eigenvalues of `K_k` (with multiplicity). These `m - 1` eigenvalues and `lambda*_k` are all the
  eigenvalues of `K_k`.

**Fixed-`C` form.** For a fixed real `C` and every `rho in (0, rho_0)` there is `T` such that for every `t > T`: `K = C + tJ`
has a real simple eigenvalue `lambda*(t)` with `|lambda*(t) - m t - kappa| < rho` (`kappa := (1/m) 1^T C 1`), and each
`B(z_l, rho)` contains exactly `m_l` further eigenvalues. So as `t -> +inf` exactly one eigenvalue is unbounded (real,
`~ m t`) and the other `m - 1` converge, with multiplicity, to the eigenvalues of `C|perp`, i.e. to the roots of
`1^T adj(lambda I - C) 1`. *Proof:* if this failed for some `rho`, there would be `t_k -> inf` violating it, contradicting
the sequence form with `C_k = C`. ∎

*Proof of the sequence form.*

0. **(a)** is Corollary I.1.1 with `X = C`, `tau = t`.
   Write `P_k := chi_(K_k)`, `chi_k := chi_(C_k) = lambda^m + sum_(j<m) c_(k,j) lambda^j`, and
   `g_k := chi_(C_k|perp) = lambda^(m-1) + sum_(j<m-1) g_(k,j) lambda^j`; put `g_(k,m-1) := 1`. The coefficients of
   `chi_k` are polynomials in the entries of `C_k`, and by (I.1a) so are those of `g_k = (1/m) 1^T adj(lambda I - C_k) 1`.
   Hence they converge to the coefficients `c_(inf,j)`, `g_(inf,j)` of `chi_(C_inf)` and `g_inf := chi_(C_inf|perp)`,
   and are bounded in `k`. By (a), `P_k = chi_k - t_k m g_k`, a monic polynomial of degree `m`.
1. **Locating the unbounded root (rescale by `t_k`).** `Q_k(mu) := t_k^(-m) P_k(t_k mu)` is monic of degree `m`:

       Q_k(mu) = mu^m - m mu^(m-1) + sum_(j<m) c_(k,j) t_k^(j-m) mu^j - m sum_(j<m-1) g_(k,j) t_k^(j+1-m) mu^j.

   Every term in the two sums is a bounded quantity times a negative power of `t_k` (exponents `j - m <= -1` and
   `j + 1 - m <= -1`), so the coefficients of `Q_k` converge to those of `Q_inf(mu) := mu^(m-1)(mu - m)`, whose roots are `m` (simple) and `0` (multiplicity
   `m - 1`), at distance `m`. By (R) with `rho = m/2`: for large `k`, `Q_k` has exactly one root `mu_k` in `B(m, m/2)`
   and exactly `m - 1` roots in `B(0, m/2)`. `Q_k` has real coefficients, so `conj(mu_k)` is also a root in
   `B(m, m/2)`; by uniqueness `mu_k` is real, and it is simple. Applying (R) with any `rho' < m/2`, the unique root of
   `Q_k` in `B(m, rho') ⊂ B(m, m/2)` is `mu_k` for large `k`, so `mu_k -> m`. The roots of `P_k` are `t_k` times those
   of `Q_k`; put `lambda*_k := t_k mu_k`. Then `lambda*_k` is a real simple root of `P_k`, the only one with
   `|lambda/t_k - m| < m/2`, and `lambda*_k / t_k -> m`; in particular `lambda*_k -> +inf`.
2. **Deflation.** For large `k`, `P_k(lambda) = (lambda - lambda*_k) R_k(lambda)` with `R_k` real and monic of degree
   `m - 1`, `R_k = sum_(j=0..m-1) r_(k,j) lambda^j`, `r_(k,m-1) = 1`. Comparing coefficients of `lambda^j`,
   `0 <= j <= m-1` (with `r_(k,-1) := 0`):

       r_(k,j-1) - lambda*_k r_(k,j) = c_(k,j) - t_k m g_(k,j),
       so  r_(k,j) = (m t_k / lambda*_k) g_(k,j) - c_(k,j) / lambda*_k + r_(k,j-1) / lambda*_k.

   Here `m t_k / lambda*_k -> 1`, `c_(k,j)/lambda*_k -> 0` (bounded over unbounded), and `g_(k,j) -> g_(inf,j)`.
   Induction on `j = 0, 1, ..., m-1`: if `r_(k,j-1)` converges (true for `j = 0`), it is bounded, so
   `r_(k,j-1)/lambda*_k -> 0` and `r_(k,j) -> g_(inf,j)`. Hence `R_k -> g_inf` coefficientwise.
3. **The bounded roots.** `R_k` and `g_inf` are monic of degree `m - 1`, so (R) applies: for `rho in (0, rho_0)` and large
   `k`, each `B(z_l, rho)` contains exactly `m_l` roots of `R_k`, which accounts for all `m - 1` of them. The roots of
   `P_k` are `lambda*_k` and the roots of `R_k`. Since `lambda*_k -> +inf`, for large `k` it lies outside every
   `B(z_l, rho)` and is not a root of `R_k`. This is (c), and it re-confirms that `lambda*_k` is simple.
4. **The `O(1)` term.** The roots of `P_k` sum to `tr K_k = tr C_k + m t_k`, and those of `R_k` to `-r_(k,m-2)`, which
   tends to `-g_(inf,m-2) = tr(C_inf|perp)`. Hence `lambda*_k - m t_k = tr C_k + r_(k,m-2) -> tr C_inf - tr(C_inf|perp)`.
   Finally `tr C - tr(U^T C U) = tr(C (I - U U^T)) = tr(C J)/m = (1/m) 1^T C 1`. ∎

**Specialisation.** For `C = I + beta N`: `kappa = 1 + beta e/m` (`e = 1^T N 1`), so the unbounded eigenvalue is
`m t + 1 + beta e/m + o(1)`, as in lemmas.md C'(c). The eigenvalues of `K` are the negatives of the roots of
`p_t(mu) = det(mu I + K)`; so (b), (c) are Lemma C(a) of lemmas.md, there proved with Rouché's theorem. Here the only
analytic input is (R), proved in §1.

---

## 4. Lemma I.3 (limit stability) -- PROVED

Throughout `m >= 2` and `alpha := alpha(A|perp) = max Re spec(A|perp)`, a number that depends only on the graph
`G[sigma]` (not on any parameter). By I.5 below, `alpha >= -1`, with equality iff `G[sigma]` is complete bidirected.
(I.5 uses only a trace computation, so there is no circularity.)

**Lemma I.3.** Let `beta > 0` and `gamma = 1 - 1/beta`.

- **(a) Limit matrix.** `C|perp = (1 - beta) I_(m-1) - beta (A|perp)^T`, and

      spec(C|perp) = { 1 - beta - beta mu : mu in spec(A|perp) }     (with algebraic multiplicities).

- **(b) Limit stability.** Every eigenvalue of `C|perp` has real part `> 0` iff `alpha < -gamma`.
  Some eigenvalue of `C|perp` has real part `< 0` iff `alpha > -gamma`. When `alpha = -gamma` some eigenvalue of
  `C|perp` lies on `iR` and none has negative real part.
- **(c) Along the ray `delta = (beta - 1) eps`, `eps -> 0+` (fixed `beta > 1`).**
  - (c1) If `alpha < -gamma`, there is `T` such that `K(beta, t)` is positive stable for all `t > T`; i.e. condition (3)
    holds at `(eps, delta) = (1/(1+t), (beta-1)/(1+t))` for all `t > T`. If `beta > 2` these points are legal once
    also `t > 1/(beta - 2)`.
  - (c2) If `alpha > -gamma`, there is `T` such that for all `t > T`, `K(beta, t)` has an eigenvalue with negative real
    part; (3) fails.
  - (c3) If `alpha = -gamma`, Lemma I.2 decides nothing (a limit eigenvalue lies on `iR`). This is a boundary case in the
    sense of CLAUDE.md §3 and must be treated as undecided by any procedure based on the limit.

  So, for fixed `beta`: "(3) holds for all sufficiently small `eps` on the ray" is implied by `alpha < -gamma` and
  implies `alpha <= -gamma`; the same `<=` follows already from "(3) holds along some sequence `eps_k -> 0` on the
  ray" (by (c2)). The strict `<` cannot be claimed in the converse direction at `beta = beta_c` without further work.
- **(d) Existence over legal rays.** There is a `beta > 2` such that (3) holds at all sufficiently small `eps` on the
  ray `delta = (beta-1) eps` (these points being legal) **iff** `alpha < -1/2`. The same equivalence holds with "along
  some sequence `eps_k -> 0` on the ray" in place of "at all sufficiently small `eps`".
- **(e) The limit-stable `beta`-set.** Let `L := { beta > 2 : alpha < -gamma(beta) }`. Then

      L = (2, beta_c),   beta_c := 1/(1 + alpha)  if alpha > -1,   beta_c := +inf  if alpha = -1,

  nonempty iff `alpha < -1/2`. For `beta > beta_c` (when `beta_c < inf`), (3) fails for all large `t` by (c2).
  At `beta = beta_c` the limit decides nothing (c3). In terms of the limiting `q = beta - 1`:
  `q_c = beta_c - 1 = -alpha/(1 + alpha)`. `alpha = -1` (so `beta_c = inf`) happens only for the complete
  bidirected graph (I.5).
- **(f) Uniformity on compact `beta`-sets.** Let `[b1, b2] ⊂ (1, inf)`.
  - (f1) If `alpha < -gamma(b2)` (equivalently `alpha < -gamma(beta)` for every `beta in [b1, b2]`, since `gamma` is
    increasing), there is `T` such that `K(beta, t)` is positive stable for every `beta in [b1, b2]` and `t > T`.
  - (f2) If `alpha > -gamma(b1)` (equivalently for every `beta in [b1, b2]`), there is `T` such that `K(beta, t)` has an
    eigenvalue with negative real part for every `beta in [b1, b2]` and `t > T`.

  In particular, on every compact subinterval of `(2, beta_c)` the stability threshold `T` is uniform. It is not
  explicit: compactness gives existence only.
- **(g) Exact finite-`t` form (every `t`, not a limit).** With `tau := 1 + t/beta`,

      K = (1 - beta) I - beta (A^T - tau J),     spec K = { 1 - beta - beta nu : nu in spec(A - tau J) },

  so for `beta > 0`: `K` is positive stable iff `alpha(A - tau J) < -gamma`. As `t -> inf` (so `tau -> inf`), `m - 1`
  eigenvalues of `A - tau J` converge with multiplicity to `spec(A|perp)` and one is real, `= -m tau + O(1)`.
  In the coordinates `gamma = delta/(eps + delta)`, `tau = (1 + delta)/(eps + delta)`, legality reads
  `gamma > 1/2` and `(tau - 1)(2 gamma - 1) > (1 - gamma)^2`.

*Proof.*

- (a) By (S4), `C = (1 - beta) I + beta J - beta A^T`. Compress: `I|perp = I`, `J|perp = 0`, `(A^T)|perp = (A|perp)^T`.
  For the spectrum: `chi_(C|perp)(lambda) = det((lambda - 1 + beta) I + beta (A|perp)^T)
  = (-beta)^(m-1) chi_((A|perp)^T)(nu)` with `nu := (1 - beta - lambda)/beta`, and `chi_(M^T) = chi_M`. The affine
  bijection `lambda <-> nu` (`beta != 0`) carries roots to roots with multiplicities. ∎
- (b) `Re(1 - beta - beta mu) > 0` iff `Re mu < (1 - beta)/beta = -gamma` (as `beta > 0`); `< 0` iff `Re mu > -gamma`;
  `= 0` iff `Re mu = -gamma`. Take the maximum over `mu`. ∎
- (c) Lemma I.2, fixed-`C` form, with `C = I + beta N`. (c1): the eigenvalues `z_l` of `C|perp` satisfy `Re z_l > 0` by
  (b); choose `rho < min(rho_0, min_l Re z_l)`. For `t > T` the `m - 1` bounded eigenvalues lie in disks contained in
  `{Re > 0}`, and `lambda*(t) > m t + kappa - rho > 0` for large `t`. Legality: (S3). (c2): some `z_l` has
  `Re z_l < 0`; choose `rho < min(rho_0, -Re z_l)`; for `t > T` the disk `B(z_l, rho) ⊂ {Re < 0}` contains `m_l >= 1`
  eigenvalues. (c3): nothing to prove. ∎
- (d) (<=) If `alpha < -1/2`, pick `beta > 2` so close to 2 that `-gamma(beta) = -1 + 1/beta > alpha` (possible since
  `-gamma(beta) -> -1/2` as `beta -> 2+`); apply (c1); legality from `t > 1/(beta-2)`. (=>) If (3) holds along a sequence
  `eps_k -> 0` on the ray of some `beta > 2`, then (c2) is excluded, so `alpha <= -gamma(beta) < -1/2` because
  `gamma(beta) > 1/2` strictly for `beta > 2`. ∎ (The strictness lost in (c) is recovered here because `gamma > 1/2`
  is strict on the legal range.)
- (e) `alpha < -gamma(beta)` iff `1 + alpha < 1/beta`. If `alpha = -1` this holds for all `beta > 0`. If `alpha > -1` it
  holds iff `beta < 1/(1 + alpha)`. Intersect with `beta > 2`: the result is `(2, beta_c)`, nonempty iff `beta_c > 2`
  iff `alpha < -1/2`. The remaining statements are (c2), (c3) and I.5. ∎
- (f) (f1) Suppose not: there are `beta_k in [b1, b2]`, `t_k -> inf` with `K(beta_k, t_k)` not positive stable. Pass to a
  subsequence with `beta_k -> beta_inf in [b1, b2]`. `C_k := I + beta_k N -> C_inf := I + beta_inf N`, and by (a), (b)
  every eigenvalue of `C_inf|perp` has positive real part (`alpha < -gamma(b2) <= -gamma(beta_inf)`). Lemma I.2
  (sequence form) with `rho < min(rho_0, min Re spec(C_inf|perp))` puts, for large `k`, `m - 1` eigenvalues of
  `K(beta_k, t_k)` in `{Re > 0}` and the last one `lambda*_k -> +inf`: contradiction. (f2) is the same argument with an
  eigenvalue `z` of `C_inf|perp` with `Re z < 0`. ∎
- (g) (S4) gives the matrix identity; `chi_(M^T) = chi_M` and the affine map `nu -> 1 - beta - beta nu` give the spectrum;
  positive stability iff `Re(1 - beta - beta nu) > 0` for all `nu` iff `Re nu < -gamma`. For the limit: by
  Corollary I.1.1, `chi_(A - tau J)(nu) = chi_A(nu) + tau m chi_(A|perp)(nu)`; apply Lemma I.2 to `C = -A` and
  `t = tau` (the eigenvalues of `-A + tau J` are the negatives of those of `A - tau J`, and `(-A)|perp = -(A|perp)`).
  Coordinates: `gamma = 1 - eps/(eps + delta) = delta/(eps + delta)`; `tau = 1 + t/beta = 1 + (1-eps)/(eps+delta)`.
  Legality `t(beta - 2) > 1` with `beta = 1/(1-gamma)`, `t = beta (tau - 1)`, `beta - 2 = (2 gamma - 1)/(1 - gamma)` reads
  `(tau - 1)(2 gamma - 1) > (1 - gamma)^2`; since `t > 0` gives `tau > 1`, this forces `gamma > 1/2`. ∎

**Combined finite-`t` statement (with I.4 and lemmas.md Lemma C).** Full support (`sigma = [n]`, so (2) is vacuous).
Let `beta in (2, beta_c)` with `A^T + gamma I` nonsingular and `w(gamma) := (A^T + gamma I)^-1 1 > 0`. Then there is
`T(beta)` such that the full support is a stable fixed point at `(eps, delta) = (1/(1+t), (beta-1)/(1+t))` for every
`t > max(T(beta), 1/(beta-2))`, and these points are legal. *Proof.* (1) holds at every `t > 0` by I.4 (positivity is
exactly `t`-free); (3) holds for `t > T(beta)` by (c1). ∎ Compared with Lemma C(b) of lemmas.md, the threshold
`t > |d/s|` there is unnecessary (by I.4 / Lemma H(H3), `d` and `s` have the same sign whenever (1) holds, so
`D = d + t s` never vanishes for `t > 0`), and by (f1) the threshold `T` is uniform on compact subsets of the open set
`{ beta in (2, beta_c) : w(gamma(beta)) > 0 }`. The non-uniformity flagged in lemmas.md §3 ("Uniformity") therefore
comes only from `beta` approaching `beta_c` or the boundary of the positivity set, where (3) resp. (1) degenerates.
For a support of size `m < n` add (I2) of lemmas.md (`in_k(w) < 1`), which is also exactly `t`-free (I.4(f)).

**Remark I.3.5 (what is NOT covered: `beta -> inf` jointly with `t`) -- OPEN.** All statements above are at fixed
`beta`, or uniform for `beta` in a compact set. Small `eps` does not force bounded `beta`: on the ray `t = r beta`
(`r > 0` fixed) one has `K = I + beta (N + r J)`, legal for large `beta`, with `eps = 1/(1 + r beta) -> 0` and `q -> inf`.
There (3) is governed by `spec(N + r J)`, equivalently by `alpha(A - tau J)` at the finite `tau = 1 + r` with `gamma -> 1`
(by (g)), **not** by `alpha(A|perp)`. So "(3) at some legal point with `eps` small implies `alpha < -1/2`" is
**not proved** here, and I do not claim it. For symmetric non-complete `A` the conclusion holds at *every* parameter
for a simpler reason (PROVED): such `A` has a pair `i != j` with no edge either way; `A - tau J` is symmetric, so
`alpha(A - tau J) = lambda_max(A - tau J) >= x^T (A - tau J) x / |x|^2` for `x = e_i - e_j`, and this Rayleigh quotient
is `(A_ii + A_jj - A_ij - A_ji)/2 - tau (1^T x)^2/2 = 0 > -gamma`; so (3) fails for all `beta > 1`, `t > 0` by (g).
For the full SFP question the large-`q` end is covered by Geneson's Thm 4.1 (`q >= n - 2 - (n-3) eps/2`; CITED in
CLAUDE.md), but that theorem concerns SFP, not (3) alone.

---

## 5. Lemma I.4 (positivity) -- PROVED, and exact for every `t > 0`, not only in the limit

Let `beta > 1`, `gamma = 1 - 1/beta in (0, 1)`, `m >= 1`, and put

    M_gamma := A^T + gamma I,      (M_gamma v)_i = gamma v_i + in_i(v),   in_i(v) := sum_(j : sigma_j -> sigma_i) v_j

(`in_i` sums over the in-neighbours of `sigma_i` inside sigma). By (S4), `C = -beta M_gamma + beta J` and
`K = -beta M_gamma + (beta + t) J`. Let `a := adj(C) 1`, `d := det C`, `s := 1^T a` (as in lemmas.md), and, when
`M_gamma` is nonsingular, `w := M_gamma^-1 1` and `W1 := 1^T w`. So `w` is the unique solution of

    gamma w_i + sum_(j -> i) w_j = 1     for every i in sigma.

**(I.4.0) Two identities valid for every `t` and `beta` (no invertibility).**
(i) `adj(C + t J) 1 = adj(C) 1` for every scalar `t` (Lemma A(ii) of lemmas.md).
(ii) `a = (-beta)^(m-1) adj(M_gamma) 1 = (-1)^(m-1) adj(beta A^T + (beta - 1) I) 1`.
*Proof.* (i) By (T1), `(adj(K) 1)_i = det K[i<-1]`. In `K[i<-1]` every column `j != i` is `C_*j + t 1`; subtracting `t`
times column `i` (which is `1`) from each leaves the determinant unchanged and produces `C[i<-1]`. (ii) Apply (i) with
`t = -beta`: `C - beta J = -beta M_gamma`, and `adj(c Y) = c^(m-1) adj(Y)` for a scalar `c`. ∎

**Lemma I.4.**

- **(a) Adjugate vs inverse.** If `d != 0`, then `a = d C^-1 1`; so `a` is one-signed (all entries nonzero, one
  common sign) iff `C^-1 1` is.
- **(b) Solutions of `C u = 1`.** Let `u` satisfy `C u = 1`. Then

      M_gamma u = kappa 1,     kappa := 1^T u - 1/beta = (beta 1^T u - 1)/beta.

  - (b1) If `M_gamma` is nonsingular, then `W1 != 1` and

        u = c w,    c = kappa = 1/(beta (W1 - 1)) != 0.

    Conversely, if `M_gamma` is nonsingular then `C u = 1` has exactly one solution when `W1 != 1` (namely
    `w/(beta(W1 - 1))`) and none when `W1 = 1`; and `d = (-beta)^m det(M_gamma) (1 - W1)`, so `d = 0` iff `W1 = 1`.
  - (b2) If `d != 0` (so `u = C^-1 1`), then `det(-beta M_gamma) = d (1 - beta 1^T u)`. Hence the two "singular cases"
    coincide: `M_gamma` is singular iff `1 - beta 1^T u = 0` iff `kappa = 0`. In that case `M_gamma u = 0` with
    `u != 0`, and `u` is **not** one-signed (so neither is `a`).
  - (b3) If `d = 0` and `M_gamma` is singular, nothing is asserted about `u` (it need not exist or be unique); but `a` is
    then not one-signed, by (c) below ((P2) => `M_gamma` nonsingular).
- **(c) Positivity is `t`-free.** Fix `beta > 1`, `theta > 0` and any `t > 0` (i.e. any `eps in (0,1)` with
  `delta = (beta - 1) eps`). The following are equivalent:
  - (P1) `K(beta, t)` is nonsingular and the (unique) solution `x` of `(I - W_sigma) x = theta 1` is `> 0` entrywise;
  - (P2) `a(beta)` is one-signed;
  - (P3) `M_gamma` is nonsingular and `w > 0` entrywise;
  - (P3') `M_gamma` is nonsingular and `w` is one-signed.

  When they hold: `W1 > 1`, `C w = beta (W1 - 1) 1`, `K w = Den 1` with `Den := beta (W1 - 1) + t W1 > 0`, and

      x = theta (1 + t) w / Den,
      a = (-beta)^(m-1) det(M_gamma) w,
      det C = (-beta)^(m-1) det(M_gamma) * beta (W1 - 1),
      det K = (-beta)^(m-1) det(M_gamma) * Den.

  In particular `sign a = sign det C = sign det K = sign((-1)^(m-1) det M_gamma)`, and `d != 0`.
- **(d) The limit.** Under (P3), `x(t) -> x_inf := theta w / W1 > 0` as `t -> inf` (fixed `beta`). Conversely, if
  `x(t)` (defined for large `t` with `det K != 0`) has a limit `x_inf > 0`, then (P2) holds. So: **in the limit the
  fixed point is positive iff `w(gamma) > 0` entrywise** -- and by (c) this is not only a limit statement: positivity
  holds at one `t > 0` iff at every `t > 0` iff `w(gamma) > 0`.
- **(e) Signs: reconciliation with Lemma A.2 (lemmas.md) and Lemma H(H3) (lemma_F_reduction.md).**

      (1) holds and det K > 0   <=>   a > 0 entrywise   <=>   w > 0  and  (-1)^(m-1) det(A^T + gamma I) > 0.

  A.2 says SFP iff `a > 0`, `o_k < 0` and `p` Hurwitz. Its "`a > 0`" is (1) *together with* `det K > 0`, and the second
  part is implied by (3) (a positive stable real matrix has positive determinant). So under (3): `a > 0` iff `a`
  one-signed iff `w > 0`; the statement "positivity iff `w > 0`" of limit_spectral.md Claim 3 is the sign-free (1),
  and there is no conflict. The case `w > 0` with `(-1)^(m-1) det M_gamma < 0` does occur in principle; there
  `a < 0`, `det K < 0` for **every** `t > 0`, so `K` has a negative real eigenvalue and (3) fails at every `t`
  (a certificate of instability on the whole ray).
- **(f) Off-support, also `t`-free (for supports with `m < n`).** Under (P3), for `k` outside sigma,

      y_k(x) = theta beta (in_k(w) - 1) / Den,     in_k(w) := sum_(b : sigma_b -> k) w_b,

  so condition (2) at `k` holds iff `in_k(w) < 1`, at every `t > 0`. (This is (I.4) of lemmas.md, re-derived directly.)

*Proof.*

- (a) `C adj(C) = d I`, so `a = adj(C) 1 = d C^-1 1` when `d != 0`; multiplying by the nonzero scalar `d` preserves
  one-signedness. ∎
- (b) From `C = -beta M_gamma + beta J`: `1 = C u = -beta M_gamma u + beta (1^T u) 1`, i.e. `M_gamma u = kappa 1`.
  - (b1) `u = kappa w`. `u != 0` since `C u = 1`, so `kappa != 0`. Summing, `1^T u = kappa W1`, and the definition of
    `kappa` gives `kappa = kappa W1 - 1/beta`, i.e. `kappa (W1 - 1) = 1/beta`. So `W1 != 1` and `kappa = 1/(beta(W1 - 1))`.
    Conversely `C w = -beta M_gamma w + beta W1 1 = beta (W1 - 1) 1`, so for `W1 != 1` the vector `w/(beta(W1-1))`
    solves `C u = 1`, and any solution equals it by the above; for `W1 = 1` a solution would force `W1 != 1`, so there
    is none. Determinant: (T2), invertible form, with `M = -beta M_gamma`, `u = 1`, `tau = beta`:
    `d = det(-beta M_gamma) (1 + beta 1^T (-beta M_gamma)^-1 1) = (-beta)^m det(M_gamma) (1 - W1)`. ∎
  - (b2) (T2), invertible form, with `M = C`, `tau = -beta`: `det(C - beta J) = d (1 - beta 1^T C^-1 1)`, and
    `C - beta J = -beta M_gamma`. As `d != 0` and `beta != 0`: `M_gamma` singular iff `1 - beta 1^T u = 0` iff
    `kappa = 0`, and then `M_gamma u = 0`. If `u` were one-signed, row `i` of `M_gamma u` would be
    `gamma u_i + in_i(u)`, nonzero with the sign of `u` (as `gamma > 0` and `A >= 0`), not 0. By (a), `a = d u` is not
    one-signed either. ∎
  - (b3) Covered by (c). ∎
- (c)
  - (P3) => (P1) and the formulas. `C w = beta (W1 - 1) 1` (computed in (b1)). Also `C w = w + beta N w >= w > 0`
    entrywise because `N >= 0` and `beta > 0`. Hence `beta (W1 - 1) > 0`, `W1 > 1`. Then
    `K w = C w + t W1 1 = Den 1` with `Den > 0`, so `x := (theta/eps) w / Den = theta (1+t) w / Den > 0` satisfies
    `eps K x = theta 1`, i.e. `(I - W_sigma) x = theta 1` by (S1). `det K`: (T2), invertible form, with
    `M = -beta M_gamma`, `tau = beta + t`:
    `det K = (-beta)^m det(M_gamma)(1 - (beta + t) W1/beta) = (-beta)^(m-1) det(M_gamma) ((beta+t) W1 - beta)`, and
    `(beta + t) W1 - beta = Den != 0`; so `K` is nonsingular and `x` is the unique solution. The formula for `d` is
    (b1) (`1 - W1 = -(W1 - 1)`); for `a` use (I.4.0)(ii) and `adj(M_gamma) 1 = det(M_gamma) w`. ∎
  - (P1) => (P3). `K x = (theta/eps) 1` with `x > 0` and `K` nonsingular. By (S4),
    `-beta M_gamma x + (beta + t)(1^T x) 1 = (theta/eps) 1`, so `M_gamma x = kappa' 1` for a scalar `kappa'`. Row `i` is
    `gamma x_i + in_i(x) > 0`, so `kappa' > 0`. Suppose `M_gamma v = 0` with `v != 0`. Then `K v = (beta + t)(1^T v) 1`.
    If `1^T v = 0` then `K v = 0`, contradicting nonsingularity. Otherwise `K v` is a nonzero multiple of
    `1 = (eps/theta) K x`, and injectivity of `K` gives `v = c x` with `c != 0`; then `0 = M_gamma v = c kappa' 1 != 0`,
    a contradiction. So `M_gamma` is nonsingular and `w = x/kappa' > 0`. ∎
  - (P3) <=> (P3'). If all `w_i < 0` then every row `gamma w_i + in_i(w)` of `M_gamma w` is `< 0`, not `1`. ∎
  - (P1) <=> (P2) (this is Lemma H(H3), reproduced). (<=) Let `a` be one-signed with sign `s_a`. `C a = d 1`, and
    `(C a)_i = a_i + beta sum_j N_ij a_j` is nonzero with sign `s_a`; so `d != 0` with sign `s_a`. Also `s = 1^T a` has
    sign `s_a`. By (T2) (`M = C`, `u = 1`, `tau = t`), `det K = d + t s`, nonzero with sign `s_a` for `t > 0`. By (I.4.0)(i),
    `K^-1 1 = adj(K) 1 / det K = a / det K > 0`, so `x = (theta/eps) K^-1 1 > 0`. (=>) `a = adj(K) 1 = det(K) K^-1 1
    = det(K) (eps/theta) x`, one-signed. ∎
- (d) Under (P3), `x(t) = theta (1+t) w / (beta (W1 - 1) + t W1) -> theta w / W1`. Conversely, for `t` with
  `det K = d + t s != 0`, `x(t) = theta (1+t) a / (d + t s)` (by (I.4.0)(i)). If `s = 0` this is `theta (1+t) a/d`, which
  has a finite limit only if `a = 0`, and then the limit is `0`, not `> 0`. If `s != 0` the limit is `theta a / s`; it is
  `> 0` only if all `a_i` have the sign of `s`, i.e. (P2). ∎
- (e) By (c), under (P1)-(P3): `sign a = sign det K = sign((-1)^(m-1) det M_gamma)`. So `a > 0` iff [(P2) and
  `det K > 0`] iff [(P3) and `(-1)^(m-1) det M_gamma > 0`]; and (P2) iff (1) (with `K` nonsingular) by (c). If `K` is
  positive stable, its real eigenvalues are `> 0` and its non-real ones pair to `|lambda|^2 > 0`, so `det K > 0`. If
  `det K < 0`, the product of the real eigenvalues is negative, so one of them is negative. ∎
- (f) For `k` outside sigma, `W[k][sigma_b] = -(1 - eps) - (eps + delta) N_kb` with `N_kb = 1` iff `sigma_b -/-> k`.
  With `x = theta (1+t) w / Den`, `(1 - eps)(1 + t) = t` and `(eps + delta)(1 + t) = 1 + (beta - 1) = beta`:
  `sum_b W[k][sigma_b] x_b = -theta (t W1 + beta (W1 - in_k(w))) / Den`, so
  `y_k = theta (Den - t W1 - beta W1 + beta in_k(w)) / Den = theta beta (in_k(w) - 1) / Den`. ∎

**Remark (the weaker reading of (1)).** If (1) is read as "`(I - W_sigma) x = theta 1` has *some* solution `x > 0`"
without nonsingularity, the implication to (P3) can fail when `K` is singular (then `x + ker K` are all solutions).
This never matters for SFP: (3) forces `K` nonsingular (lemmas.md §0).

**(g) Remark: at rational `beta`, `M_gamma` is never singular -- PROVED.** The eigenvalues of the integer matrix `A^T`
are algebraic integers, and a rational algebraic integer is an integer. For rational `beta > 1`, `-gamma = -1 + 1/beta`
is a rational number in `(-1, 0)`, hence not an eigenvalue of `A^T`; so `det M_gamma != 0`. Thus at every rational
sample point (all of Stage A) `w` is defined and positivity is exactly "`w > 0`". The singular case (b2) occurs only at
irrational `beta` with `-gamma(beta)` an eigenvalue of `A`, and then (1) fails at every `t > 0` by (c).

**Geneson's graph (sanity, not proof; script `s_I4_geneson.py`).** At `beta = 54/25` (`gamma = 29/54`): `w > 0`,
`W1 = 451746072/272447123 > 1`, `w` is proportional to the paper's printed vector `v` (test_ctln.py), and
`(-1)^(m-1) det M_gamma > 0`, so `a > 0`. Numbers in the Appendix.

---

## 6. Lemma I.5 (density bound) -- PROVED

**Lemma I.5.** Let `m >= 2`, let `mu_1, ..., mu_(m-1)` be the eigenvalues of `A|perp` with multiplicity, and `|E| = 1^T A 1`.

- **(a)** `tr(A|perp) = sum_i mu_i = -|E|/m`.
- **(b)** `alpha(A|perp) >= (1/(m-1)) sum_i Re mu_i = -|E|/(m(m-1)) >= -1`, and `alpha(A|perp) = -1` iff `G[sigma]` is complete
  bidirected.
- **(c)** If `alpha(A|perp) < -gamma(beta)` (limit stability at `beta`, I.3(b)), then `|E| > gamma(beta) m (m-1)`;
  equivalently `beta < m(m-1)/e` when `e > 0` (so `beta_c <= m(m-1)/e`), equivalently `q_lim = beta - 1 < |E|/e`.
  At any legal `beta` (`> 2`): `|E| > m(m-1)/2`, i.e. more than half of the `m(m-1)` ordered pairs are edges.
  For the full support: `n = 6` needs `|E| >= 16` of 30; `n = 7` needs `|E| >= 22` of 42.
- **(d) Refinement with known eigenvalues.** If `mu_1, ..., mu_r` (`r < m - 1`) are known eigenvalues of `A|perp`
  (counted with algebraic multiplicity, e.g. from I.7), then `alpha(A|perp) >= (-|E|/m - sum_(j<=r) Re mu_j)/(m - 1 - r)`.
  E.g. if `-1` has algebraic multiplicity `>= r`, limit stability at a legal `beta` needs `|E| > m(m - 1 + r)/2`.

*Proof.*
- (a) `tr(U^T A U) = tr(A U U^T) = tr(A (I - J/m)) = tr A - (1^T A 1)/m = 0 - |E|/m`. The trace is the sum of the
  eigenvalues with multiplicity. ∎
- (b) The sum `-|E|/m` is real, so `sum Re mu_i = -|E|/m`, and a maximum is at least the mean. `|E| <= m(m-1)` gives
  `>= -1`. If `alpha = -1` then the mean is `<= -1`, so `|E| >= m(m-1)`: complete. Conversely for the complete graph
  `A = J - I`, `A|perp = -I_(m-1)`, `alpha = -1`. ∎
- (c) `-|E|/(m(m-1)) <= alpha < -gamma` gives `|E| > gamma m(m-1)`. With `|E| = m(m-1) - e` and `gamma = 1 - 1/beta`:
  `1 - e/(m(m-1)) > 1 - 1/beta` iff `beta < m(m-1)/e`; and `beta_c = 1/(1 + alpha) <= 1/(1 - |E|/(m(m-1))) = m(m-1)/e`.
  For `beta > 2`, `gamma > 1/2`. `n = 6`: `|E| > 15`; `n = 7`: `|E| > 21`. ∎
- (d) The remaining `m - 1 - r` eigenvalues have real parts summing to `-|E|/m - sum_(j<=r) Re mu_j`; take the mean.
  With `mu_j = -1` (`j <= r`) the bound is `(r - |E|/m)/(m - 1 - r)`, and `< -1/2` iff `|E| > m(m - 1 + r)/2`. ∎

This is C'(f) of lemmas.md in the `A`-picture. Geneson's graph: `|E| = 19 > 15`, mean real part `-19/30 ~ -0.633`,
`alpha ~ -0.538`. It has an eigenvalue `-1` (bidirected twins, I.7), and (d) with `r = 1` gives `alpha >= -13/24 ~ -0.5417`,
close to the true value.

---

## 7. Lemma I.6 (symmetric part / Bendixson) -- PROVED; both directions stated carefully

**(a) Bendixson's inequality (PROVED).** For a real `k x k` matrix `M` and every eigenvalue `lambda`,

    lambda_min(sym M) <= Re lambda <= lambda_max(sym M).

*Proof.* Let `M v = lambda v` with `v in C^k`, `v* v = 1`. Then `lambda = v* M v` and, as `M` is real,
`conj(lambda) = v* M^T v`; so `Re lambda = v* H v` with `H := sym M` real symmetric. Writing `v = p + i q` (`p, q` real,
`|p|^2 + |q|^2 = 1`), `v* H v = p^T H p + q^T H q` (the cross terms cancel because `H = H^T`), which lies in
`[lambda_min(H), lambda_max(H)]` by the Rayleigh bounds for real symmetric matrices. ∎
The same computation shows `{Re z : z in W(M)} = [lambda_min(sym M), lambda_max(sym M)]`, where
`W(M) := { y* M y : y in C^k, |y| = 1 }` is the numerical range; `W(M)` contains `spec M` but for non-normal `M` can be
much larger than its convex hull.

**(b) Application to `A|perp`.** `sym(A|perp) = U^T sym(A) U`, and for `y in R^(m-1)`, `x := U y` ranges over `1^perp`
with `|x| = |y|` and `y^T (A|perp) y = x^T A x`. Hence

    lambda_max(sym(A|perp)) = max { x^T A x / |x|^2 : 0 != x in 1^perp },   lambda_min likewise with min,
    lambda_min(sym(A|perp))  <=  alpha(A|perp)  <=  lambda_max(sym(A|perp)).

(The lower bound: `alpha >= Re mu >= lambda_min` for any eigenvalue `mu`.)

**(c) The two directions, stated correctly.**
- **Sufficient:** `lambda_max(sym(A|perp)) < -gamma`, i.e. `x^T A x < -gamma |x|^2` for **every** nonzero
  `x in 1^perp`, implies `alpha(A|perp) < -gamma` (limit stability at `beta`).
- **Necessary:** `alpha(A|perp) < -gamma` implies `lambda_min(sym(A|perp)) < -gamma`, i.e. `x^T A x < -gamma |x|^2` for
  **some** `x in 1^perp`. Contrapositive: if `x^T A x >= -|x|^2/2` for **all** `x in 1^perp`, then `alpha >= -1/2` and
  no legal `beta` is limit-stable. A single vector with `x^T A x >= -|x|^2/2` implies nothing about `alpha` from below.
- **The necessary direction is redundant given I.5:** `tr sym(A|perp) = tr(A|perp)`, so
  `lambda_min(sym(A|perp)) <= -|E|/(m(m-1)) <= alpha(A|perp)`; the trace test of I.5 is at least as strong.

**(d) Graph corollaries (pairs).** For `i != j` let `x := e_i - e_j in 1^perp`, `|x|^2 = 2`. Exactly:

    x^T A x = A_ii + A_jj - A_ij - A_ji = -(A_ij + A_ji),     Rayleigh quotient  -(A_ij + A_ji)/2  in {0, -1/2, -1}

for {no edge, exactly one edge, edges both ways} between `i` and `j`. Hence:
- If some pair has **no edge in either direction**, `lambda_max(sym(A|perp)) >= 0`; equivalently the numerical range
  `W(A|perp)` contains the point `0` (take `y = U^T x/sqrt 2`, `y* (A|perp) y = x^T A x / 2 = 0`). This kills the
  **sufficient** test at every `beta`, but **not** the necessary one, and it says **nothing** about `alpha`: Geneson's graph
  has the non-adjacent pair {paper nodes 4, 5}, `lambda_max(sym(A|perp)) ~ +0.207`, yet `alpha ~ -0.538 < -1/2`
  (Appendix). This is how a non-adjacent pair is "survivable" (limit_spectral.md): `A|perp` is far from normal.
- If some pair has **exactly one edge**, `lambda_max(sym(A|perp)) >= -1/2 > -gamma` at every legal `beta`.
- **Hence the sufficient test certifies only the complete graph:** there is a legal `beta` with
  `lambda_max(sym(A|perp)) < -gamma(beta)` iff `G[sigma]` is complete bidirected. (`<=`: `A|perp = -I`,
  `lambda_max = -1 < -gamma` since `gamma < 1`. `=>`: a non-complete graph has a pair with `A_ij + A_ji <= 1`, so
  `lambda_max >= -1/2 > -gamma` for `beta > 2`.)

**Corollary I.6.1 (oriented graphs and tournaments).** If `G[sigma]` is oriented (no pair with edges both ways), then
`A + A^T = J - I - S`, where `S` is the 0/1 symmetric matrix of the pairs with no edge. Hence
`sym(A|perp) = -(1/2) I - (1/2) S|perp`.
- For a **tournament** (`S = 0`): `sym(A|perp) = -(1/2) I`, so by (a) **every eigenvalue of `A|perp` has real part exactly
  `-1/2`**, and `alpha(A|perp) = -1/2`. So no tournament is limit-stable at any legal `beta` (I.3(c2): `-1/2 > -gamma`), and
  the trace bound I.5(b) is attained (all real parts equal the mean `-1/2`).
- For any oriented graph, I.5 already gives `alpha >= -|E|/(m(m-1)) >= -1/2` (`|E| <= m(m-1)/2`), so oriented graphs are
  never limit-stable at a legal `beta`. (Consistent with CLAUDE.md §6: oriented graphs are already settled.)

*Proof.* Oriented: `A_ij + A_ji = 1` for an adjacent pair and `0` for a non-adjacent one, i.e. `A + A^T = J - I - S`.
Compress (`J|perp = 0`) and halve. For `S = 0`, `lambda_min = lambda_max = -1/2` in (a). ∎

This explains part of the exploration count in the Appendix: an exact Hurwitz test of `chi_(A|perp)(z - 1/2)` returns
"undecided" (a zero coefficient or minor) for 738 non-complete `n = 5` graphs; by a float classification (exploration
only) 167 of these have `alpha = -1/2` (the 12 tournaments and 155 mixed graphs) and 571 have `alpha > -1/2` (the
undecided verdict there comes from a zero coefficient). For the legal range both are dead in the limit:
`alpha >= -1/2 > -gamma` whenever `beta > 2` (I.3(c2)); mathematically a zero coefficient or minor already means "not
Hurwitz", i.e. `alpha >= -1/2` (lemmas.md (b0), (b1)).

**What I.6 does and does not imply.** With the standard inner product, the symmetric part adds nothing to I.5 for the
limit problem: its necessary condition is implied by the trace bound, and its sufficient condition holds only for the
complete graph. It does **not** imply that a non-adjacent pair (or any local pattern) prevents limit stability. It does
not bound `alpha` from above usefully for any non-complete graph.

**(e) Remark (CITED; not used).** Lyapunov's theorem: `alpha(M) < -gamma` iff there is a symmetric `P > 0` (positive
definite) with `(M + gamma I)^T P + P (M + gamma I) < 0` (negative definite). (Standard; e.g. Horn & Johnson, *Topics in
Matrix Analysis*, Ch. 2 -- chapter from memory.) Bendixson's sufficient test is `P = I`. For a non-complete graph any
certificate of this kind needs `P != I`; an exact rational `P` (definiteness checked by Sylvester minors) would be an
alternative exact certificate to Routh-Hurwitz for a fixed rational `gamma`. Not pursued.

---

## 8. Lemma I.7 (exact `t`-independent eigenvalues; twins and modules) -- PROVED

This is the "useful direction" corollary that I.6 could not give: an exact, finite-parameter obstruction read off the
graph. It upgrades "dead in the limit" to "dead at every legal parameter" for a recognisable class of graphs.

**Lemma I.7.** Let `m >= 1`, `beta != 0`, `t` real, `tau := 1 + t/beta`, and `nu := (1 - beta - lambda)/beta`.

- **(a)** `chi_(K(beta,t))(lambda) = (-beta)^m [ chi_A(nu) + tau m chi_(A|perp)(nu) ]` (for `m = 1` read `chi_(A|perp) := 1`).
- **(b)** Let `h := gcd(chi_A, chi_(A|perp)) = gcd(chi_A, chi_(A+J))` (monic, in `Q[nu]`). For every root `mu` of `h`,
  `1 - beta - beta mu` is an eigenvalue of `K(beta, t)` for **every** `beta != 0` and **every** `t`. Conversely, if for a
  fixed `beta != 0` some `lambda` is an eigenvalue of `K(beta, t)` for two distinct values of `t`, then
  `lambda = 1 - beta - beta mu` with `h(mu) = 0`.
- **(c) Eigenvector sources.** If `x != 0`, `1^T x = 0` and `A x = mu x`, then `x^T K(beta, t) = (1 - beta - beta mu) x^T`
  for all `beta, t`, and `h(mu) = 0`. If `x != 0`, `1^T x = 0` and `A^T x = mu x`, then
  `K(beta, t) x = (1 - beta - beta mu) x` for all `beta, t`, and `h(mu) = 0`.
- **(d) Twins.** Call `i != j` (in sigma) *in-twins* if `A_ki = A_kj` for every `k` in sigma other than `i, j` (same
  in-neighbours apart from each other), and *out-twins* if `A_ik = A_jk` for every such `k`. If `i, j` are in-twins or
  out-twins **and** the pair is symmetric (`A_ij = A_ji =: c`), then `x = e_i - e_j` is as in (c) with `mu = -c`:
  - `c = 0` (**non-adjacent twins**): `K(beta, t)` has the eigenvalue `1 - beta` for every `beta, t`; equivalently
    `-I + W_sigma` has the eigenvalue `+delta`. So sigma fails (3) at **every** `(eps, delta)` with `eps > 0, delta > 0`
    (legal or not).
  - `c = 1` (bidirected twins): `K` has the eigenvalue `1` for every `beta, t` (harmless). Geneson's paper nodes 2, 3
    are bidirected in-twins; this is the eigenvalue `-1` of `A|perp` (the zero eigenvalue of `G` in lemmas.md §4).
  - A one-way pair (`A_ij != A_ji`) gives nothing: `e_i - e_j` is then not an eigenvector.
- **(e) Modules.** Call `S ⊆ sigma` an *in-module* if every `k` in `sigma \ S` sends to all of `S` or to none of `S`
  (`A_ki` constant in `i in S`), and an *out-module* if every `k` in `sigma \ S` receives from all of `S` or from none
  (`A_ik` constant in `i in S`). If `S` is an in-module and `z in C^S` satisfies `A_S z = mu z`, `1_S^T z = 0`, `z != 0`
  (`A_S` the adjacency induced on `S`), then `x := (z on S, 0 elsewhere)` satisfies `A x = mu x`, `1^T x = 0`. Same for
  out-modules with `A_S^T`. Twins are the case `|S| = 2`. Example: an in- or out-module inducing a **directed 3-cycle**
  gives `mu = e^(±2 pi i/3)`, `Re mu = -1/2`.
- **(f) The prune.** Assume `0 < eps < 1`, `delta > 0` (so `beta > 1`).
  - If `h` has a root `mu` with `Re mu >= -gamma(beta)`, then (3) fails at `(beta, t)` for every `t > 0`.
  - If `h` has a root with `Re mu >= -1/2`, then (3) fails at **every legal** `(eps, delta)`; so sigma is never an SFP
    support anywhere in the legal region. (Directed 3-cycle modules are of this kind.)
  - If `h` has a root with `Re mu >= 0` (e.g. non-adjacent twins), (3) fails at every `(eps, delta)` with `delta > 0`.
  - Every root of `h` is an eigenvalue of `A|perp`, so a sigma killed by the second bullet already has
    `alpha(A|perp) >= -1/2` (it is dead in the limit, I.3(d)). The content of (f) is that for these sigma the limit
    verdict holds **at every legal parameter**, with no finite-`t` analysis.

*Proof.*
- (a) By I.3(g), `K = (1 - beta) I - beta (A^T - tau J)`. So
  `lambda I - K = beta [ (lambda - 1 + beta)/beta I + A^T - tau J ] = beta [ -nu I + A^T - tau J ]`, and
  `det(lambda I - K) = beta^m (-1)^m det(nu I - (A^T - tau J)) = (-beta)^m chi_(A - tau J)(nu)` (a matrix and its
  transpose have the same characteristic polynomial, and `J^T = J`). By Corollary I.1.1 with `-tau` in place of `tau`,
  `chi_(A - tau J) = chi_A + tau m chi_(A|perp)`. (For `m = 1`: `A = (0)`, `K = (1 + t)`, and the right side is
  `-beta (nu + tau) = -(1 - beta - lambda) - beta - t = lambda - 1 - t = chi_K(lambda)`.) ∎
- (b) `m chi_(A|perp) = chi_A - chi_(A+J)` (Corollary I.1.1), so the two gcds agree. If `h(mu) = 0` then by (a) the
  bracket vanishes at `nu = mu` for every `tau`, i.e. `lambda = 1 - beta - beta mu` is a root of `chi_K` for every `t`.
  Conversely, for fixed `beta` the bracket is affine in `tau` (and `tau` is injective in `t`); if it vanishes at the
  same `nu` for two values of `tau`, both `chi_A(nu)` and `chi_(A|perp)(nu)` vanish. ∎
- (c) By I.3(g), `x^T K = (1 - beta) x^T - beta x^T A^T + beta tau (1^T x) 1^T = (1 - beta - beta mu) x^T` when `A x = mu x`
  and `1^T x = 0`; similarly `K x = (1 - beta) x - beta A^T x + beta tau 1 (1^T x)`. For `h(mu) = 0`: `mu` is an
  eigenvalue of `A` (resp. `A^T`, same characteristic polynomial); and with `y := U^T x != 0` (`x in 1^perp`, so
  `U y = U U^T x = x`): `(A|perp) y = U^T A x = mu y` (resp. `(A|perp)^T y = U^T A^T x = mu y`), so `mu` is an eigenvalue of
  `A|perp` too. ∎
- (d) For in-twins and `x = e_i - e_j`: `(A x)_k = A_ki - A_kj = 0` for `k != i, j`; `(A x)_i = A_ii - A_ij = -c`;
  `(A x)_j = A_ji - A_jj = c`. So `A x = -c x`. Out-twins: the same computation with `A^T`. For `c = 0` the eigenvalue
  `1 - beta < 0` (as `beta > 1` iff `delta > 0`), and `-I + W_sigma = -eps K` has the eigenvalue `-eps(1 - beta) = delta`. ∎
- (e) For `k` outside `S`: `(A x)_k = sum_(i in S) A_ki z_i = A_(k,S) (1_S^T z) = 0` (constant row segment); for `i in S`:
  `(A x)_i = (A_S z)_i = mu z_i`. And `1^T x = 1_S^T z = 0`. Out-modules: apply this to `A^T` (an out-module of `G` is an
  in-module of the reversed graph). For the 3-cycle `A_S = ` cyclic permutation, `z = (1, w, w^2)` with `w = e^(2 pi i/3)`
  (or its conjugate) is an eigenvector with `1 + w + w^2 = 0`. ∎
- (f) By (b), `1 - beta - beta mu` is an eigenvalue of `K` at every `t`, and `Re(1 - beta - beta mu) <= 0` iff
  `Re mu >= -gamma`; by (S2) (3) fails. Legal points have `beta > 2`, so `-gamma < -1/2 <= Re mu` gives strict `< 0`.
  `Re mu >= 0` and `beta > 1` give `Re(1 - beta - beta mu) <= 1 - beta < 0`. The last bullet: roots of `h` are roots of
  `chi_(A|perp)`. ∎

**Deciding (f) exactly.** `h` is computed over `Q` from two integer characteristic polynomials. Whether some root has
`Re mu >= -1/2` is "some irreducible factor `f` of `h` has `f(z - 1/2)` not Hurwitz", decided exactly by Routh-Hurwitz
(lemmas.md B.1 (b0), both directions, CITED). *Policy note (for the orchestrator, not a change to CLAUDE.md §3):* here a
vanishing Hurwitz minor of `f(z - 1/2)` also means "not Hurwitz" (lemmas.md (b0), (b1)), hence a root with
`Re mu >= -1/2`, and the conclusion of (f) only needs the non-strict inequality, because legality makes `gamma > 1/2`
strict. So the boundary case is mathematically decisive for this prune. Whether production code may use it is a policy
decision; the conservative option is to prune only on a strictly negative minor/coefficient, or on an exact root of a
low-degree factor of `h` (typically `nu`, `nu + 1`, `nu^2 + nu + 1`).

**Literature status.** I have not checked whether the twin/module statements are already in the CTLN literature
(they resemble known "graph rules"). Treat novelty as UNVERIFIED; the proofs above stand on their own.

---

## 9. Summary: the limit rule, and how this file relates to the existing notes

**The limit rule for the full support (PROVED; I.3, I.4).** Fix `beta > 2`, `gamma = 1 - 1/beta`, and let
`eps -> 0+` along `delta = (beta - 1) eps` (so `q -> beta - 1`). Let `alpha = alpha(A|perp)` and
`w = (A^T + gamma I)^-1 1` (defined at every rational `beta`, I.4(g)).

| condition | conclusion |
|---|---|
| `w > 0` and `alpha < -gamma` | SFP at every legal point of the ray with `t > T(beta)`; `T` uniform on compact subsets (I.3(f)). A **counterexample** iff the graph is not complete. |
| `w` not `> 0` (or `M_gamma` singular) | (1) fails at **every** `t > 0` on the ray (I.4(c)); no SFP anywhere on it. |
| `alpha > -gamma` | (3) fails for all large `t` (I.3(c2)); no conclusion at moderate `t`. |
| `alpha = -gamma` (i.e. `beta = beta_c`) | undecided by the limit (I.3(c3)). |

Equivalently: the limit-SFP `beta`-set is `{ beta in (2, beta_c) : w(gamma(beta)) > 0 }`, `beta_c = 1/(1 + alpha)`.

**Graph-level necessary conditions for limit stability at a legal `beta`** (each PROVED):
`|E| > m(m-1)/2` (I.5); not oriented (I.6.1); no root `mu` of `gcd(chi_A, chi_(A+J))` with `Re mu >= -1/2`, in particular no
non-adjacent twins and no 3-cycle module (I.7 -- and these are fatal at *every* legal parameter, not only in the limit).
**Sufficient conditions** from the symmetric part exist only for the complete graph (I.6(d)).

**Relation to `notes/limit_spectral.md`.**
- Claim 1 = I.2 (with I.1). Proved, in a sequence form that also gives uniformity.
- Claim 2 = I.3(a), (b), (e). Proved as stated, plus the boundary case (c3) and the existence form (d).
- Claim 3 = I.4. Proved; the nondegeneracy hypothesis "C invertible and `1 - beta 1^T u != 0`" is unnecessary (I.4(b2),
  (c)), and positivity is exactly `t`-free, not just a limit statement. "`adj(C) 1` one-signed iff `C^-1 1` one-signed" is
  correct as stated: one-signedness of either forces `det C != 0` (I.4(c)).
- "Consequences": the density bound = I.5 (correct); the complement form "all eigenvalues of `Abar|perp` have
  `Re > -1/beta`" is correct (`A|perp = -I - Abar|perp` since `J|perp = 0`); "a non-adjacent pair is survivable in a larger
  graph" is correct and is made precise in both directions: non-adjacent **twins** are fatal at every parameter (I.7(d)),
  while a non-adjacent pair in general only defeats the Bendixson sufficient test (I.6(d)).

**Relation to `notes/lemmas.md`.** Nothing there is contradicted.
- I.1 = C'(a) (same idea). I.2 = C(a) + C'(c), proved here with the elementary root-continuity lemma (R) instead of
  Rouché, and in sequence form. I.3(a)-(e) = C'(d), (e) + (I.2) there. I.5 = C'(f) and step 4 of (I.2) there.
- I.4 = (I.3), (I.4) there, sharpened: explicit `t`-free formulas for `x`, `det C`, `det K` in terms of `w`; the sign rule
  (I.4(e)); the case analysis of `C u = 1` (I.4(b)); nonsingularity of `M_gamma` at rational `beta` (I.4(g)).
- Lemma C(b) there uses a threshold `t > |d/s|`; by I.4(c) it is unnecessary (whenever (1) holds, `d` and `s` share a sign,
  so `D = d + t s != 0` for all `t > 0`). Lemma C(d)'s last sentence ("(1) holds for all large `t` iff all `a_i` nonzero with
  a common sign") is correct under the reading (P1) (`K` nonsingular); if `s = d = 0` then `K` is singular for every `t`
  and the weak reading of (1) is not covered by it -- irrelevant for SFP since (3) fails (I.4, Remark).
- The caveat "Uniformity" in lemmas.md §3 is sharpened by I.3(f) and the combined statement after I.3: `T` is uniform on
  compact subsets of the open limit-SFP `beta`-set; it is still not explicit.

**Relation to `notes/lemma_F_reduction.md`.** I.4(c) (P1) <=> (P2) is Lemma H(H3), reproduced; I.4(f) is (H4) in `w`-form.

---

## 10. Not proved here / open

- **`beta -> inf` jointly with `eps -> 0`** (Remark I.3.5): whether (3) at legal points with small `eps` but large `beta`
  forces anything about `alpha(A|perp)`. OPEN. On rays `t = r beta` the relevant quantity is `alpha(A - (1 + r) J)`, not
  `alpha(A|perp)`.
- **An explicit `T(beta)`** (an explicit `eps_0`) for the limit-stable cells. I.3(f) gives uniformity on compacts but no
  number. The finite-`t` route is lemmas.md Lemma B at fixed rational `beta`.
- **`beta = beta_c`.** Undecided by the limit (I.3(c3)).
- **Single-interval question** (lemmas.md §9: is the stable `t`-set at fixed `beta` an interval `(t*, inf)`?). By I.3(g),
  (3) at `(beta, t)` iff `phi(1 + t/beta) < -gamma` with `phi(tau) := alpha(A - tau J)`. If `phi` were nonincreasing the
  answer would be yes for (3), hence for SFP (I.4 makes (1), (2) `t`-free). But `phi` is **not** monotone in general.
  Exact example: `A = [[0,0,0],[1,0,1],[1,1,0]]` has `chi_(A - tau J)(l) = (l + 1)(l^2 + (3 tau - 1) l + tau)`, so
  `phi(tau) = (1 - 3 tau + sqrt((3 tau - 1)^2 - 4 tau))/2`, increasing from `-1` at `tau = 1` to `-1/3` as `tau -> inf`.
  For this graph the stable part `{phi < -gamma}` lies entirely outside the legal range (checked exactly at sample
  points; it is also a 3-node graph, settled by Geneson Thm 4.1). A float scan finds `phi` increasing somewhere on
  `[1, 10^4]` for 3283 of 9608 `n = 5` graphs, the largest rises just above `tau = 1` (small `t`, mostly illegal). Whether
  the stable `t`-set is an up-interval **on the legal range** remains OPEN.
- **Novelty of I.7** (twins/modules) against the CTLN literature: not checked.
- **Lean.** Good first targets, all purely algebraic (no analysis): (T1), (T2), I.1, Corollary I.1.1, (S4), I.4(c) and its
  formulas, I.7(a)-(e). I.2 needs (R) (continuity of roots); I have not checked what Mathlib provides for it.
- **Policy question** (I.7, "Deciding (f) exactly"): may a vanishing Hurwitz minor of `f(z - 1/2)` be used to prune?
  Mathematically yes; CLAUDE.md §3 is a policy I have not changed.

---

## Appendix: sanity checks (tests, not proofs) and exploration

Scratch scripts are in `/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverI/v2/`
(session-specific; copy into the repo if they are to be kept). All `s_*` checks use exact rationals except `s_I2.py`
(mpmath, 60 digits). `x_*` files are explorations, labelled as such; nothing above depends on them.

| script | what it checks | result |
|---|---|---|
| `s_I1.py` | I.1 symbolically (generic `X`, `m = 2, 3`, three unit vectors `u`); I.1a and Cor. I.1.1 exactly for `m = 4, 5` | 26/26 pass |
| `s_I2.py` | I.2: bounded eigenvalues of `C + tJ` vs `spec(C|perp)`, and `lambda* - m t - kappa`, 60 random graphs, `t = 1e6, 1e10` | max errors `4.9e-4 / 7.3e-6` at `1e6`, `4.9e-6 / 7.3e-10` at `1e10` (the `t^(-1/2)` rate of the first reflects a defective eigenvalue of `C|perp` in some sample; no rate is claimed in I.2) |
| `s_I3.py` | I.3(a) and (g) identities exactly; limit verdict (`chi_(A|perp)(z - gamma)` Hurwitz) vs exact Hurwitz of `det(mu I + K)` at `t = 1e9` | 250 cases; 249 agree, 0 disagree, 1 limit-boundary skipped; 84 limit-stable |
| `s_I4.py` | I.4(I.4.0), (b1), (c) P1<=>P2<=>P3<=>P3' against `ctln.analyze_support`, all formulas, (e), (f) | 600 cases pass; 326 positive; 34 cases with `w > 0`, `a < 0`, each with `det K < 0` and `ctln` stable = False; 652 off-support values exact |
| `s_I4_geneson.py` | Geneson's graph: `chi_(A|perp)`, `alpha`, `beta_c`, `w(29/54)`, pairs, `spec sym(A|perp)` | below |
| `s_I7.py` | I.7(a) symbolically in `beta, t` (`m <= 4`); (b) at random rational `(beta, t)` for every root of `h`; (f) vs `ctln` | 300 graphs (103 with planted twins); 193 symbolic identities; `h` nontrivial in 255; 113 killed graphs, 452 legal points, never stable |
| `x_I7_prune.py` | EXPLORATION: prune counts, full support | see below |
| `x_half.py`, `x_monotone.py`, `x_monotone_exact.py` | EXPLORATION: §7 boundary counts; non-monotonicity of `phi` | see §7 (I.6.1) and §10 |

**Geneson's graph** (test_ctln.py `GENESON_N`; `A[b][a] = 1 - N[a][b]`; 0-based nodes): `|E| = 19`, `e = 11`.

- `chi_(A|perp)(l) = (l + 1)(6 l^4 + 13 l^3 + 19 l^2 + 13 l + 4)/6`. Eigenvalues (decimals for orientation only):
  `-1`, `-0.544938 ± 0.358105 i`, `-0.538395 ± 1.130501 i`. So `alpha = -0.538395144311666...`,
  `beta_c = 1/(1 + alpha) = 2.16635502784914...` (`q_c = 1.16635...`); the lemmas.md sextic vanishes there (`~1e-26`).
  This agrees with lemmas.md §4 (there `G = -I - (A|perp)^T`; its quartic `6z^4 + 11z^3 + 16z^2 + 10z + 3` is the image of
  the quartic above under `z = -1 - l`).
- The eigenvalue `-1` comes from bidirected in-twins (paper nodes 2, 3 = 0-based 1, 2), I.7(d).
- Non-adjacent pair: 0-based (3, 4) = paper nodes {4, 5}. `spec sym(A|perp) ~ {-1.2071, -1.2007, -0.7274, -0.2385, +0.2071}`:
  `lambda_max >= 0` as I.6(d) predicts, while `alpha < -1/2`.
- At `beta = 54/25`: `w = (46649574, 73946250, 73946250, 99502074, 12637350, 145064574)/272447123 > 0`,
  `W1 = 451746072/272447123`, proportional to the paper's `v`; `(-1)^5 det M_gamma > 0`;
  `x_inf/theta = w/W1 = (863881, 1369375, 1369375, 1842631, 234025, 2686381)/8365668`.

**Prune exploration (`x_I7_prune.py`; counts only, not certificates; full support; "kill" = some root of
`h = gcd(chi_A, chi_(A+J))` has `Re >= -1/2`, decided by exact Hurwitz on each irreducible factor, counting a vanishing
minor as a kill per the policy note in §8).**

| set | graphs | domination-free | I.7 kills | of which via vanishing minor only | I.7 kills among domination-free | limit-dead (`alpha >= -1/2`) | domination-free and limit-alive |
|---|---|---|---|---|---|---|---|
| `n = 5`, all non-complete | 9607 | 1433 | 5295 | 389 | 896 | 9514 | 0 |
| `n = 6`, random sample | 3000 | 400 | 1502 | 69 | 211 | 2997 | 0 |

So at `n = 5` the I.7 prune removes, at every legal parameter, 896 of the 1433 domination-free graphs (4906 kills overall
without using vanishing minors). The `n = 5` row is consistent with the certified `n = 5` result (only `K5` is ever
stable). In the `n = 6` sample only 3 of 3000 graphs are limit-alive.

**Standalone re-run (sympy only) of the Geneson numbers and the I.7 test:**

```python
import sympy as sp
Nl = [[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,0,0,1],[1,1,0,0,1,0],[0,0,1,1,0,0],[0,0,0,1,1,0]]   # test_ctln.GENESON_N
m = 6; l = sp.symbols('l'); I, J = sp.eye(m), sp.ones(m, m)
A = sp.Matrix(m, m, lambda i, j: 0 if i == j else 1 - Nl[j][i])         # A[i][j] = 1 iff i -> j
chiA, chiAJ = (l*I - A).det(), (l*I - A - J).det()
chiperp = sp.factor(sp.expand((chiA - chiAJ)/m))                          # Cor. I.1.1: m chi_{A|perp} = chi_A - chi_{A+J}
print(chiperp)                                                            # (l+1)(6l^4+13l^3+19l^2+13l+4)/6
alpha = max(sp.re(r) for r in sp.Poly(sp.expand(chiperp), l).nroots(n=30)); print(alpha, 1/(1 + alpha))
b = sp.Rational(54, 25); g = 1 - 1/b
w = (A.T + g*I).LUsolve(sp.ones(m, 1)); print(list(w), sum(w))            # I.4: w > 0, W1 > 1
h = sp.gcd(sp.Poly(chiA, l), sp.Poly(chiAJ, l)); print(h)                 # I.7: h = l + 1 (bidirected twins), harmless
```
