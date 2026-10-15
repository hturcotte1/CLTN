
---

## 12. A second mechanism at the same point: k-fold self-joins (PROVED; instances CERTIFIED)

Added in the continuation session. Script `selfjoin.py` writes `selfjoin_certificate.json` and `cert_selfjoin_R_eps1_21.json`.

**Definition.** For a digraph G on m nodes and k >= 1, `G^{vk}` is the disjoint union of k copies of G, with every pair of
nodes in different copies joined in both directions; it has `n = km` nodes. It is not a blow-up of a fixed core, since each
copy keeps all of G's internal non-edges. So Lemmas U and U'' do not apply to it. It has no clique class that sends to
everything, and the non-edge degree of every node stays bounded in k. Clique-class families are different here: in `R_k`, node 2
has k non-edges.

**Lemma SJ.** Let G have non-edge matrix D and let `w_G(beta) = (A^T + (1 - 1/beta) I)^{-1} 1`.
- (a) `spec(D' + aJ_{km}) = spec(D + akJ_m) u spec(D)^{(k-1)}`, counted with multiplicity.
- (b) The full support of `G^{vk}` is positive at `(eps, delta)` iff `w_G(beta) > 0`, for every k and t.

*Proof.*
- (a) `D' = I_k (x) D` and `J_{km} = J_k (x) J_m`. Take an orthogonal U with `U^T J_k U = diag(k, 0, ..., 0)`. Conjugating
  `D' + aJ` by `U (x) I_m` gives `blockdiag(D + akJ_m, D, ..., D)`. `selfjoin.py` (J3) checks the identity
  `det(zI - aJ - D') = det(zI - D - akJ) chi_D(z)^{k-1}` exactly on the explicit graphs for k = 2, 3.
- (b) `A'^T + gamma I` (`gamma = 1 - 1/beta`) is nonsingular at rational beta. Its determinant vanishes only if `-gamma` is a
  rational eigenvalue of an integer matrix, hence an integer, and `gamma` lies in `(1/2, 1)`.
  - By Lemma I.4 (as in X2'), positivity is equivalent to `w' := (A'^T + gamma I)^{-1} 1 > 0`.
  - Permuting copies is an automorphism, so by uniqueness w' is the same vector `w_1` on every copy.
  - Every node of another copy sends to every node of copy i. So the rows of copy i read
    `(A^T + gamma I) w_1 + (k-1) s 1 = 1`, with `s = 1^T w_1`. Hence `w_1 = (1 - (k-1)s) w_G`, and
    `s = W_G / (1 + (k-1) W_G)`, so `w_1 = w_G / (1 + (k-1) W_G)`.
  - If `w_G > 0`, this is positive. Conversely, `w_1 > 0` forces `w_G` to be one-signed.
  - `w_G < 0` is impossible. For `v = -w_G > 0`, `(D^T + rI) v = (1 + V) 1` with `V = 1^T v`, and summing gives
    `sum_j outdeg_D(j) v_j = m + (m - r) V > (m-1) V`. But every out-degree is at most `m - 1`, so the left side is at most
    `(m-1) V`, which is a contradiction. ∎

**Theorem SJ.** Put `g(z, c) := z^5 - 3z^2 - 2z - 1 - c (7z^4 + 12z^3 + 21z^2 + 13z + 5)`, so that
`det(zI - D_R - cJ) = z^2 g(z, c)` for `R = &F^U|m^QTq?`.
- For every real `c >= 0`, every root of `g(., c)` has `Re > -10/21`. The five Hurwitz minors of `-g(-10/21 - x, c)` are
  polynomials in c with **all coefficients positive**; they are listed in `selfjoin_certificate.json`.
- Consequently, for every legal `(eps, delta)` with `delta <= (11/10) eps` (i.e. `beta <= 21/10`) and every `k >= 1`, the full
  support of `R^{vk}` (n = 7k) is a stable positive fixed point. This includes `(1/21, 11/210)`.

*Proof.*
1. Legal with `beta <= 21/10` gives `r = 1/beta >= 10/21`.
2. By Lemma SJ(a), the spectrum is the roots of `z^2 g(z, ak)` together with those of `z^2 g(z, 0)` (the latter only when k >= 2).
   All of these have `Re > -10/21 >= -r`, and the roots `z = 0` are harmless.
3. Lemma S then gives stability.
4. Positivity: `w_R > 0` on `[2, 21/10]` (exact Sturm counts, J2), then Lemma SJ(b).
5. R is an induced subgraph, so no member is a clique. ∎

At `beta1 = 43/20` and `11/5` the locus condition fails. This is already the case at `c = 0`, because
`min Re spec(D_R) ~ -0.4710 < -20/43`. So the self-join's region is smaller than the clique-class region of Theorem A+.

CERTIFIED instances (`cert_selfjoin_R_eps1_21.json`), via `ctln.analyze_support` at `(1/21, 11/210)`. Each is legal, has
`x > 0`, has every Hurwitz minor of the exact characteristic polynomial positive, and is not a clique:

| k | 1 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|
| n | 7 | 14 | 21 | 35 | 70 |

(These are additional CANDIDATE counterexamples pending human hand-verification. Nothing was built on them.)

**Remark (one polynomial, two families).** Theorem A's polynomial is `g_k(z) = g(z, a) - a k (z^4 + z^3 + 3z^2 + z + 1)`. The
self-join uses `g(z, ak)` together with `spec(D)`. Both families are governed by the same core pencil `chi_D - c psi_D`. The
clique class adds `ak chi_L / z^3`; the self-join instead rescales `c = a -> ak`.

**Lemma M2 (multiple eigenvalues survive; PROVED).** If `lambda` is an eigenvalue of `D'` with geometric multiplicity at least 2,
then `lambda` is an eigenvalue of `D' + aJ` for every a. *Proof:* the eigenspace meets `1^perp` nontrivially, and there
`(D' + a 1 1^T) v = lambda v`. ∎ Consequences:
- (i) A stable full support needs every eigenvalue of D' of geometric multiplicity >= 2 to have `Re > -r`.
- (ii) Self-joins: for k >= 2 the copies of `spec(D)` survive (Lemma SJ(a)). Positivity forces D to contain a directed cycle
  (§8), so `rho(D) >= 1`, `tr D = 0`, and `min Re spec(D) <= -1/(m-1)`. So **self-join families of an m-node core exist only
  at `beta < m - 1`**, the same cap as joins with `K_k`.
- (iii) Symmetric constructions do not raise ceilings. For a sink-class family over a core `G^{v2}`, the surviving copy of
  `spec(D_G)` caps the ceiling at G's join ceiling. New ceilings above `beta_Y` need asymmetric cores.

**Structural summary of every proved fixed-parameter family (R, Q, Q2, Y, joins, self-joins).**
- The nonzero spectrum of the non-edge matrix `D'_k` is a fixed finite set, independent of k up to multiplicity. For sink
  classes `spec(D') = spec(D) u {0}`; for self-joins it is `spec(D)` repeated.
- `D'_k + aJ` has exactly one eigenvalue escaping to `+inf` (the Perron direction, carried by `1`). All other eigenvalues
  converge to a k-independent finite set.
- The failing families are the opposite case. Uniform blow-ups have `spec(D'_k) = N spec(...)` growing linearly (Lemma U),
  and F has its fixed limit set on `Re = -1/2` (FP(ii)).
