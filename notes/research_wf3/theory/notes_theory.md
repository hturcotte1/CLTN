# General-n theory of full-support counterexamples (workflow wf3, label `theory`)

Date 2026-10-14. Conventions are exactly CLAUDE.md §2: `A[i][j] = 1` iff `i -> j`; `N = J - I - A^T`;
`K = I + beta N + t J`, `t = (1-eps)/eps`, `beta = 1 + delta/eps`, `gamma = 1 - 1/beta`, `r := 1/beta`;
legal iff `t(beta-2) > 1`. Nothing in the repository was modified. Every script named here lives in
`/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf3/theory/` (below: `THEORY/`)
and re-runs standalone with `python3 <script>` from that directory (they import the repo read-only).

Labels: **PROVED** (proof written here, or cited from a repo note that marks it proved); **CERTIFIED** (exact
rational/integer computation, script and output file named); **HEURISTIC** (floating point, used only to steer
searches; never a decision); **OPEN**; **REFUTED**.

> **Hand-verification status (CLAUDE.md §3).** Every graph on n >= 7 nodes named below as a counterexample is a
> **CANDIDATE pending human hand-verification**. The infinite-family statements are proved by exact algebra plus
> lemmas cited from `notes/lemma_I_spectral.md` (I.2, I.3, I.4) and `notes/lemma_M_refuted.md` (identity (E),
> legality formula (2), Lemma AS). Their conclusions are still candidate counterexamples, and nothing here should be
> built on until a human has checked at least the base instances. The n=7 bases used below (`&FNDL~[uds?`,
> `&F^S{^ReL|?`, `&F^V\K[epa?`, the 123-graph list) are themselves unverified candidates. The proofs below do **not**
> use the fact that they are counterexamples; they use them only as explicit matrices.

---

## 0. Summary

| # | statement | status | where |
|---|---|---|---|
| 1 | Non-edge picture: with `D := J - I - A`, `A|perp = -I - D|perp`; limit-stable at beta iff every eigenvalue of `D|perp` has Re > -1/beta; positivity iff `u = (D^T + rI)^{-1} 1 > 0` and `1^T u > 1` | PROVED | §1 |
| 2 | Lemma W: closed-twin blow-ups. Spectrum `n chi_{D'|perp} = z^{n-m} mu^T adj(zI - DM) 1`; fixed point class-constant with `((A^T+I)M + (gamma-1)I) y = 1` | PROVED; CHECKED on explicit graphs | §2, `weighted.py` |
| 3 | Lemma X1/X2: clique-class extension `G[S,T;k]`: `Phi_k = z psi_D + k chi_L`, `L = D - b1^T - 1c^T`; for T = all, positivity is k-free: `w_G > 0` and `sum_{i in S} w_G,i > 1` | PROVED; CHECKED (30 random + all families) | §2, `classext.py` |
| 4 | **Theorem F.** `F_k` (Geneson's graph + bidirected clique class of size k), n = 6 + k: a non-clique full-support SFP at legal parameters for **every** n >= 6. Limit q-sup 1.1664 (n=6), 1.2747 (n=7), **the plastic number 1.3247 for n = 8..14**, then decreasing to 1 | PROVED (closed-form Hurwitz data, all k); instances n = 7..16 CERTIFIED | §3, `family_F_proof.py` |
| 5 | Catalogue: 9 infinite families on the three n=6 hits and **617** on 117 of the 124 n=7 bases, each proved for every k >= 1 (all but one also for k = 0) | PROVED by exact computation | §3.4, `family_prover.py` |
| 6 | Mechanism of `&FNFL~[udo?`: it is the n=6 hit 1387528 with one node doubled; alpha = -1/2 is an exact root-locus crossing of the cube root of unity; first-order coefficient R* = 1/4 > 1/8 | PROVED + CERTIFIED | §4 |
| 7 | Lemma B (boundary mechanism): alpha = -1/2 from a simple moving cube-root pair, R* > 1/8, positivity at beta = 2  =>  SFP at legal points with eps -> 0, never limit-stable | PROVED (from Lemma AS, (E), (2)) | §4.2 |
| 8 | **Theorem E.** `E_a` (n = 7 + a, a >= 0; `E_0 = &FNFL~[udo?`): finite-eps-only counterexample for **every** n >= 7. Closed forms: `Phi = z (z^2+z+1)((a+7) z (z^2+z+1) + 2)`, R* = 1/4 for all a, positivity vector at beta = 2 equal to `(6,2,2,18,6,2,10)/(18a+17)` | PROVED; instances n = 7..12 CERTIFIED | §4.3, `family_E.py` |
| 9 | 15 further n=8 finite-eps-only candidates (one-node extensions of n=7 candidates), alpha = -1/2 exactly, R* in {1/6, 1/4, 1/3, 1/2}, each with an exact legal SFP witness | CERTIFIED | §4.4, `cube_n8_witnesses.json` |
| 10 | Lower bounds on full-support q_n* at exact rational points: n=8: 1709/1001 (1.7073); n=9: 17691/10010 (1.7673); n=10: 18418/10001 (1.8416); n=11: 17220/9091 (1.8942) | CERTIFIED | §5.1 |
| 11 | **Theorem Q.** `Q_k = B[S,all;k]` (B = `&F^S{^ReL|?`): positive at every beta > 2; limit-stable at beta = 5/2 for all k >= 25, at 27/10 for all k >= 79, at 14/5 for all k >= 286. Hence for every n >= 32 (86, 293) a full-support counterexample with q arbitrarily close to 3/2 (1.7, 1.8) exists. Limit q of the family = 1.84416... | PROVED; instance n = 32 CERTIFIED | §5.2, `family_Q_proof.py` |
| 11b | **Theorem Q2.** Two sink classes (sizes k, 3k) on the 7-node core `&F^U|m^QTq?`, n = 7 + 4k: positive at every beta > 2; limit-stable at beta = 3 (resp. 7/2, 19/5) for all k >= 3 (10, 33). So for every n = 7 + 4k >= 19 (47, 139) there are counterexamples with q arbitrarily close to 2 (5/2, 14/5). Instances n = 19 (q = 20000/10001) and n = 47 (q = 25000/10001) | PROVED; CERTIFIED | §5.2b, `family_Q2_proof.py` |
| 11c | **Theorem Y.** Two sink classes (sizes 4k, k) on the 8-node core `&G^ZZZr|KjqTo`, n = 8 + 5k: positive at every beta > 2; limit-stable at beta = 4 (9/2, 24/5) for all k >= 46 (117, 363): q arbitrarily close to 3 (3.5, 3.8) for every n = 8 + 5k >= 238 (593, 1823); family limit q -> 3.96. Instance n = 73 (beta = 3, q = 20000/10001) | PROVED; CERTIFIED | §5.2c, `family_Y_proof.py` |
| 12 | Scale law for weighted blow-ups: (mu, beta) feasible iff (N mu, beta/N) feasible | PROVED | §5.3 |
| 13 | Conjecture P: if `D != 0` and `(D^T + rI)^{-1}1 > 0` for some r in (0, 1/2), then D has an eigenvalue with Re <= -1/2. Consequences: joins `G v K_k` fail for large k, `G v G` never works | OPEN; HEURISTIC evidence (38,452 random + adversarial ascent, sup exactly -1/2) | §5.4 |
| 14 | "Sink-class limit spectrum always touches Re = -1/2" (conjecture SL) | REFUTED: 55 of the 135 distinct limit polynomials of the proved families have min Re > -1/2 (e.g. -0.3516) | §5.4 |
| 15 | Is sup_n q_n* finite? | OPEN; it is at least 3.8. Certified: 1.166, 1.325, 1.707, 1.767, 1.842, 1.894 (n = 6..11), 1.9998 (n = 19), 2.4998 (n = 47). Proved: q_n* >= 3.8 for every n = 8 + 5k >= 1823 (q arbitrarily close to 3.8 at legal points). Geneson's upper bound is n - 2 | §5.5 |

---

## 1. The non-edge picture (PROVED)

Let `D := J - I - A` (`D[i][j] = 1` iff `i != j` and `i -/-> j`; `D = N^T`), `e = 1^T D 1`.

**Lemma D1.** `A|perp = -I_{m-1} - D|perp`. Hence `spec(A|perp) = {-1 - zeta : zeta in spec(D|perp)}`, and
`alpha(A|perp) < -gamma` iff every eigenvalue of `D|perp` has `Re zeta > -(1 - gamma) = -1/beta`.
*Proof.* `J|perp = 0` (notes/lemma_I_spectral.md §0), `I|perp = I`. ∎

So the limit criterion (Lemma I.3) reads: **limit-stable at beta iff every root of
`psi_D(z) := 1^T adj(zI - D) 1 = m chi_{D|perp}(z)` (Lemma I.1) lies in `Re z > -r`, r = 1/beta.** The roots of `psi_D`
are the zeros of `T_D(z) = 1^T (zI - D)^{-1} 1` (moving) together with the common roots of `psi_D` and `chi_D` (fixed; e.g.
eigenvalues of D with a right eigenvector orthogonal to 1).

**Lemma D2 (positivity).** For rational beta > 1 (so `D^T + rI` is nonsingular, Lemma I.4(g)): `w = (A^T + gamma I)^{-1} 1 > 0`
iff `u := (D^T + rI)^{-1} 1 > 0` and `U := 1^T u > 1`; then `w = u/(U - 1)`.
*Proof.* `A^T + gamma I = J - (D^T + rI)`. If `(J - M) w = 1` with `M = D^T + rI` then `M w = (1^T w - 1) 1`, so
`w = (W - 1) u` and `W = (W - 1) U`, i.e. `W - 1 = 1/(U - 1)`. `u < 0` is impossible (`M >= 0` entrywise off the
diagonal and `r > 0`, so `Mu = 1` forbids `u < 0`). ∎

Consequences used below (PROVED, elementary): every node with a D-out-edge has `u_i < 1` and D-in-mass
`sum_{j -> i in D} u_j = 1 - r u_i > 1 - r`; a node with no D-in-edge has `u = 1/r`, so its D-out-neighbours get
`u < 0` unless it has no D-out-edge (this is graphical domination in the D language).

---

## 2. Blow-ups and clique-class extensions (PROVED)

### 2.1 Lemma W (closed-twin blow-up)

H on m nodes, weights `mu in Z_{>0}^m`, `M = diag(mu)`. `H^mu`: node a replaced by a bidirected clique `V_a` of size
`mu_a`; for a != b every node of `V_a` sends to every node of `V_b` iff `a -> b` in H. n = sum mu. Let `E` be the
n x m class-indicator matrix. Then `D' = E D E^T`, `E^T E = M`, `1_n = E 1_m`, `J_n = E J_m E^T`,
`A' = E (A + I) E^T - I_n`.

(W1) `n chi_{D'|perp}(z) = 1^T adj(zI - D') 1 = z^{n-m} mu^T adj(zI - DM) 1 =: z^{n-m} Phi_mu(z)`.

(W2) `det(zI + A') = (z-1)^{n-m} det(zI + A_Q)`, `1^T adj(zI + A') 1 = (z-1)^{n-m} mu^T adj(zI + A_Q) 1`,
`A_Q := (A + I) M - I`.

(W3) For gamma with `A'^T + gamma I` nonsingular: `w' = E y`, where `((A^T + I) M + (gamma - 1) I) y = 1`.

*Proof.* (W1) Push-through: `(zI - E D E^T)^{-1} E = E (zI - D E^T E)^{-1} = E (zI - DM)^{-1}` for z off the spectra, so
`1^T (zI - D')^{-1} 1 = 1^T E^T E (zI - DM)^{-1} 1 = mu^T (zI - DM)^{-1} 1`. Sylvester:
`det(zI_n - E D E^T) = z^{n-m} det(zI_m - D E^T E)`. Multiply; both sides are polynomials, so the identity holds for all z.
(W2) The same computation with `zI + A' = (z-1) I + E (A+I) E^T`. (W3) `(A'^T + gamma I) E y = E[(A^T+I) M y + (gamma-1) y]
= E 1 = 1`; uniqueness. ∎

CHECKED exactly against explicit blow-ups (`weighted.py` main; `family_E.py` (E1) for a = 1..8).

### 2.2 Lemma X1 / X2 (clique-class extension)

`G[S,T;k]`: add a class V of k nodes, bidirected among themselves; `i -> v` (all v in V) iff `S_i = 1`; `v -> i` iff
`T_i = 1`. `b := 1 - S`, `c := 1 - T`. This is Lemma W for `G + v` with `mu = (1,...,1,k)`.

**(X1)** `n chi_{D'|perp} = z^{k-1} Phi_k`, `Phi_k(z) = z psi_D(z) + k chi_L(z)`, `L := D - b 1^T - 1 c^T`.

*Proof.* With `D_+ = [[D, b],[c^T, 0]]`, `D_+ M = [[D, kb],[c^T, 0]]`, `R = (zI - D)^{-1}`. Schur:
`det(zI - D_+ M) = chi_D (z - k c^T R b)`. Solving `(zI - D_+M)(x; y) = (1; 1)` gives
`y = (1 + c^T R 1)/(z - k c^T R b)` and `mu^T (x; y) = 1^T R 1 + k y (1 + 1^T R b)`. Multiplying,
`mu^T adj(zI - D_+M) 1 = z psi_D + k chi_D [(1 + c^T R 1)(1 + 1^T R b) - (c^T R b)(1^T R 1)]`, and the bracket times
`chi_D` is `det(zI - D + b 1^T + 1 c^T)` by the rank-2 determinant lemma with `U = [b, 1]`, `V = [1, c]`. ∎
CHECKED: `classext.check_X1` on 30 random (G, S, T), k = 1..5, and on every family below.

**(X2) Positivity for T = all (c = 0).** For every k >= 1 and rational beta: the full support of `G[S,all;k]` is positive
iff `w_G(beta) > 0` and `sigma_S(beta) := sum_{i in S} w_G,i(beta) > 1`. In particular it does not depend on k.

*Proof.* By symmetry and uniqueness `w'` is constant (= y) on V. Rows of V: `(gamma + k - 1) y + sigma = 1` with
`sigma = sum_{i in S} w'_i`; rows of G: `(A^T + gamma I) w'_G = (1 - k y) 1`, so `w'_G = (1 - ky) w_G`. Solving,
`y = (1 - sigma_S)/(gamma - 1 + k(1 - sigma_S))` and `1 - ky = (gamma - 1)/(gamma - 1 + k(1 - sigma_S))`, where now
`sigma_S` is computed from `w_G`. As `gamma < 1`, `w' > 0` iff `w_G > 0` and `sigma_S > 1` (a one-signed `w_G` is positive,
Lemma I.4(P3')). ∎
Interpretation: `sigma_S > 1` says exactly that V would be **activated** by G's own fixed point (Lemma I.4(f)): the old
support fails its off-support condition at V, and the enlarged support takes over.

**(X3) The k -> infinity limit.** `Phi_k / k -> chi_L` coefficientwise with constant degree, so by root continuity
(notes/lemma_I_spectral.md (R)) the moving spectrum of `G[S,T;k]` converges to `spec(L)`. Several sink classes with in-sets
`S_j` and sizes `k_j` give `Phi = z psi_D + p chi_D + 1^T adj(zI - D) B`, `p = sum k_j`, `B = sum k_j b_j`; their limit is
`spec(D - bbar 1^T)` with `bbar` in the convex hull of the `b_j` (same proof).

---

## 3. Question 1: an infinite family of limit counterexamples

### 3.1 Theorem F (PROVED; conclusions are candidates pending hand-verification)

`G1` = Geneson's graph in the paper's labelling (notes/geneson_2607.21396.md; `n7/geneson.d6 = &EZjtwWE`):

    G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]

`F_k := G1[S,all;k]`, `S = (1,1,1,1,1,0)`: k bidirected clique nodes, each sending to all six G1 nodes and receiving from
paper vertices 1..5, not from vertex 6. n = 6 + k. `F_0 = G1`; `F_1 = &FNFL|\u|o?` (one of the 123 n=7 candidates).

**Theorem F.** For every k >= 0, `alpha(F_k|perp) < -1/2`, and `F_k` is positive at every beta in `(2, b_p)`,
`b_p - 1` = the plastic number (real root of `q^3 = q + 1`). Hence for every k >= 0 and every rational
`beta in (2, min(b_p, beta_c(k)))` there is T such that for all t > T the full support of `F_k` is a stable positive fixed
point at the legal point `eps = 1/(1+t)`, `delta = (beta-1)/(1+t)`. `F_k` is not a clique. So **for every n >= 6 there is a
full-support counterexample on n nodes** (candidate status as above).

*Proof.*
1. (X1) with the closed forms, verified exactly: `z psi_D = z^2 (6z^4 + 11z^3 + 16z^2 + 10z + 3)` and
   `chi_L = z^4 (z^2 + z + 1)`, so `Phi_k = z^2 R_k` with `R_k = (6+k) z^4 + (11+k) z^3 + (16+k) z^2 + 10 z + 3`.
2. `16 R_k(-1/2 - x) = 16(k+6) x^4 + 16(k+1) x^3 + 8(2k+17) x^2 + 12(k+1) x + (3k+16)` (exact), with
   Hurwitz data `a_1 = 16(k+1)`, `Delta_2 = 64(k+1)(k+16)`, `Delta_3 = 8192 (k+1)^2`, `a_4 = 3k+16`: all positive for
   k > -1. So every root of `Phi_k` has `Re zeta > -1/2` (zeros of the factor `z^2` and of `z^{k-1}` have Re 0), and by
   Lemma D1 `alpha(F_k|perp) < -1/2`. By Lemma I.3(e) the limit-stable set is `(2, beta_c(k))`, `beta_c(k) > 2`.
3. Positivity (X2): with `Mb = beta A^T + (beta-1) I`, `det Mb = -2b^5 + 7b^4 - 12b^3 + 10b^2 - 6b + 1`,
   `adj(Mb) 1 = (-(b-1)(b^3-2b^2+b-1), -(b^2-b+1), -(b^2-b+1), -(b-1)^2(b^2-b+1), b^3-3b^2+2b-1, -b^4+3b^3-4b^2+2b-1)`,
   `w = b adj(Mb)1 / det Mb`, `sigma_S - 1 = -(b^2-b+1)/det Mb`. Exact root counting on `[2, 93/40]`: `det` and every
   `adj_i` except `adj_4` have no root and are negative at 2; `adj_4 = b^3-3b^2+2b-1` has exactly one root `b_p`
   (`b_p - 1` solves `q^3 - q - 1 = 0`). So `w > 0` and `sigma_S > 1` on `(2, b_p)`, and positivity fails on
   `(b_p, 93/40]` (`adj_4 > 0 > adj_1`).
4. Lemma I.3 combined statement (lemma_I_spectral.md, after I.3): positivity is t-free (I.4), stability holds for all
   large t on compact subsets of the open limit-stable set. ∎

CERTIFIED (`family_F_proof.py` -> `family_F_certificate.json`): every assertion above, the explicit-graph check
`alpha < -1/2` for k = 0..15 (`tlib.alpha_less_than`, exact Hurwitz), and exact SFP witnesses by `ctln.analyze_support`:

| k | n | digraph6 | beta | t | eps | delta | q |
|---|---|---|---|---|---|---|---|
| 1 | 7 | `&FZt|^DaK~?` | 21/10 | 20 | 1/21 | 11/210 | 22/21 |
| 2 | 8 | `&GZy~VwuKW~^w` | 58/25 | 20 | 1/21 | 11/175 | 44/35 |
| 3 | 9 | `&HZ|^t~FWwo~^|~o` | 58/25 | 20 | 1/21 | 11/175 | 44/35 |
| 4 | 10 | `&IZ}n|^w|bp_~^}~|~w` | 58/25 | 20 | 1/21 | 11/175 | 44/35 |
| 5 | 11 | `&JZ~V~V~FuNb?~^~^~^~^~?` | 58/25 | 20 | 1/21 | 11/175 | 44/35 |
| 6, 8, 10 | 12, 14, 16 | see certificate | 21/10 | 20 | 1/21 | 11/210 | 22/21 |

All: legal, x > 0, all Hurwitz minors positive, not a clique.

### 3.2 q along the family (PROVED values; decimals for display only)

`beta_c(k)` = smallest root > 2 of the binding Hurwitz minor of `R_k(-1/beta - x)` (minimal polynomials in the
certificate, F5). Limit-sup of q = `min(b_p, beta_c(k)) - 1`:

| k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12 | -> inf |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 18 | |
| beta_c(k) | 2.1664 | 2.2747 | 2.3391 | 2.3720 | 2.3833 | 2.3809 | 2.3701 | 2.3549 | 2.3376 | 2.3197 | 2.3022 | 2.2700 | 2 |
| q-sup | 1.1664 | 1.2747 | **1.3247** | 1.3247 | 1.3247 | 1.3247 | 1.3247 | 1.3247 | 1.3247 | 1.3197 | 1.3022 | 1.2700 | 1 |

k = 0 reproduces Geneson's degree-6 polynomial `797b^6 - 6468b^5 + ...` exactly. The k -> infinity behaviour is forced by
`chi_L = z^4 (z^2+z+1)`: the limit spectrum touches `Re = -1/2` at the cube roots of unity, so `beta_c(k) -> 2`
(numerically `beta_c(k) ~ 2 + 4/k`).

### 3.3 How the n=7 candidates arise from n=6 hits (CERTIFIED by direct computation, `sub6.py`, `ext_patterns.py`)

Only 21 of the 123 n=7 candidates contain an n=6 hit as an induced subgraph (8 contain G1, 7 contain G2, 5 contain H3 =
idx 1387528, 1 contains G2 and H3); none contains K6. In 21 of the 23 embeddings the added node sends to **all** six
nodes (T = all); the in-pattern varies. `&FNFL|\u|w?` and `&FNFL~[u|w?` are `G1 v K1` and `G2 v K1` (bidirected universal
node). Joining larger cliques fails: exactly, `alpha(G1 v K_1) < -1/2` but `alpha(G1 v K_2) >= -1/2` (Hurwitz minor 4
negative). By (X3) the join limit is `spec(D)`, which for G1 contains `-0.5545 +- 1.1562 i` (display only); the general
statement is Conjecture P, §5.4.

### 3.4 Catalogue of provable infinite families (PROVED by exact computation, `family_prover.py`, `run_family_prover.py`)

For a base G and in-set S (T = all), `family_prover.prove_family` proves "for every integer k >= k0,
`alpha(G[S,all;k]|perp) < -1/2`" by showing that every Hurwitz minor of `(-1)^m Phi_k(-1/2 - x)` (a polynomial in k) has
positive leading coefficient and no real root in `[k0, inf)` (exact real-root isolation), and checks `w_G(2) > 0`,
`sigma_S(2) > 1` exactly (then positivity holds near beta = 2 by continuity, so for each k some beta in
`(2, min(beta_c(k), 2+eta))` gives an SFP for small eps).

- n=6 bases (`families_catalogue_n6.json`): 3 families per hit (G1, G2, H3), all proved for **all k >= 0**. Limit
  polynomials `z^4(z^2+z+1)`, `z^2(z^4+2z^3+4z^2+3z+1)`, `z^2(z^2+z+1)^2`.
- n=7 bases (`families_catalogue_n7.json`): of 128 x 124 (G,S) pairs, 617 pass a float prefilter, and **all 617 are then
  proved exactly** for every k >= 1 (616 also for k = 0); they live on 117 distinct bases. 132 distinct limit polynomials
  `chi_L` (135 together with the n=6 catalogue).

So "add a large bidirected clique that sends everywhere and is activated by the base" is a robust generator of
counterexamples at every n.

---

## 4. Question 2: the alpha = -1/2 boundary graph `&FNFL~[udo?` and finite-eps-only families

### 4.1 What FNFL is (CERTIFIED)

- Nodes 2 and 3 of `&FNFL~[udo?` are closed twins (bidirected, same in- and out-neighbours). Deleting either gives the
  n=6 hit idx 1387528 (H3). So **FNFL = H3 with one node doubled** (Lemma W with `mu = (1,1,2,1,1,1)`).
- `7 chi_perp = (l+1)(l^2+l+1)(7l^3+14l^2+14l+5)`: `l+1` is the closed-twin eigenvalue; `l^2+l+1` is a **moving**
  eigenvalue (not a root of `chi_A`; the root `z0 = 1 + omega` of `psi_A` is simple). In the D-picture the moving
  polynomial of the node-blow-up family of H3 is
  `Phi_k = z[(6+k) z^5 + (12+2k) z^4 + (18+3k) z^3 + (13+3k) z^2 + (6+3k) z + (1+k)]`, and the cube-root factor appears
  at exactly k = 1:  `Phi_1 = z (z^2+z+1)(7z^3+7z^2+7z+2)`. Modulo `z^2+z+1`, `z psi_D == -z^2` and `chi_L == +z^2`, so
  `Phi_k(omega) = (k-1) omega^2`: the root locus of the twin-class size passes through omega exactly at class size 2.
- There is **no 3-cycle module**: a 3-cycle module would give a FIXED eigenvalue `omega` of A, i.e. an eigenvalue
  `1 - beta/2 -+ i beta sqrt(3)/2` of K for every t, unstable for every beta > 2. FNFL's cube root is moving, which is
  what allows finite-t stability. Its eigenvector (D-picture, computed exactly) has full support.
- First-order coefficient (Lemma AS, notes/lemma_M_refuted.md §4): `R* = Re(chi_A(z0)/psi_A'(z0)) = 1/4` (`rstar.py`).

### 4.2 Lemma B (boundary mechanism) — PROVED

Let G be a digraph with `chi(z) = det(zI + A)`, `psi(z) = 1^T adj(zI + A) 1`. Assume
(i) `psi` has a conjugate pair of simple roots `z0, z0bar` with `Re z0 = 1/2` and `chi(z0) != 0`;
(ii) every other root of `psi` has `Re > 1/2`;
(iii) `R := Re(chi(z0)/psi'(z0)) > 1/8`;
(iv) `w(2) = (A^T + I/2)^{-1} 1 > 0`.
Then the full support is a stable positive fixed point at legal `(eps, delta)` with eps arbitrarily small, and
`alpha(A|perp) = -1/2`, so G is never limit-stable at a legal beta (Lemma I.3(d)).

*Proof.* `nu_A(c) := min Re spec(cJ - A)`. By (E) of lemma_M_refuted.md, `K = beta(cJ - A^T - gamma I)` with
`c = 1 + t/beta`, so K is positive stable iff `nu_A(c) > gamma`; by its eq. (2), for c > 1 the point is legal iff
`gamma > nu_P2(c) := c - sqrt(c^2 - c) = 1/2 + 1/(8c) + O(c^-2)`. As c -> inf, one eigenvalue of `cJ - A` tends to +inf and
the others to the roots of psi (Lemma I.2 with C = -A). By (ii) all but the two branches at `z0, z0bar` stay in
`Re > 1/2 + delta0`; by Lemma AS and (i) those branches are `z0 + 1/(c s'(z0)) + O(c^-2)` with `1/s'(z0) = chi(z0)/psi'(z0)`.
Hence `nu_A(c) = 1/2 + R/c + O(c^-2)` and, by (iii), `nu_A(c) - nu_P2(c) >= (R - 1/8)/(2c) > 0` for c >= c1. Take
`gamma(c) := (nu_A(c) + nu_P2(c))/2`: legal and stable, and `gamma(c) -> 1/2`, i.e. `beta(c) -> 2+`. By (iv) and
continuity (`A^T + gamma I` is nonsingular at the rational point beta = 2, Lemma I.4(g)), `w > 0` on `[2, 2 + eta)`;
positivity is t-free (I.4). For c large, `beta(c) < 2 + eta`: SFP at `t = beta(c)(c - 1)`, `eps = 1/(1+t) -> 0`. Finally,
the roots of `psi` are the numbers `-mu`, mu in `spec(A|perp)`, so `-z0` is an eigenvalue of `A|perp` with real part -1/2
and by (ii) all others have real part < -1/2: `alpha(A|perp) = -1/2`. ∎

This is exactly why the stable window of FNFL is bounded above in t and sits at beta slightly above 2.

### 4.3 Theorem E: an infinite family of finite-eps-only counterexamples (PROVED; candidates)

`H = &FNDL~[uds?` (an n=7 candidate) in the labelling

    H = [[0,0,1,1,1,1,0],[0,0,1,0,1,0,0],[1,1,0,1,1,1,1],[1,1,1,0,1,1,1],[0,0,1,1,0,1,1],[0,1,0,0,1,0,1],[1,1,0,1,0,0,0]]

`E_a := H^{mu(a)}`, `mu(a) = (1,1,2,a,1,1,1)` (a = 0 deletes node 3), n = 7 + a. `E_0 = &FNFL~[udo?`. Nodes 2 and 3 of H
send to everybody (D-sinks), so `E_a` is also the sink-class extension `FNFL[S, all; a]`.

**Theorem E.** For every a >= 0, `E_a` is not a clique, `alpha(E_a|perp) = -1/2` exactly (never limit-stable at a legal
beta), and the full support of `E_a` is a stable positive fixed point at legal `(eps, delta)` with eps arbitrarily small.
So **for every n >= 7 there is a finite-eps-only full-support counterexample** (candidate status as above).

*Proof (all closed forms CERTIFIED by `family_E.py` -> `family_E_certificate_part1.json`).*
1. (E2) Lemma W: `Phi_{mu(a)}(z) = z (z^2+z+1) ((a+7) z (z^2+z+1) + 2)` identically in a; the moving spectrum of `E_a` is
   `{0} u {omega, omegabar} u roots of the cubic`, plus `0` with multiplicity `n - 7` (closed twins).
2. (E3) With `c = a+7`, the cubic shifted by `z = -1/2 - y` has Hurwitz data `c, c/2, 3c/8 - 2` and `Delta_2 = 2c`: all
   its roots have `Re z > -1/2` iff `c > 16/3`, true for every a >= 0. In the A-picture (`z_A = 1 + zeta`) this is (ii) of
   Lemma B, and the cube-root pair gives `Re z0 = 1/2`.
3. (E4) `det(z0 I + A_Q(a)) = omega` (independent of a, nonzero), `z0 = 1/2 - i sqrt(3)/2`; so `chi(z0) != 0` (W2), and
   omega is a simple root of `Phi` (the cubic equals 2 at omega): (i).
4. (E5) `1/s'(z0) = 1/4 - i sqrt(3)/12` for every a (computed from (W2) as `-det(X)^2 / (mu^T adj(X)^2 1)`, `X = z0 I + A_Q`):
   `R* = 1/4 > 1/8`: (iii). The float check `c (nu_A(c) - 1/2) -> 0.2500` at c = 10^4 for a = 0, 1, 3 agrees (display only).
5. (E6) At beta = 2 the class vector is `y = (6,2,2,18,6,2,10)/(18a+17) > 0`, so by (W3) `w' > 0`: (iv).
6. Lemma B. ∎

CERTIFIED instances (`family_E_witnesses.py` -> `family_E_witnesses.json`; `ctln.analyze_support`, all at
`t = 20, beta = 821/400, eps = 1/21, delta = 421/8400, q = 421/420`, legal, x > 0, all Hurwitz minors positive; exact
check that `l^2+l+1` divides `n chi_perp` once and the cofactor is Hurwitz at shift -1/2):

| a | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| n | 7 | 8 | 9 | 10 | 11 | 12 |
| digraph6 | `&FNFL~[udo?` | `&GNbR^z~[zP[_` | `&HNpcv~^|~vNYDr?` | `&INwpL~z~v~n~[}sDr_` | `&JN{WR^~^~^~^~^~[~Y@[{?` | `&KN}KCv~z~|~}~~^~n~vNzODrw` |

### 4.4 Further n=8 finite-eps-only candidates (CERTIFIED, `scan_cube_n8.py`, `witness_cube.py`)

Among all 2,031,616 (= 124 x 2^14) one-node extensions of the 124 n=7 graphs (123 candidates + FNFL), floats flagged 9,808 with an
eigenvalue of `A|perp` near a primitive cube root and alpha near -1/2 (9,520 distinct). Exact triage: 798 have
`(l^2+l+1)^2`, 1,269 have another root with Re >= -1/2, 528 have the cube root as a FIXED eigenvalue (root of chi; never
finite-stable), 6,902 fail positivity at beta = 2 + 10^-4, 8 have `R* <= 1/8`, and **15 satisfy every hypothesis of
Lemma B**. All 15 have exact legal SFP witnesses (`cube_n8_witnesses.json`):
`&GM@RYJvJzp^w` (R*=1/6), `&GMPRYJvJzp^w` (1/2), `&GMPVYJvNzpNw` (1/3), `&GMPVZJvJzpNw` (1/3), `&GMRF^JtJzpNw` (1/3),
`&GNRf^jv\YoNw` (1/6), `&GNAb^rb\Yw^w` (1/2), `&GNQb^rb\Yw^w` (1/2), `&GNQf^rb\ZwNw` (1/4), `&GNQf^rr\YwNw` (1/4),
`&GNQf^rbXyo^w` (1/6), `&GNQf^zr\IsNw` (1/4, = E_1), `&GNQf^zbXYo^w` (1/6), `&GNRf^zr\ioNw` (1/6), `&GNRb^zrXIo^w` (1/6).
Witnesses at `(t, beta) = (20, 821/400)` (q = 421/420) or `(48, 647/320)` (q = 981/980). All have bidirected twin pairs.

Every one of the 16 boundary graphs lies on a 2- or 3-parameter integer "omega-locus" of closed-twin blow-up weights
(`omega_weights.py`): blowing up certain nodes arbitrarily keeps omega a moving root. Along these lines only the line of
Theorem E kept R* > 1/8 and positivity for all tested parameters (`omega_lines.py`, a = 1..5); one other line
(`&GMRF^JtJzpNw`, node 5) works for a = 1, 2, 3 and fails R* > 1/8 at a = 4 (R* = 1/12).

**Answer to question 2.** The cube roots of unity are not a symmetry of the graph. They come from an exact crossing of
a root locus in the twin-class-size parameter (`Phi_k(omega) = (k-1) omega^2` for the H3 blow-up). Theorem E is an
infinite family in which the crossing persists for every class size (`Phi` factors as `z(z^2+z+1)(...)` identically in a).
The finite-eps mechanism is Lemma B: a simple moving boundary pair whose first-order drift `R*` beats the legality curve's
drift 1/8.

---

## 5. Question 3: how large can q be?

`q_n*` := sup of q over legal points at which some non-clique n-node digraph has a stable positive full-support fixed
point. In the limit sense, `q_n^lim = sup (min(beta_c, sup Pos) - 1)`; finite points have `q = (beta-1) t/(t+1) < beta - 1`.

### 5.1 Certified lower bounds (CERTIFIED; search was HEURISTIC)

Search: one-node extensions of the 124 n=7 graphs (`run_extend_n8.py`, 2.0M graphs), then a beam search
(`beam.py`; nauty `labelg` for dedup) n = 8 -> 9 -> 10 -> 11. Every reported value is an exact point checked by
`ctln.analyze_support` and by the exact limit test (`certify_q.py`, `certify_top.py`):

| n | digraph6 | beta | t | eps | delta | q (exact) | q |
|---|---|---|---|---|---|---|---|
| 6 | Geneson (repo) | | | | | limit sup | 1.1664 |
| 7 | n=7 candidates (repo) | | | | | limit sup (plastic) | 1.3247 |
| 8 | `&G^ZnRba\ba^w` | 2709/1000 | 1000 | 1/1001 | 1709/1001000 | 1709/1001 | 1.7073 |
| 9 | `&H^\us{\HrEa^|~_` | 27691/10000 | 1000 | 1/1001 | 17691/10010000 | 17691/10010 | 1.7673 |
| 10 | `&I^Qz|NbkfuNa^}~|tw` | 14209/5000 | 10000 | 1/10001 | 9209/50005000 | 18418/10001 | 1.8416 |
| 11 | `&J^R\~R{\a^k^a^~^~]n^~?` | 14471/5000 | 100000 | 1/100001 | 861/45455000 | 17220/9091 | 1.8942 |

The n=8 graph is `&F^V\K[epa?` (n=7 candidate, q-sup about 1.17) plus one node that sends to all and receives from all
but one node. Its own class family `B7[S,all;k]` has `chi_L = z^3 (z^2+1)(z^2+z+2)` and positivity up to the root 2.7693 of
`b^3 - 3b^2 + b - 1`. The n=9..11 graphs add more such nodes and one gadget node that lifts the positivity bound to 3.023
(root of `b^5 - 5b^4 + 8b^3 - 7b^2 + 3b - 1`). The n=11 graph doubles a sink of the n=10 graph.

### 5.2 Theorem Q: q bounded away from 1 for all large n (PROVED; instance CERTIFIED)

`B = &F^S{^ReL|?` (an n=7 candidate), `S = (1,1,1,1,0,1,1)`, `Q_k = B[S,all;k]`, n = 7 + k.
- `z psi_D = z^3 (7z^4+12z^3+19z^2+11z+4)`, `chi_L = z^3 (z^4+z^3+3z^2+z+1)`.
- **Positivity at every beta > 2**: `det Mb = 3b^5 - 9b^4 + 16b^3 - 12b^2 + 7b - 1`, the components of `adj(Mb)1` are
  `b^2-b+1` (three times), `b^4-3b^3+5b^2-3b+1`, `(b-1)^2(b^2-b+1)`, `b^4-2b^3+3b^2-b+1`, `b^4-3b^3+4b^2-2b+1`, and the
  `sigma_S - 1` numerator is `b^4-b^3+3b^2-b+1`. None has a root in `[2, inf)` and all are positive at 3 (exact counts).
  So by (X2) every `Q_k` is positive at every beta > 2.
- **Stability**: for `beta* = 5/2, 27/10, 14/5`, every Hurwitz minor of `-Phi_k(-1/beta* - x)` is a polynomial in k with
  positive leading coefficient and no real root `>= k0 = 25, 79, 286` (exact).

**Theorem Q.** For every n >= 32 (resp. 86, 293) there is a non-clique n-node digraph whose full support is an SFP at
legal points with q arbitrarily close to 3/2 (resp. 17/10, 9/5). As k -> inf, `beta_c(k) -> beta_inf = -1/rho = 2.84416...`
(rho = the smaller real part of the roots of `z^4+z^3+3z^2+z+1`). So `liminf_n q_n^lim >= 1.84416...` along this family.
CERTIFIED instance: k = 25 (n = 32), beta = 5/2, t = 10^4, eps = 1/10001, delta = 3/20002, q = 15000/10001, legal,
x > 0, all Hurwitz minors positive (`family_Q_certificate.json`, Q5).

This **refutes** the tempting picture "clique blow-ups push q back to 1": the conjecture SL below is false.

### 5.2b Theorem Q2: two sink classes push q past 2 (PROVED; instances CERTIFIED)

**Lemma X1' / X2' (several sink classes; PROVED).** Let classes `V_j` (sizes `k_j`, in-sets `S_j`, `b_j = 1 - S_j`) all send to
every other node and be mutually bidirected. Then `n chi_{D'|perp} = z^{p-1} [z psi_D + p chi_D + 1^T adj(zI - D) B]`, with
`p = sum k_j` and `B = sum k_j b_j`. *Proof:* the class rows of `D'` are zero, so
`1^T (zI - D')^{-1} 1 = 1^T R 1 + sum_j k_j (1 + 1^T R b_j)/z` and `det(zI - D') = z^p chi_D`. Positivity: let
`s = sum_j k_j y_j`. Then `w_G = (1 - s) w_0` and `y_j = (1 - s)(sigma_j - 1)/(1 - gamma)`, so `s = kappa/(1 + kappa)` with
`kappa = sum k_j (sigma_j - 1)/(1 - gamma)`. Hence positivity holds iff `w_0 > 0` and `sigma_{S_j} > 1` for every j; this is
again k-free. The limit `p -> inf` with fixed ratios is `spec(D - bbar 1^T)`, `bbar = B/p`. CHECKED on explicit graphs
(`multiclass.check`).

`B' = &F^U|m^QTq?`, classes `V1` (size k, in-set all but node 2) and `V2` (size 3k, in-set all but node 3); n = 7 + 4k.
- Limit polynomial `det(zI - D + bbar 1^T) = z^3 (2z^2+z+2)^2 / 4` (bbar = (e_2 + 3 e_3)/4), min Re = -1/4, so the
  ceiling is beta = 4 (q -> 3).
- `B'` is positive at every beta > 2, and `sigma_{S_j} > 1` at every beta > 2 for both in-sets used. Exact: no root of
  any numerator in `[2, inf)`, and all are positive at 3. (An exploratory exact check shows the same for all seven
  one-node-deleted in-sets.)
- Exact thresholds (all Hurwitz minors in k positive for k >= k0): beta* = 5/2 -> k0 = 1; 3 -> 3; 7/2 -> 10; 19/5 -> 33.
- CERTIFIED: n = 19 (k = 3), beta = 3, t = 10^4, q = 20000/10001; n = 47 (k = 10), beta = 7/2, t = 10^4,
  q = 25000/10001. Both legal, x > 0, all Hurwitz minors positive (`family_Q2_certificate.json`).

### 5.2c Theorem Y: an 8-node core reaches q -> 3.96 (PROVED)

`C8 = &G^ZZZr|KjqTo` was found by a HEURISTIC core-design search (`core_design.py`, `core_grow.py`). Classes: `V1` (size 4k,
in-set all but node 0) and `V2` (size k, in-set all but node 2); n = 8 + 5k.
- Limit polynomial `z^2 (5z^6+5z^5+12z^4+7z^3+7z^2+2z+1)/5`, min Re about -0.2016 (display), ceiling beta about 4.96.
- Positive at every beta > 2. Exact: every positivity numerator has the sign of det on `[2, inf)`.
- Exact thresholds: beta* = 3 -> k0 = 13; 7/2 -> 24; 4 -> 46; 9/2 -> 117; 24/5 -> 363.
- So for every n = 8 + 5k >= 1823 there is a full-support counterexample with q arbitrarily close to 3.8.
- CERTIFIED instance (Y5): k = 13, n = 73, beta = 3, t = 10^4, eps = 1/10001, delta = 2/10001, q = 20000/10001;
  legal, x > 0, all Hurwitz minors positive.

The family ceilings grow with the core: the n=6 bases reach plastic + 1 (about 2.32), the 7-node cores 2.84 (one class)
and 4 (two classes), the 8-node core about 4.96 (two classes). HEURISTIC core growth to 9..12 nodes (`core_grow.log`) with a
weak inner optimiser (pairs of one-node classes only) did not beat 4.96.

### 5.3 Scale law (PROVED)

For real weights (Lemma W with `mu in R_{>0}^m`), `T_mu(z) = 1^T (z M^{-1} - D)^{-1} 1` satisfies `T_{N mu}(z) = T_mu(z/N)`,
and `u_{N mu}(r) := (D^T N M + rI)^{-1} 1 = u_mu(r/N)/N`. Hence `(mu, beta)` is limit-SFP iff `(N mu, beta/N)` is.
Multiplying all class sizes by N divides the attainable beta by N. A HEURISTIC integer-weight search on the n=8 record
graph (`weight_search.py`) found no blow-up better than doubling one sink (q 1.7051 -> 1.7111). Growth in q has to come
from new non-twin structure, or from a few huge classes as in Theorem Q, not from uniform blow-ups.

### 5.4 Conjectures about limits of clique classes

**Conjecture SL (REFUTED).** "For a sink class (T = all) with positivity, `spec(D - b1^T)` always has an eigenvalue with
Re <= -1/2". It passed 4,409 random positive (G,S,beta) tests (`conj_sink_limit.py`) and 6,330 random fractional mixtures
(`conj_sink_mix.py`). It is nevertheless false. Each proved family has every moving root at Re > -1/2 for every k, so
its limit polynomial has every root at Re >= -1/2 (root continuity). For 55 of the 135 distinct limit polynomials
(n6 + n7 catalogues), an exact gcd test (`Re` and `Im` parts of `chi_L(-1/2 + iy)`) shows there is no root on
`Re z = -1/2`, so min Re > -1/2 strictly. Example: the 12 families with limit polynomial `z^3(z^4+z^3+3z^2+z+1)` have
min Re -0.3516 (Theorem Q). Lesson recorded: random sampling missed the structured counterexamples.

**Conjecture P (OPEN).** If `D` is a nonzero 0/1 matrix with zero diagonal and `(D^T + rI)^{-1} 1 > 0` for some
`r in (0, 1/2)`, then D has an eigenvalue with `Re <= -1/2`.
- Evidence (HEURISTIC): 38,452 random positive cases (m = 3..12) and an adversarial single-flip ascent (150 restarts) never
  exceed min Re = -1/2. The supremum -1/2 is attained, for example by graphs with cube-root eigenvalues.
- Positivity is essential. Without it, non-nilpotent 0/1 matrices reach min Re = -0.452 (`conj_cycle.py`). Even
  `min in-degree >= 1` and `min out-degree >= 1` together reach -0.412 (`conj_degree.py`).
- PROVED special case: D normal, for example any abelian Cayley digraph. Then Re spec(D) = spec((D+D^T)/2), which
  contains a 2x2 principal block `[[0,w],[w,0]]` with `w >= 1/2`, so by interlacing `lambda_min <= -1/2`.
- Consequences if true: (a) joins `G v K_k` are limit-unstable for large k, since by (X3) with b = c = 0 the limit is
  `spec(D)`; (b) `G v G` is never limit-stable, because each right eigenvector x of D gives the eigenvector (x, -x),
  orthogonal to 1, of the doubled non-edge matrix. (b) is PROVED unconditionally whenever spec(D) has Re <= -1/2. For
  Geneson's G1, `spec(D)` contains `-0.5545 +- 1.1562 i`.

### 5.5 Is q_n* bounded? (OPEN)

- Upper bound: Geneson's Theorem 4.1 gives `q < n - 2 - (n-3) eps/2`. The trace bound gives
  `beta_c <= n(n-1)/e`, where e is the number of ordered non-edges (Lemma I.5).
- Lower bounds (CERTIFIED): 1.1664, 1.3247, 1.7073, 1.7673, 1.8416, 1.8942 for n = 6..11 (the n = 6, 7 values are limit
  suprema from the repo's exhaustive maps; the n >= 8 values are exact finite points).
- Theorem Q gives `q_n^lim >= 3/2` for all n >= 32 and `>= 1.8` for n >= 293, with family limit 1.8442.
- Theorem Q2 / Y: q arbitrarily close to 2 at n = 19, 5/2 at n = 47, 14/5 for n = 7 + 4k >= 139, and 3.8 for
  n = 8 + 5k >= 1823. Family ceilings are 4 (7-node core) and about 4.96 (8-node core).
- Ceiling of a fixed core (PROVED, trace): take sink classes that are not joins (each `b_j != 0`). The limit matrix
  `D - bbar 1^T` has trace `-|bbar|_1 <= -1`, so its min real part is `<= -1/m`. Hence the many-class limit of an m-node
  core has `beta <= m`, i.e. `q <= m - 1`. The 7-node core of Q2 reaches 3 of a possible 6; the 8-node core of Y reaches
  3.96 of 7. So q_n* -> inf would need cores of unbounded size whose ceilings grow. (Join classes, `b = 0`, escape this
  bound; for them Conjecture P is the relevant obstruction.)
- Mechanism seen in every record graph: a 7-node core whose positivity range sets the ceiling, plus sink nodes and sink
  classes that move the stability boundary up to that ceiling. Growth so far comes from cores with larger positivity
  ranges (2.3247 -> 2.7693 -> 3.023 -> all beta for `&F^S{^ReL|?`). I found neither a proof of boundedness nor a family
  with `q -> inf`. The evidence points towards unboundedness: the ceiling grows with the core size (2.32, 2.84/4, 4.96 for
  6, 7, 8-node cores) and nothing caps it below `n - 2`. The quantity to study is
  `sup { -1/min Re spec(D - bbar 1^T) : bbar in conv(light sets of a core D), D positive }`.
  Conjecture P controls only the case bbar = 0.

---

## 6. Files (all in THEORY/; every exact claim re-runs)

Exact certificate files:
- `family_F_certificate.json` (from `family_F_proof.py`): Theorem F data and witnesses.
- `family_E_certificate_part1.json` (from `family_E.py`): Theorem E closed forms. `family_E_witnesses.json` (from
  `family_E_witnesses.py`): instances a = 0..5.
- `family_Q_certificate.json` (from `family_Q_proof.py --instances`): Theorem Q data and the n = 32 witness.
- `families_catalogue_n6.json`, `families_catalogue_n7.json` (from `run_family_prover.py [n7]`): 9 + 617 proved families.
- `cube_n8_good.json`, `cube_n8_witnesses.json` (from `scan_cube_n8.py`, `witness_cube.py`): the 15 n=8 boundary
  candidates.
- `top_n8.json`, `cert_q_910.json`, `cert_q_11.json` (from `certify_top.py`, `certify_q.py`): high-q instances.
- `family_Q2_certificate.json` (from `family_Q2_proof.py --instances`): Theorem Q2 data, witnesses n = 19, 47.
- `family_Y_certificate.json` (from `family_Y_proof.py --instances`): Theorem Y data (and the n = 73 witness).
- `F_family_witnesses.json` (from `cert_F_instances.py`): early F witnesses.
- `verify_all.py` re-checks every instance record above with `ctln.analyze_support` and the exact limit/boundary tests.
  Last run: **ALL CHECKS PASSED (78 records)**, including the n = 19, 32, 47, 73 family instances.

Libraries: `tlib.py` (exact helpers), `classext.py` (Lemma X1), `multiclass.py` (Lemma X1'), `weighted.py` (Lemma W),
`posfam.py` (positivity polynomials), `rstar.py` (Lemma AS coefficient), `boundary_tools.py`, `family_prover.py`, `omega_weights.py`.
Heuristic search scripts (floats steer only): `heur_extend.py`, `run_extend_n8.py`, `beam.py`, `localsearch.py`,
`weight_search.py`, `heur_blowclass.py`, `mixture_limit.py`, `simplex_opt.py`, `core_design.py`, `core_grow.py`,
`conj_*.py`, `explore_*.py`. Their JSON/log outputs are
search records, not certificates.

---

## 7. Items for ANOMALIES.md (not appended here: the repository was not to be modified)

1. `&FNFL~[udo?` (the n=7 finite-eps-only candidate) is the signed-off n=6 counterexample idx 1387528 with one node
   doubled into a closed-twin pair. Its alpha = -1/2 is an exact crossing of the class-size root locus at the cube root
   of unity (`Phi_k(omega) = (k-1) omega^2`). Theorem E then makes this persistent for every n.
2. The plastic number shows up as a **positivity** boundary: `b^3 - 3b^2 + 2b - 1 = 0` (b - 1 = plastic) is the root of a
   single adjugate component (`adj_4`) in Geneson's graph, and it caps `F_k` for k = 2..8. The n=7 summary's six
   "plastic" q_c values have the same origin.
3. Conjecture SL survived 10,739 random float tests and is false. Its counterexamples are structured: 55 limit
   polynomials among the exactly proved families. Random sampling is not evidence for spectral conjectures here.
4. Single-node extensions of the n=7 candidates jump from q 1.3247 to 1.7073 at n = 8. The jump comes from one
   near-universal node added to a candidate whose own q-sup is only about 1.17. Large-q graphs are therefore not near
   the large-q graphs of the previous size.
5. Every proved clique-class family's limit polynomial `chi_L` has a factor `z^j` with j >= 1. The class nodes and the
   activated base nodes contribute the eigenvalue 0 in the D-picture.
