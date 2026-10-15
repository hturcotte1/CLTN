# One parameter point, infinitely many graphs (workflow wf7, label `prover_fixed`)

Date 2026-10-15. Conventions are exactly CLAUDE.md §2: `A[i][j] = 1` iff `i -> j`; `D := J - I - A` (non-edge digraph,
`D[i][j] = 1` iff `i != j` and `i -/-> j`; the notes' `N = J - I - A^T = D^T`); `t = (1-eps)/eps`, `beta = 1 + delta/eps`,
`r := 1/beta`, `a := t/beta`; legal iff `t(beta-2) > 1`. Nothing in the repository was modified. Every script below lives in
`/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf7/prover_fixed/` (below: `PF/`) and
re-runs standalone with `python3 <script>` from that directory (it imports the repo's `ctln.py` read-only).

Labels: **PROVED** (proof written here; computer algebra steps are exact and named), **CERTIFIED** (exact rational/integer
computation on explicit graphs, script and output file named), **HEURISTIC** (floating point, used only to choose what to
try; never a decision), **OPEN**, **REFUTED**.

> **Hand-verification status (CLAUDE.md §3).** Every graph on n >= 8 nodes named below as a stable non-clique full support is
> a **CANDIDATE counterexample pending human hand-verification**. The cores used (`&F^U|m^QTq?`, `&F^S{^ReL|?`) are among the
> 124 signed-off n=7 graphs; the 8-node core `&G^ZZZr|KjqTo` is not signed off. The proofs use the cores only as explicit
> matrices. Nothing beyond the requested certificates was built on the new graphs.

---

## 0. Answer and summary

**Question.** Is there ONE legal `(eps, delta)` and an infinite sequence of non-clique digraphs `G_k`, `|G_k| -> inf`, whose full
support is a stable positive fixed point at that same `(eps, delta)`?

**Answer: yes (PROVED), and robustly.**

- At `(eps, delta) = (1/21, 11/210)` (`beta = 21/10`, `t = 20`, `q = 22/21`), which is exactly the point where Geneson's
  clique-class family `F_k` dies after `k = 40`, the family `R_k` is an SFP for **every** `k >= 0`. Here `n = 7 + k`,
  `R = &F^U|m^QTq?` (signed off), and the k new nodes form a bidirected clique that sends to every node of R and receives from
  every node except node 2.
- The same holds at **every** legal `(eps, delta)` with `delta <= (23/20) eps` (i.e. `beta <= 43/20`), uniformly in `k >= 0`.
- Every legal `(eps, delta)` with `beta = 1 + delta/eps <= 24/5` carries an infinite family. This includes `eps = 1/2`.

| # | statement | status | where |
|---|---|---|---|
| S | K positive stable iff every eigenvalue of `D + aJ` has `Re > -1/beta` | PROVED | §1 |
| P | sink-class families: `det(zI - aJ - D') = z^(p-1) Psi`, `Psi = z det(zI - D - aJ) - a sum_j k_j det(zI - D + b_j 1^T)` | PROVED; CHECKED on 24 explicit families | §1, `check_lemmaP.py` |
| X2' | positivity of every member iff `w_G(beta) > 0` and `sigma_{S_j}(beta) > 1` (k-free, t-free) | PROVED (re-derived); CHECKED on 140 cases (61 positive) | §1 |
| FP | **t drops out.** For fixed legal `(beta, t)`, as `k -> inf` the spectrum converges to `spec(Lbar) u {0} u {+inf}`, `Lbar = D - bbar 1^T`, independent of t. Eventually SFP iff positivity at beta and `spec(Lbar)` in `Re > -1/beta`; eventually unstable if some root has `Re < -1/beta`. (iv): one k0 works uniformly for all legal points with `beta <= beta_hi < beta_inf` | PROVED | §2 |
| A | `R_k` is an SFP at `(1/21, 11/210)` for every `k >= 0`; all five Hurwitz minors are integer-coefficient polynomials in k with all coefficients positive | PROVED (exact) | §3, `theorem_A.py` |
| A+ | `R_k` is an SFP at every legal `(eps, delta)` with `beta <= 43/20`, for every `k >= 0` | PROVED (exact sign-pattern certificate in (a, k)) | §4, `theorem_B.py` |
| B | region theorems for R, Q, Q2, Y: SFP at every legal point with `beta <= beta1`, for all `k >= k0(beta1)`, up to `beta1 = 24/5` (Y, k0 = 1764) | PROVED (exact) | §4 |
| B' | every legal `(eps, delta)` with `beta < beta_Y = 4.9607...` carries an infinite family (Y). With rational `beta <= 24/5` an explicit k0 is given | PROVED | §4 |
| T | question (a): `T_k(beta)` is bounded for R, Q (it equals the legal boundary for `k >= k0`). For F it is infinite for `k >= 39` at beta = 21/10 | PROVED/CERTIFIED | §5, `Tk_table.py` |
| N1 | F_k at `(1/21, 11/210)`: stable exactly for `k = 0..40`, unstable for every `k >= 41`. F fails at **every** fixed legal point for large k (cube roots in `chi_L`) | PROVED (exact) | §6, `negative.py` |
| U | uniform blow-ups `H^(N mu)` of any non-clique H fail at every fixed point for large N. For blow-ups of a fixed core, the fastest-growing classes must induce a complete digraph | PROVED | §6 |
| P-refuted | Conjecture P of notes_theory §5.4 is **false**: `D` of `&F^U|m^QTq?` has every eigenvalue with `Re > -10/21` and `(D^T + (10/21) I)^{-1} 1 > 0`. The join `R v K_k` is an SFP at `(1/21, 11/210)` for every k | REFUTED (exact) | §7, `conjP_refutation.py` |
| n-set | at `(1/21, 11/210)` a non-clique n-node digraph with a stable full support exists **iff n >= 6** | PROVED (n <= 5 from the project's exact n<=5 map) | §3 |
| O | Does every legal `(eps, delta)` carry an infinite family? Open for `beta >= beta_Y`; for sink/join-class families this reduces to a t-free question about cores | OPEN | §8 |

---

## 1. Three identities (PROVED)

**Lemma S (stability in non-edge coordinates).** For a digraph G' on n nodes with non-edge matrix D' and any `eps, delta > 0`,
condition (3) for the full support holds iff every eigenvalue `mu` of `D' + aJ` satisfies `Re mu > -r`.

*Proof.* By (S1)/(S2) of notes/lemma_I_spectral.md, `I - W = eps K` with `K = I + beta N + tJ`, `N = J - I - A^T = D'^T`, and (3)
iff K is positive stable. `K^T = I + beta D' + tJ = I + beta (D' + aJ)` has the same spectrum as K, namely
`{1 + beta mu : mu in spec(D' + aJ)}`, and `Re(1 + beta mu) > 0` iff `Re mu > -1/beta`. ∎

This is identity (E) of notes/lemma_M_refuted.md read in the D-picture: `nu_A(c) = 1 + min Re spec(D + aJ)`, `c = 1 + a`.

**Sink-class families.** Take a core G (m nodes, non-edge matrix D) and classes `V_1..V_s` of sizes `k_1..k_s`, `p = sum k_j`.
Each class is a bidirected clique, the classes are mutually bidirected, and every class node sends to every other node. Node i of
G sends to `V_j` iff `S_j[i] = 1`; put `b_j := 1 - S_j`. (s = 1 is `G[S, all; k]` of notes_theory §2.2; `S = all` is the join
`G v K_k`.)

**Lemma P.** With `chi_a(z) := det(zI - D - aJ)` and `chi_{L_j}(z) := det(zI - D + b_j 1^T)`:

    det(zI_n - aJ_n - D') = z^(p-1) Psi(z),    Psi(z) := z chi_a(z) - a sum_j k_j chi_{L_j}(z).

Psi is monic of degree m+1, and its coefficients are affine in `(a, a k_1, ..., a k_s)`. Since `b -> det(zI - D + b 1^T)` is affine
(matrix determinant lemma), `sum_j k_j chi_{L_j} = p chi_{Lbar}` with `Lbar := D - bbar 1^T`, `bbar := sum_j k_j b_j / p`.

*Proof.* The class rows of D' are zero (class nodes send everywhere). Ordering core first, `D' = [[D, B],[0, 0]]` with
`B = [b_1 1^T_{k_1}, ..., b_s 1^T_{k_s}]`, so `det(zI - D') = z^p chi_D(z)`, and with `R = (zI - D)^{-1}`,
`(zI - D')^{-1} = [[R, R B / z],[0, I/z]]`. Hence `1^T (zI - D')^{-1} 1 = 1^T R 1 + (1^T R B 1 + p)/z`, with `B 1 = sum_j k_j b_j`.
By the matrix determinant lemma, `det(zI - D' - aJ) = det(zI - D') (1 - a 1^T (zI - D')^{-1} 1)`, so it equals
`z^(p-1) [ z chi_D - a z 1^T adj(zI-D) 1 - a sum_j k_j (chi_D + 1^T adj(zI-D) b_j) ]`. Now
`chi_D - a 1^T adj(zI-D) 1 = chi_a` (determinant lemma again), and `chi_D + 1^T adj(zI-D) b_j = det(zI - D + b_j 1^T)`. Both
sides are polynomials that agree off a finite set, so they agree everywhere. ∎

CHECKED exactly (`check_lemmaP.py` -> `check_lemmaP.json`, part C1) on the explicit graphs of R, Q, F (k = 1, 2, 3), Q2
(k = 1, 2), Y (k = 1) and 12 random (core, in-sets, sizes), at `a in {3/7, 5, 200/21}`. Both sides are monic of degree n in z
and agree at n + 2 exact rational points, which proves the identity for each tested graph and a. A negative control (wrong k
in Psi) fails as it should. `theorem_A.py` also checks R for k = 1..4 symbolically in z. Part C2 checks Lemma S against
`ctln.analyze_support` on 40 random (graph, beta, t) (4 stable, 36 unstable; all agree).

**Lemma X2' (positivity is k-free and t-free).** For rational `beta > 1` and all sizes `k_j >= 1`: the full support of the
family member is positive iff `w_G(beta) := (A^T + gamma I)^{-1} 1 > 0` and `sigma_j := sum_{i: S_j[i]=1} w_G,i > 1` for every j.

*Proof.* (This is notes_theory X2/X2'; it is re-derived here because everything rests on it.) By Lemma I.4(c, g) of
lemma_I_spectral.md, positivity at any t > 0 is equivalent to `w' := (A'^T + gamma I)^{-1} 1 > 0`, and `A'^T + gamma I` is
nonsingular at rational beta. By uniqueness and symmetry, w' is constant (`y_j`) on `V_j`.
- Core rows: every class node sends to every core node, so `(A^T + gamma I) w'_G + s 1 = 1` with `s = sum k_j y_j`. Hence
  `w'_G = (1 - s) w_G`.
- Class-`V_j` rows: `(gamma - 1) y_j + s + (1 - s) sigma_j = 1`. Hence `y_j = (1 - s)(sigma_j - 1)/(1 - gamma)`.
- Summing gives `s = (1 - s) kappa` with `kappa = sum_j k_j (sigma_j - 1)/(1 - gamma)`. So `1 + kappa != 0` and
  `1 - s = 1/(1 + kappa)`.

If `w_G > 0` and every `sigma_j > 1`, then `kappa > 0` and all entries are positive. Conversely, `w' > 0` forces `w_G` to be
one-signed, hence `> 0` (Lemma I.4(P3')). That forces `1 + kappa > 0`, and then every `sigma_j > 1`. ∎

CHECKED against ctln.
- `check_lemmaP.json` part C3: 60 random cases, all non-positive and all agreeing.
- `check_X2.py` -> `check_X2.json`: 80 draws from the catalogue families and Q2, Y, with `beta in (2, 3]` and `k in [1, 4]`.
  61 are positive and 19 are not, and all agree.

*Remark (irrational beta).* The "if" direction needs only `det(A^T + gamma I) != 0` (true wherever the exact positivity
polynomials show `w_G > 0`), not rationality. The vector `w'` built above satisfies `(A'^T + gamma I) w' = 1` by direct
substitution. As in the proof of I.4(c), `C w' >= w' > 0` gives `W1' > 1`, and then `K w' = Den 1` with `Den > 0`. So
`x = theta (1+t) w'/Den > 0` solves the fixed-point equation, and it is the unique solution whenever K is nonsingular, which
stability guarantees.

---

## 2. Theorem FP: at a fixed parameter point, t drops out of the large-family limit (PROVED)

Fix a core with in-sets `S_j` and proportions `rho_j` (positive integers, `k_j = rho_j k`), and fix legal `(beta, t)`, i.e.
`r = 1/beta`, `a = t/beta`. Put `bbar = sum rho_j b_j / sum rho_j`, `Lbar = D - bbar 1^T`.

**Theorem FP.**
- (i) If `w_G(beta) > 0`, `sigma_j(beta) > 1` for all j, and every root of `chi_{Lbar}` has `Re > -r`, then there is `k0` such
  that for every `k >= k0` the full support of `G_k` is a stable positive fixed point at `(eps, delta) = (1/(1+t), (beta-1)/(1+t))`.
- (ii) If `chi_{Lbar}` has a root with `Re < -r`, then there is `k1` such that `G_k` is unstable at `(eps, delta)` for every
  `k >= k1`.
- (iii) If `w_G(beta) > 0` fails, or some `sigma_j(beta) <= 1`, then no member with `k >= 1` is positive.
- (iv) **Uniform version.** Let `beta_hi < beta_inf` and assume positivity on `(2, beta_hi]`. Then there is one `k0` such that
  for every legal `(eps, delta)` with `beta <= beta_hi` and every `k >= k0`, the full support of `G_k` is an SFP. In particular
  `T_k(beta) = 1/(beta - 2)` (the legal boundary) for all `k >= k0(beta)`, which answers question (a) in general. Theorem B
  makes k0 explicit for rational beta_hi.

So the set of fixed legal points that carry the family for all large k is, up to the boundary `min Re spec(Lbar) = -1/beta`,
`{legal (beta, t) : beta in Pos(G, S) and beta < beta_inf := -1/min Re spec(Lbar)}`. It is a union of full t-strips.

*Proof.*
1. By Lemmas S and P, `G_k` is stable iff every root of `Psi_k = z chi_a - a P chi_{Lbar}` (`P = k sum rho_j`) has `Re > -r`.
   The extra roots `z = 0` of `z^(p-1)` are harmless because `0 > -r`.
2. Write `Psi_k = -aP (chi_{Lbar} - z chi_a/(aP))`. `chi_{Lbar}` is monic of degree m, and `z chi_a` is a fixed polynomial of
   degree m+1.
3. (i) Let `zeta_1..zeta_m` be the roots of `chi_{Lbar}`. Choose disjoint closed disks `B_i` around the distinct `zeta_i` inside
   `{Re > -r}`, and a circle `|z| = R0` enclosing them. On each `dB_i` and on `|z| = R0` we have `|chi_{Lbar}| >= c > 0` and
   `|z chi_a| <= C`. For `aP > C/c`, Rouche's theorem gives Psi_k exactly as many roots in each `B_i` as `chi_{Lbar}` has
   (m in total), and m roots in `|z| < R0`.
4. The remaining root lies outside `|z| < R0`. It is real, because a non-real root would bring its conjugate, also outside.
   The coefficient of `z^m` in Psi_k is `-a(m + P)`, so the sum of the roots is `a(m + P)`. Hence the remaining root is
   `>= a(m + P) - m R0 > 0` for P large. So all roots have `Re > -r`.
5. Positivity is Lemma X2'.
6. (ii) A disk around a root with `Re < -r`, inside `{Re < -r}`, contains a root of Psi_k for all large P (Rouche again).
7. (iii) Lemma X2' (converse direction).
8. (iv) Write `Psi_k = z chi_D - a (z psi_D + P chi_{Lbar})` (`chi_a = chi_D - a psi_D`, `psi_D = 1^T adj(zI - D) 1`), so
   `Psi_k/(aP) = -chi_{Lbar} + [z chi_D/(aP) - z psi_D/P]`.
   - Legal points with `beta <= beta_hi` have `r >= r_hi := 1/beta_hi` and `a > a_hi := r_hi^2/(1 - 2 r_hi) > 0`.
   - On the fixed circles of step 3, drawn inside `{Re > -r_hi}` (which lies inside `{Re > -r}`), the bracket is at most
     `C (1/a_hi + 1)/P`, uniformly in `a >= a_hi`. So Rouche applies with one P0 for all such a.
   - The escaping root is `>= a_hi (m + P) - m R0 > 0`.
   - Positivity is k-free and t-free (Lemma X2' and its remark). ∎

Consequences (exact limit data in `ceilings.py` -> `ceilings.json`):

| family | `chi_{Lbar}` | min Re (display) | `beta_inf` | positivity |
|---|---|---|---|---|
| F (Geneson + class) | `z^4 (z^2+z+1)` | **-1/2 exactly** | 2 | `(2, b_p)` |
| R, Q (one class) | `z^3 (z^4+z^3+3z^2+z+1)` | -0.351597 (`2 rho` is a root of `s^4+2s^3+6s^2+5s+1`) | 2.844165 | all `beta > 2` |
| Q2 (classes 1:3) | `z^3 (2z^2+z+2)^2 / 4` | -1/4 exactly | 4 | all `beta > 2` |
| Y (classes 4:1) | `z^2 (5z^6+5z^5+12z^4+7z^3+7z^2+2z+1)/5` | -0.201586 (degree-15 minimal polynomial for `2 rho` in ceilings.json) | 4.960659 | all `beta > 2` |

For each of R, Q, Q2, Y, ceilings.py checks exactly that `chi_{Lbar}(-1/beta1 - x)` is Hurwitz at every beta1 used in §4. It
also checks that it is **not** Hurwitz at a rational just above `beta_inf` (1423/500, 4001/1000, 2481/500).

**F is impossible at every fixed point.** `z^2 + z + 1` divides `chi_L`, so `chi_L` has roots with `Re = -1/2 < -1/beta` for
every legal beta (> 2). By FP(ii), F_k is unstable for all large k at every fixed legal `(eps, delta)`. This is the structural
reason behind the KNOWN data, and it explains why the F family needs a k-dependent schedule (`beta -> 2`). Remark: R's quartic
is palindromic; in `w = z + 1/z` it reads `w^2 + w + 1`. So the cube roots of unity govern both limits. For F it is `z = omega`
(Re -1/2, fatal). For R it is `z + 1/z = omega` (Re -0.3516 and -0.1484, harmless below `beta = 2.844`).

---

## 3. Theorem A: one point, every k (PROVED; instances CERTIFIED)

`R = &F^U|m^QTq?` (signed-off n=7 counterexample), rows

    0111110 / 1011011 / 1101101 / 1100111 / 1101001 / 0010101 / 1100100     (row i lists the j with i -> j)

(12 one-way pairs, no non-adjacent pair; rows as stored in `theorem_A_certificate.json`), and `S = (1,1,0,1,1,1,1)`. `R_k := R[S, all; k]`, `n = 7 + k`.

**Theorem A.** At `(eps, delta) = (1/21, 11/210)` (`t = 20`, `beta = 21/10`, `q = 22/21`, legal), for every integer `k >= 0` the
full support of `R_k` is a stable positive fixed point. No `R_k` is a clique.

*Proof (all steps exact, `theorem_A.py` -> `theorem_A_certificate.json`).*
1. Lemma P: `Psi_k = z^3 g_k`, with
   `g_k(z) = z^5 - 3z^2 - 2z - 1 - a[(k+7) z^4 + (k+12) z^3 + (3k+21) z^2 + (k+13) z + (k+5)]`.
   Here `z^3 chi_{L}` with `chi_L = z^3 (z^4+z^3+3z^2+z+1)`, and `z chi_a = z^3 (z^5 - 7a z^4 - ...)`.
2. At `a = 200/21`, `r = 10/21`, `h_k(x) := -g_k(-10/21 - x)` is monic in x, with coefficients
   `1, 50(4k+29)/21, 200(19k+33)/441, (258600k+1217983)/9261, 2(1489900k+1831049)/194481, 13(3433400k+8122137)/4084101`.
3. Its Hurwitz minors, cleared of positive denominators, are
   - `21 Delta_1 = 200k + 1450`
   - `9261 Delta_2 = 760000k^2 + 6571400k + 8352017`
   - `85766121 Delta_3 = 77344000000k^3 + 759190040000k^2 + 1860475888600k + 2626155959161`
   - `16679880978201 Delta_4 = 101566081600000000k^4 + 928678436864000000k^3 + 2091830212667400000k^2 + 5008829616746115200k + 4346419052022427717`
   - `68122318582951682301 Delta_5 = 4533320799350720000000000k^5 + ... + 458928742999171704870925977` (all six coefficients positive)

   Every coefficient is positive, so every `Delta_j(k) > 0` for every real `k >= 0`. Hence (Routh–Hurwitz) every root of `g_k`
   has `Re > -10/21`; the roots `z = 0` are harmless; by Lemma S, `R_k` is stable. (theorem_A.py also verifies the
   undivided degree-8 version: all eight minors have only positive coefficients in k.)
4. `w_R(21/10) = (2741781, 2741781, 841071, 1282071, 695100, 1767171, 695100)/7691294 > 0` and `sigma_S = 381654/295819 > 1`, so
   by Lemma X2' every `R_k` is positive.
5. R has one-way pairs (e.g. 0 -> 5, 5 -/-> 0), so it is not a clique, and neither is any `R_k` (R is an induced subgraph). ∎

CERTIFIED instances (`certify_members.py R_eps1_21` -> `cert_R_eps1_21.json`). ctln.analyze_support at `(1/21, 11/210)`:
legal, x > 0, all Hurwitz minors of the exact characteristic polynomial positive, not a clique. The Psi-verdict agrees.

| k | 0 | 1 | 2 | 5 | 10 | 20 | 100 |
|---|---|---|---|---|---|---|---|
| n | 7 | 8 | 9 | 12 | 17 | 27 | 107 |
| seconds | 5 | 0.5 | 0.5 | 0.5 | 0.8 | 2.4 | 755 |

`theorem_A.py` also cross-checks the Psi-verdict against ctln for k = 0..5.

**Corollary (complete size spectrum at this point).** At `(1/21, 11/210)`, a non-clique n-node digraph whose full support is
an SFP exists **iff n >= 6**.
- n = 6: Geneson's G1 = F_0 (CERTIFIED by ctln in `verify_all.py`; also stable by negative.py, minors positive at k = 0).
- n >= 7: `R_{n-7}`.
- n <= 5: none, by the project's exact n <= 5 classification. The only full-support SFP for n <= 5 at any legal point is K_n
  (notes/lemma_M_refuted.md §2.4; Stage 3 CAD: "n=5: 705 never + K5").

---

## 4. Region theorems: uniform in t and in k (PROVED; `theorem_B.py` -> `theorem_B_certificate_{R,Q_Q2,Y}.json`)

**Reduction (B3).** `spec(D' + aJ)` does not involve r. So if every eigenvalue has `Re > -r1`, the same holds for every
`r >= r1`.

Legality is `a > r^2/(1 - 2r)`, and the right side increases on `(0, 1/2)`. Hence every legal point with `beta <= beta1` has
`r >= r1 := 1/beta1` and `a > a1 := r1^2/(1 - 2 r1)`. So it suffices that every Hurwitz minor of `-(-1)^{m} g(-r1 - x)`, a
polynomial in `(a, k)`, is positive on `{a >= a1, k >= k0}`.

**Exact certificate.** After substituting `a = a1 + u`, `k = k0 + v`, every coefficient of every minor is `>= 0` and the
constant term is `> 0`. Positivity on `[2, beta1]` is checked separately: det, every entry of `adj(beta A^T + (beta-1) I) 1`,
and every `sigma_j - 1` numerator has no root in `[2, beta1]` and has the sign of det (exact Sturm counts).

**Theorem B.** For every legal `(eps, delta)` with `1 + delta/eps <= beta1`, and every integer `k >= k0`, the full support of
the family member is a stable positive fixed point:

| family | n | beta1 | k0 |
|---|---|---|---|
| R | 7 + k | 43/20 | **0** |
| R | 7 + k | 11/5 | 4 |
| R | 7 + k | 23/10 | 10 |
| R | 7 + k | 12/5 | 19 |
| R | 7 + k | 5/2 | 33 |
| R | 7 + k | 13/5 | 60 |
| R | 7 + k | 27/10 | 127 |
| R | 7 + k | 14/5 | 627 |
| Q | 7 + k | 21/10 | 4 |
| Q | 7 + k | 5/2 | 40 |
| Q2 | 7 + 4k | 5/2 | 2 |
| Q2 | 7 + 4k | 3 | 6 |
| Q2 | 7 + 4k | 7/2 | 19 |
| Q2 | 7 + 4k | 19/5 | 73 |
| Y | 8 + 5k | 3 | 22 |
| Y | 8 + 5k | 7/2 | 53 |
| Y | 8 + 5k | 4 | 133 |
| Y | 8 + 5k | 9/2 | 502 |
| Y | 8 + 5k | 24/5 | 1764 |

(The legal boundary at beta1 is `t1 = a1 beta1 = 1/(beta1 - 2)`. Where the first try failed, the k0 search history,
e.g. 501 -> 627, is recorded in the JSON.)

In particular **Theorem A+**: `R_k` is an SFP at every legal `(eps, delta)` with `delta <= (23/20) eps`, for every `k >= 0`.
This region is `{eps/(1-eps) < delta <= 23 eps/20}`, nonempty exactly for `eps < 3/23`.

**Corollary B'.** By FP(i) with Y, every legal `(eps, delta)` with `beta < beta_Y = 4.96066...` carries an infinite family of
non-clique SFPs. Positivity of Y on `[2, 5]` is exact (`positivity_ranges.py` -> `positivity_ranges.json`; same for R, Q on
`[2, 3]` and Q2 on `[2, 4]`, so Theorem FP's window is exactly `(2, beta_inf)` for each). FP(i) applies at irrational beta too:
see the remark after Lemma X2'. With `t = 1` this includes **eps = 1/2**, every `delta in (1, 1.98)`.
- CERTIFIED at `(1/2, 9/8)` (`beta = 13/4`, `t = 1`, `q = 9/8`): Q2 members `k = 9, 10, 12` (n = 43, 47, 55)
  (`cert_Q2_eps1_2.json`). The exact pointwise threshold there is `k0 = 9` (`k_threshold`; Q2 at `(1/2, 5/4)` has `k0 = 16`).
- `(1/2, 9/8)` has `beta = 13/4 <= 7/2`, so the Q2 region theorem gives every `k >= 19`. The exact pointwise threshold at
  this point is `k0 = 9`.
- Exact pointwise thresholds (`pointwise_thresholds.py` -> `pointwise_thresholds.json`):

  | family | `(eps, delta)` | every k >= | n >= |
  |---|---|---|---|
  | R | (1/21, 11/210) | 0 | 7 |
  | R | (1/8, 23/160) | 0 | 7 |
  | Q | (1/21, 11/210) | 3 | 10 |
  | Q2 | (1/2, 9/8) | 9 | 43 |
  | Q2 | (1/2, 5/4) | 16 | 71 |
  | Y | (1/2, 9/8) | 31 | 163 |
  | Y | (1/2, 3/2) | 90 | 458 |

  The Y instances are not instance-certified, because of cost (n >= 163).
- Largest eps covered: `eps < (beta_Y - 2)/(beta_Y - 1) = 0.7475` (display).

Randomized exact cross-check of the region theorems (tests, not proofs): `crosscheck_region.py` -> `crosscheck_region.json`.
All passed.
- X1: 25 random legal points with `beta <= 43/20`, random `k in [0, 30]`; ctln says SFP every time.
- X2: 10 random points for Q2 with `beta <= 5/2`, `k in [2, 6]`; ctln says SFP every time.
- X3: 30 random points with `beta in [11/5, 14/5]`. The Psi-verdict equals ctln's verdict at every point (14 stable,
  16 unstable).

---

## 5. Question (a): T_k(beta) (PROVED/CERTIFIED; `Tk_table.py` -> `Tk_table.json`)

At fixed beta and k, Psi_k is linear in a. So the stable t-set is computed exactly from the real roots (in a) of the Hurwitz
minors, with one exact rational sample per cell.
- **Q at beta = 5/2** (legal t > 2): unstable for every t when `k <= 24`. Then `T_25 = 111.50`, `T_26 = 23.31`, `T_30 = 5.65`,
  `T_35 = 2.92`, `T_39 = 2.117`. For every `k >= 40`, `T_k = 2`, the legal boundary (and Theorem B proves this for all k >= 40).
- **R at beta = 5/2**: `T_17 = 1047.9`, `T_20 = 10.22`, `T_25 = 3.875`, `T_30 = 2.396`, and `T_k = 2` for every `k >= 33`.
- **F at beta = 21/10** (legal t > 10): `T_k = 10` (every legal t is stable) for `k <= 38`. For `k = 39, 40, 41` the stable set
  is a **bounded** window: `(10, 87.98)`, `(10, 28.25)`, `(10, 16.83)` (finite-eps-only behaviour at this beta). For
  `k = 45, 60` it is empty.

So the answer to (a) is: bounded, indeed eventually equal to the legal boundary, for every family whose limit spectrum is
strictly inside `Re > -1/beta`. This is PROVED in general by Theorem FP(iv); the tables show it exactly. For F the sequence is unbounded (infinite) from k = 39 on. Decimals are display-only; the JSON
endpoints are roots of the stated minors.

---

## 6. Negative results (PROVED; `negative.py` -> `negative_certificate.json`)

**N1 (F at the KNOWN point).** At `(1/21, 11/210)`, the seven Hurwitz minors of F's `h_k` are positive at k = 0.
- Minors 1-4 have no root in `[0, inf)`.
- Minors 5, 6, 7 have negative leading coefficient. Their smallest roots `>= 0` are 55.99, 40.6079, 40.6079.
- Minor 6 is `< 0` on all of `[41, inf)` (exact Sturm count).

Hence **F_k is stable exactly for k = 0..40 and unstable for every k >= 41.** This extends the KNOWN "k <= 38 stable; 41, 45, 60
unstable". ctln confirms k = 40 (stable) and k = 41 (unstable, minor 6 negative) on the explicit 46- and 47-node graphs.

**Lemma U (uniform blow-ups fail).** Let H be a non-complete digraph on m >= 2 nodes, `mu in Z_{>0}^m`, and let `H^(N mu)` be
the closed-twin blow-up (node i becomes a bidirected clique of size `N mu_i`). Fix any `a, r > 0`. Then `H^(N mu)` is unstable
at `(beta, t) = (1/r, a/r)` for all large N. Moreover, stability holds iff `N < r/|Re lambda|` for the worst eigenvalue
`lambda` of `Y := (aJ + D) diag(mu)`, so it is monotone in N.

*Proof.* `D' + aJ_n = E (D + aJ) E^T` (E = class-indicator matrix, `E^T E = N diag(mu)`). So its nonzero eigenvalues are those
of `N Y` (`PQ` and `QP` share their nonzero spectra).
- Y is a positive matrix with `tr Y = a |mu|`, since D has zero diagonal.
- `Y >= aJ diag(mu) =: Y0`, with `Y != Y0` because `D != 0`, and Y0 is irreducible. By strict monotonicity of the Perron root,
  `rho(Y) > rho(Y0) = a|mu|`.
- So the other m-1 eigenvalues sum to `a|mu| - rho(Y) < 0`, and some eigenvalue lambda has `Re lambda < 0`.
- Then `N lambda` is an eigenvalue of `D' + aJ_n` with real part `-> -inf`. ∎

CERTIFIED: R itself is stable at `(1/21, 11/210)`, but its uniform double `R^(2)` (n = 14) is not (ctln: Hurwitz minor 4
negative). So N* = 1.

**Lemma U'' (fastest classes must be complete).** Let `mu^(k)` have `|mu^(k)| -> inf` and `mu^(k)/|mu^(k)| -> pi`, and put
`P = supp(pi)`. If the induced digraph `H[P]` is not complete, then `H^(mu^(k))` is unstable at any fixed `(beta, t)` for all
large k.

*Proof.*
- `D' + aJ` has nonzero spectrum `|mu^(k)| spec((aJ + D) Pi_k)` with `Pi_k -> Pi = diag(pi)`.
- The columns of `(aJ + D) Pi` outside P vanish, so its nonzero spectrum is that of the positive matrix `(aJ + D_PP) Pi_P`.
  This matrix has trace a and Perron root `> a` (as in Lemma U, since `D_PP != 0`), so it has an eigenvalue `lambda` with
  `Re lambda < 0`.
- By continuity of roots, `(aJ + D) Pi_k` has an eigenvalue `lambda_k -> lambda`, so `|mu^(k)| lambda_k` has real part
  `-> -inf`. ∎

So in any fixed-parameter infinite family obtained by blowing up a fixed core, the classes of maximal growth order induce a
complete digraph. They can be merged into one growing clique (or, as in Q2 and Y, several mutually bidirected cliques). The
next-order behaviour is then the Schur complement, which for sink classes is exactly `spec(Lbar)` (Theorem FP). For a
non-sink class (`c = 1 - T != 0`) the Schur complement is `D - b 1^T - 1 c^T - b c^T / a`, and it does depend on t.

---

## 7. Conjecture P is false (REFUTED, exact; `conjP_refutation.py` -> `conjP_refutation.json`)

notes_theory §5.4, Conjecture P: "if `D != 0` and `(D^T + rI)^{-1} 1 > 0` for some `r in (0, 1/2)`, then D has an eigenvalue
with `Re <= -1/2`." It was supported by 38,452 random and adversarial float samples.

Counterexample: D = non-edge matrix of `R = &F^U|m^QTq?`, `r = 10/21`.
- `chi_D = z^2 (z^5 - 3z^2 - 2z - 1)`.
- `(-1)^7 chi_D(-1/2 - x)` and `(-1)^7 chi_D(-10/21 - x)` are Hurwitz (exact), so every eigenvalue has `Re > -10/21 > -1/2`.
  The min Re is about -0.4710, display only.
- `u = (D^T + (10/21) I)^{-1} 1 > 0` (exact), with `U = 10764075/3072781 > 1`.

Consequence: the note's predicted corollary "joins `G v K_k` fail for large k" is false too. The join `R v K_k` (S = all, b = 0,
so `Lbar = D`) is an SFP at `(1/21, 11/210)` for **every** `k >= 0`. The minors are positive on `[0, inf)`, and ctln certifies
k = 1, 5, 20.

The exact catalogue scan (§9) found 9 other bases with the same behaviour (S = all, k0 = 0 at that point). Several have
different `chi_D`: `z^7 - 6z^4 - 8z^3 - 8z^2 - 4z - 1` (min Re about -0.411), `z^7 - 8z^4 - ...`, `z^7 - 9z^4 - ...`,
`z^7 - 10z^4 - ...`. The n7 catalogue of notes_theory already contained these join families, proved limit-stable for all k,
which contradicts P; nobody had cross-checked. This belongs in ANOMALIES.md (§10).

---

## 8. What remains open, and partial results

- **O1 (the full fixed-parameter map).** Is the set Lambda of legal `(eps, delta)` carrying an infinite non-clique SFP family
  equal to the whole legal region?
  - PROVED: `Lambda ⊇ {legal : 2 < beta < beta_Y = 4.9607...}`, and it is t-free for every sink/join-class family (Theorem FP).
  - For beta >= beta_Y, nothing is known either way.
  - No upper bound on n at any fixed legal point is known. Geneson's theorem only gives a lower bound,
    `n > (q + 2 - 3 eps/2)/(1 - eps/2)`.
  - For sink-class families over an m-node core with no join class (every `b_j != 0`, so `|bbar|_1 >= 1`),
    `tr Lbar = -|bbar|_1 <= -1` caps `beta_inf <= m` (notes_theory §5.5).
  - For joins (`bbar = 0`, so `Lbar = D`) the cap is `m - 1`. Positivity forces D to contain a directed cycle:
    `r u_i + sum_{j -> i in D} u_j = 1`, so a node without D-in-edge has `u_i = 1/r > 1` and cannot have a D-out-neighbour
    (that neighbour's row would force `u < 0`); walking backwards from any D-edge then closes a cycle. Hence `rho(D) >= 1`,
    and with `tr D = 0` the other eigenvalues sum to `<= -1`, so `min Re <= -1/(m-1)`.
  - Mixtures of join and non-join classes are not covered by either cap.
  - So Lambda = everything would need cores of unbounded size whose ceilings grow. This is the fixed-parameter form of the
    "is q_n* bounded?" question. Conjecture P, now refuted, was the only proposed obstruction for joins.
  - Nesting does not help, at least in the iterated limit (new class size `k' -> inf` first, then old size `k -> inf`;
    PROVED for that order). Put a new sink class, in-set complement `b' = (b'_G, b'_V)` constant on old classes, on top of a
    sink-class family. By Lemma W, `D_{G_k} - b' 1^T = E (D_+ - b'_+ 1^T) E^T`.
    - If some old class does not send to the new one (`b'_V != 0`), the old-class block contributes the eigenvalue
      `-k rho^T b'_V -> -inf`.
    - Otherwise the matrix is block upper-triangular and the limit is `spec(D - b'_G 1^T) u {0}`: the old classes drop out.
    So larger ceilings need new non-class structure in the core.
  - HEURISTIC: a 15-minute simulated-annealing search from Y's core (`heur_ceiling_search.py 0 15`; flips, in-set changes,
    class ratios, node additions up to 12 nodes; floats) found nothing above 4.9607.
- **O2.** Theorem FP is silent on the boundary `min Re spec(Lbar) = -1/beta` (an irrational beta for R, Q, Y; beta = 4 for Q2).
  There the k -> inf behaviour depends on the next-order term (as in Lemma AS). That is undecided here, never "unstable".
- **O3.** Families that are not blow-ups of a fixed core, i.e. whose core grows, are not covered by Lemmas U/U''.
- Partial impossibility results, all PROVED: Lemma U (uniform blow-ups), Lemma U'' (fastest classes complete), FP(ii)
  (spectral obstruction; kills F at every fixed point), and FP(iii) (positivity obstruction, independent of k and t).

---

## 9. Exact catalogue scan (CERTIFIED computations; float preselection HEURISTIC)

`heur_scan_catalogue.py` (floats) kept the 130 of the 626 proved catalogue families (notes_theory §3.4) whose limit polynomial
has float min Re > -1/2. `scan_k0_catalogue.py` -> `scan_k0_catalogue.json` then computed, exactly, positivity and the exact
`k0` (all minors > 0 on `[k0, inf)`) at four fixed points:

| point | families positive | eventually SFP (k0 exists) | SFP for all k >= 0 |
|---|---|---|---|
| (21/10, 20) = (1/21, 11/210) | 119 | 96 (32 distinct bases) | 45 |
| (11/5, 20) | 114 | 61 | 2 |
| (23/10, 20) | 113 | 43 | 2 |
| (12/5, 100) | 97 | 30 | 0 |

So fixed-parameter infinite families are abundant at small beta, not a rare accident.

---

## 10. Items for ANOMALIES.md (not appended: the repository was not to be modified)

1. **Conjecture P is false** (§7). The refuting cores were already in the proved n7 catalogue: join families `S = all`,
   proved limit-stable for all k, are incompatible with P. The 38k-sample float evidence missed structured cases again (cf. SL).
2. `F_39, F_40, F_41` have **bounded** stable t-windows at beta = 21/10, namely `(10, 87.98)`, `(10, 28.25)`, `(10, 16.83)`.
   They are finite-eps-only at this beta, a second mechanism next to Lemma B's cube-root boundary.
3. The limit polynomials of the Geneson family (`z^2+z+1`) and of R/Q (`z^4+z^3+3z^2+z+1`, which is `w^2+w+1` in
   `w = z + 1/z`) are both governed by the cube roots of unity. Only the first sits on Re = -1/2.
4. The first version of `fixlib.to_d6` wrote a wrong digraph6 header for n >= 63. Only stored d6 strings were affected (no
   computation); they were corrected by `fix_d6.py` and validated with nauty `showg`.
5. `notes/referee_wf4/lit_diff/` is empty in this checkout. The statement about the preprint's "5r + 6" family was taken from
   the task text only, and nothing here depends on it.

---

## 11. Files (all in PF/)

Libraries:
- `fixlib.py`: Lemma P/S machinery, Hurwitz minors in k or (a, k), positivity, explicit builders, d6.
- `theorem_B.py`: also holds FAMILIES.

Proof scripts and exact certificate files:
- `theorem_A.py` -> `theorem_A_certificate.json` (Theorem A)
- `theorem_B.py [R|Q Q2|Y]` -> `theorem_B_certificate_R.json`, `theorem_B_certificate_Q_Q2.json`,
  `theorem_B_certificate_Y.json` (Theorems A+, B)
- `ceilings.py` -> `ceilings.json` (limit polynomials, ceilings, Hurwitz checks)
- `positivity_ranges.py` -> `positivity_ranges.json` (exact positivity of R, Q on [2,3], Q2 on [2,4], Y on [2,5])
- `pointwise_thresholds.py` -> `pointwise_thresholds.json` (exact k0 at single points, including eps = 1/2)
- `verify_all.py [--full]` (round-trip of every instance record in cert_*.json, G1 at (1/21, 11/210), join instances)
- `negative.py` -> `negative_certificate.json` (N1, Lemma U instance)
- `Tk_table.py` -> `Tk_table.json` (question (a))
- `conjP_refutation.py` -> `conjP_refutation.json`
- `certify_members.py <job>` -> `cert_R_eps1_21.json`, `cert_R_eps1_8.json`, `cert_Q2_eps1_2.json` (instance certificates;
  `fix_d6.py` corrects stored d6 strings)

Tests of the proofs:
- `check_lemmaP.py` -> `check_lemmaP.json`
- `check_X2.py` -> `check_X2.json`
- `crosscheck_region.py` -> `crosscheck_region.json`

Exact computation with float preselection:
- `scan_k0_catalogue.py` -> `scan_k0_catalogue.json`

HEURISTIC (floats only; their outputs are not certificates):
- `heur_scan_catalogue.py` -> `heur_catalogue_good.json`
- `heur_R_region.py`, `heur_R_beta1.py` (the float optimum of the all-k region for R is about beta = 2.162)
- `heur_ceiling_search.py` -> `heur_ceiling_0.json`, `heur_ceiling_0.log`
- `theorem_B.heur_k0`, which only chooses the first k0 to try

Scratch: `t0.py` to `t6.py` (exploration; superseded by the scripts above).
