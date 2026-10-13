# Lemma B' (two-parameter sign-invariance), Corollary B'' (CAD correctness), Lemma B''' (projection sufficiency)

Prover notes, 2026-10-13. This closes the documentation gap found by the external audit. `region2d.py` decides
the verdict on each open cell of a cylindrical algebraic decomposition (CAD) of the `(beta, t)` plane, but
`notes/lemmas.md` proves Lemma B only for fixed `t`, i.e. one parameter. The conventions are CLAUDE.md §2,
unchanged. This note changes no code and no §3 policy.

Status labels are as in `notes/lemmas.md`:

- **PROVED**: the proof is here, or in the cited project note at the cited place.
- **CITED**: a classical theorem, with a reference.
- **UNPROVED / TRUSTED**: labelled as such.
- **SANITY**: an exact computation that supports a claim. It is a test, not a proof. Scripts are in
  `/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverB/`.

| item | status | rests on |
|---|---|---|
| 0.1 coordinates; N (no boundary polynomial is identically 0); Z; Conn; O' (openness in `(beta,t)`) | PROVED | Lemma A, A.2, (b0) |
| (b3) restated; R2 (root counts in disks, any parameter space) | PROVED | (R) of `notes/lemma_I_spectral.md` §1 |
| **B'** (sign-invariance on connected sets) and B'-comp (component form) | PROVED | A.1, (b3), R2 |
| **B''** (1)-(5) (CAD correctness, density, undecided samples, exact description of the stable set) | PROVED | B', O', (b0), (b1) |
| **B'''** (a)-(g) (delineability over open beta-intervals; sectors are sign-invariant cells; density) | PROVED, elementary; Collins / McCallum CITED for context only | R2, Sylvester kernel lemma |
| §4: what `region2d.py` guarantees | code reading + SANITY; the trust base is listed explicitly | |

---

## 0. Setting, and what is taken from other notes

Everything in this section is from `notes/lemmas.md` §0-§2, restated without change, except 0.1, N, Z, Conn, O'
and R2, which are proved here.

**Coordinates.** `t := (1-eps)/eps` and `beta := 1 + delta/eps`. Put

    L(beta, t) := t (beta - 2) - 1,        R := { (beta, t) in R^2 : beta > 2, t > 1/(beta - 2) }.

**0.1 (PROVED).** `Phi(eps, delta) := (1 + delta/eps, (1-eps)/eps)` is a homeomorphism of `(0,1) x (0,inf)` onto
`(1,inf) x (0,inf)`. Its inverse is `eps = 1/(1+t)`, `delta = (beta-1)/(1+t)`. It maps the legal region
`{eps > 0, delta > 0, q > 1}` onto `R`, so `R` is open and the standing hypothesis (H) of `lemmas.md` holds at every
point of `R` (for every `theta > 0`). Moreover

    R = {beta > 2} ∩ {L > 0},     {L > 0} = R ⊔ R',     R' := {beta < 2, t < 1/(beta-2)},

and `R`, `R'` are disjoint open sets.

*Proof.*
- Both `Phi` and its inverse are rational with denominators that do not vanish on the respective domains, so both are
  continuous. Mutual inverseness is the computation in `lemmas.md` §0.
- Legal forces `eps < 1`, since `eps >= 1` gives `q <= 0`. So the legal region lies in the domain of `Phi`.
- By `lemmas.md` §0, legal iff `t (beta - 2) > 1`. For `t > 0` this forces `beta > 2` and then `t > 1/(beta-2)`.
  Conversely, `beta > 2` and `t > 1/(beta-2)` give `t > 0`, `beta > 1` and `L > 0`.
- On the line `beta = 2` we have `L = -1`. So `{L > 0} = ({L>0} ∩ {beta>2}) ⊔ ({L>0} ∩ {beta<2})`. The first set is `R`.
  For the second, `beta - 2 < 0`, so `t(beta-2) > 1` iff `t < 1/(beta-2)`. Both sets are open. ∎

**Support data** (`lemmas.md` §0-§1). Support `sigma`, `m = |sigma| >= 1`, `N`, `C = I + beta N`, `K = C + t J`.

| polynomial | ring | definition |
|---|---|---|
| `a_i` (i in sigma) | `Z[beta]` | `a = adj(C) 1` |
| `o_k` (k outside sigma) | `Z[beta]` | `o_k = d - beta sum_b N_kb a_b` |
| `d`, `s` | `Z[beta]` | `d = det C`, `s = 1^T a` |
| `D` | `Z[beta, t]` | `D = det K = d + t s` |
| `p(mu; beta, t)` | `Z[beta, t][mu]` | `det(mu I + K) = sum_(l=0..m) c_l(beta,t) mu^(m-l)`, with `c_0 = 1` and `c_m = D` |
| `Delta_k` | `Z[beta, t]` | `det[c_(2i-j)]_(i,j=1..k)` (`c_l := 0` for `l` outside `0..m`), `Delta_0 := 1` |

The univariate `a_i`, `o_k` are also regarded as functions on the plane that do not depend on `t`. Further:

- `SFP(beta, t)` means: sigma supports a stable fixed point at `eps = 1/(1+t)`, `delta = (beta-1)/(1+t)` (any `theta > 0`).
- `S := { (beta,t) in R : SFP(beta,t) }` is the **stable set**.
- `Z(h)` is the real zero set of `h` in `R^2`.
- `Z_0 := Z(D) ∪ Z(Delta_(m-1)) ∪ U_i Z(a_i) ∪ U_k Z(o_k)`.
- `B+ := { beta > 2 : a_i(beta) > 0 for all i, o_k(beta) < 0 for all k }` (the positivity set), and `R_+ := R ∩ (B+ x R)`.

**Facts used, with where they are proved.**
- **(A.1)** (`lemmas.md` Cor. A.1) At `(beta,t) in R`: `SFP <=> D != 0, a_i/D > 0 (all i), o_k/D < 0 (all k), p Hurwitz`.
- **(A.2)** (`lemmas.md` Lemma A.2) If `p` is Hurwitz then `D > 0`. Hence at `(beta,t) in R`:
  `SFP <=> a_i > 0 (all i), o_k < 0 (all k), p Hurwitz`. In particular `S ⊂ R_+`.
- **(b0)** CITED (Routh-Hurwitz; `lemmas.md` B.1). For a real polynomial with positive leading coefficient: Hurwitz
  iff `Delta_1 > 0, ..., Delta_n > 0`. Both directions are used.
- **(b1)** PROVED in `lemmas.md` B.1: Hurwitz implies every coefficient is `> 0`.
- **(b2)** PROVED in `lemmas.md` B.1: `Delta_n = c_n Delta_(n-1)`.
- **(R)** PROVED in `notes/lemma_I_spectral.md` §1: continuity of the roots of monic polynomials. The proof is
  self-contained (Cauchy bound + compactness + unique factorisation) and uses no complex analysis.

**(b3) restated, PROVED (the kernel-vector form of Orlando's formula; `lemmas.md` B.1 (b3)).** Let
`f(z) = c_0 z^n + ... + c_n` be real with `c_0 != 0` and `n >= 2`. If `f(iw) = 0` for some real `w != 0`, then
`Delta_(n-1)(f) = 0`.

*Proof.*
1. Put `P(x) := sum_(l=0..n) c_l x^l = x^n f(1/x)`. Then `P(0) = c_0 != 0`.
2. `f` is real, so `f(iw) = f(-iw) = 0`. Hence `P(-i/w) = (iw)^(-n) f(iw) = 0` and likewise `P(i/w) = 0`. These two
   distinct roots are the roots of the real quadratic `1 + w^2 x^2`. So `1 + w^2 x^2` divides `P` in `C[x]`, hence
   in `R[x]` (the remainder of division in `R[x]` is the same as in `C[x]`, namely 0).
   Write `P = (1 + w^2 x^2) R` with `R` real, `deg R <= n-2`, `R(0) = c_0 != 0`.
3. Put `V(x) := x R(-x) = sum_(j=1..n-1) v_j x^j`. Then `V != 0`, because `v_1 = R(0) != 0`.
4. `P(x) V(x) = x (1 + w^2 x^2) R(x) R(-x)` is an odd polynomial: `R(x)R(-x)` and `1 + w^2 x^2` are even. So every
   even-degree coefficient of `P V` vanishes.
5. For `i = 1..n-1`: `[x^(2i)](P V) = sum_(j=1..n-1) c_(2i-j) v_j = (H_(n-1) v)_i`. Here `H_(n-1) = [c_(2i-j)]_(i,j=1..n-1)`
   (with `c_l = 0` outside `0..n`) is exactly the matrix whose determinant is `Delta_(n-1)(f)`.
6. So `H_(n-1) v = 0` with `v = (v_1..v_(n-1)) != 0`, and `Delta_(n-1)(f) = det H_(n-1) = 0`. ∎

**Lemma N (no boundary polynomial is the zero polynomial), PROVED.**
- `a_i(0) = 1` and `o_k(0) = 1`.
- `D(0,0) = 1`.
- `c_l(0,0) = binom(m, l)`.
- `Delta_k(0,0) > 0` for `k = 1..m`.

So `a_i`, `o_k`, `D`, every `c_l` and every `Delta_k` are nonzero polynomials. The same holds for fixed `t >= 0`:
`D_t(0) = 1 + m t > 0` and `Delta_(m-1)(p_t)(0) > 0`. Hence the "zero polynomial" branch of `lemmas.md` Lemma B(d)
is never taken.

*Proof.*
- At `beta = 0`, `C = I`. So `a = adj(I) 1 = 1`, `d = 1`, and `o_k = d - 0 = 1`.
- At `(0, 0)`, `K = I`, so `p = (mu + 1)^m`. Its coefficients are binomial, and its only root is `-1`, so it is
  Hurwitz. By (b0) `=>`, all `Delta_k > 0`.
- For fixed `t >= 0` and `beta = 0`: `K = I + tJ` has eigenvalues `1` (multiplicity `m-1`) and `1 + mt`. So
  `p_t = (mu+1)^(m-1) (mu + 1 + mt)`, which is Hurwitz, and `D_t(0) = 1 + mt`. ∎

(The point `(0, 0)` is not legal. That does not matter: these are identities between polynomials.)

**Lemma Z, PROVED.** A nonzero `g in R[beta, t]` vanishes identically on no nonempty open subset of `R^2`. Hence
`Z(g)` is a closed set with empty interior. In particular `Z_0 = Z(prod_i a_i prod_k o_k D Delta_(m-1))` is closed
with empty interior, by Lemma N.

*Proof.*
1. A nonempty open set contains a box `I_1 x I_2` of open intervals. Write `g = sum_k g_k(beta) t^k`.
2. Suppose `g` vanishes on the box. For each `beta in I_1`, `g(beta, .)` vanishes on the infinite set `I_2`, so
   every `g_k(beta) = 0`.
3. Then each `g_k` vanishes on the infinite set `I_1`, so `g_k = 0`, and `g = 0`. ∎

**Lemma Conn, PROVED.** A continuous real function without zeros on a connected topological space `X` has constant
sign there.
*Proof.* `X = f^-1(0, inf) ⊔ f^-1(-inf, 0)` is a union of two disjoint open subsets, so one of them is empty. ∎

**Lemma O' (the stable set is open; Lemma O in these coordinates), PROVED.**

    S = R ∩ { a_i > 0 for all i } ∩ { o_k < 0 for all k } ∩ { Delta_1 > 0, ..., Delta_m > 0 }.

In particular `S` is open.

*Proof.* By (A.2), `SFP <=> a_i > 0, o_k < 0, p Hurwitz`. Since `p` is monic, (b0) gives
`p Hurwitz <=> Delta_1, ..., Delta_m > 0`. Each set on the right is the preimage of an open set under a polynomial,
hence continuous, map, and `R` is open (0.1). ∎

This is `notes/lemma_F_reduction.md` Lemma O transported through the homeomorphism `Phi` of 0.1. The proof above is
direct and does not need the transport.

**Lemma R2 (root counts in disks, any parameter space), PROVED from (R).** Let `X` be a metric space and `x_0 in X`.
Let `q(mu; x) = mu^n + e_1(x) mu^(n-1) + ... + e_n(x)` with `n >= 1` and each `e_l : X -> C` continuous. Let `z` be a
root of `q(.; x_0)` of multiplicity `r`, and let `0 < rho < dist(z, other roots of q(.; x_0))` (`+inf` if there are
none). Then there is `eta > 0` such that for every `x in X` with `dist(x, x_0) < eta`, the open disk `B(z, rho)`
contains exactly `r` roots of `q(.; x)`, counted with multiplicity.

*Proof.*
1. Suppose not. Then there is a sequence `x_k -> x_0` such that, for every `k`, the disk does not contain exactly `r`
   roots of `q(.; x_k)`.
2. The coefficients `e_l(x_k)` converge to `e_l(x_0)` by continuity.
3. The "Consequently" clause of (R) then says that for all large `k` the disk contains exactly `r` roots of
   `q(.; x_k)`. This contradicts step 1. ∎

(`lemmas.md` B.2 states the one-parameter case by citing Rouché's theorem. R2 needs neither one parameter nor
Rouché.)

---

## 1. Lemma B' (two-parameter sign-invariance)

**Lemma B' (PROVED).** Let `m >= 1` and fix the support `sigma`. Let `U ⊂ R` be connected. Assume that none of the
following vanishes at any point of `U`:

- `a_i(beta)` (i in sigma);
- `o_k(beta)` (k outside sigma);
- `D(beta, t)`;
- `Delta_(m-1)(beta, t)`.

Then the truth value of `SFP(beta, t)` is constant on `U`.

The task statement asks for `U` connected and open. **Openness is not used.** The lemma holds for every connected
`U ⊂ R`. Two special cases:
- `U = I x {t}` is `lemmas.md` Lemma B(d).
- `U = {beta} x I` is the "roles swapped" remark after Lemma B.

*Proof.*

**(i) The first three clauses of (A.1) are constant on `U`.**
- `D` is a polynomial, hence continuous on `U`, and has no zero there. By Conn its sign is constant on `U`, and
  `D != 0` holds on all of `U`.
- For each `i`, `(beta, t) -> a_i(beta)/D(beta, t)` is continuous and nonvanishing on `U`. By Conn, the truth value of
  `a_i/D > 0` is constant on `U`.
- The same argument applies to each `o_k/D`.

**(ii) Stability: "`p(.; beta, t)` is Hurwitz" is constant on `U`.**
- *No root on the imaginary axis.* Let `(beta, t) in U` and let `p(.; beta, t)` have a root `zeta` with
  `Re zeta = 0`.
  - If `zeta = 0`, then `D(beta, t) = p(0; beta, t) = 0`. This is excluded.
  - If `zeta = iw` with `w != 0`, then `m >= 2`: for `m = 1`, `p` has degree 1 and its only root is real. Then
    (b3) applied to the real polynomial `p(.; beta, t)` (leading coefficient 1, degree `m >= 2`) gives
    `Delta_(m-1)(beta, t) = 0`. This is excluded too.

  So for every `(beta, t) in U`, `p(.; beta, t)` has no root on `iR`.
- *The count `n_+(beta, t)` of roots with positive real part (with multiplicity) is locally constant on `U`.*
  1. Fix `x_0 = (beta_0, t_0) in U`, and let `z_1, ..., z_r` be the distinct roots of `p(.; x_0)`. By the previous
     point, `Re z_l != 0` for all `l`.
  2. Choose `rho > 0` smaller than half of every distance `|z_l - z_l'|` (`l != l'`), and smaller than every
     `|Re z_l|`. Then the disks `B(z_l, rho)` are pairwise disjoint, and each lies in one open half-plane
     `{Re < 0}` or `{Re > 0}`.
  3. `p` is monic of degree `m` in `mu`, with coefficients `c_l(beta, t)` that are polynomials and hence continuous on
     `R^2` (`lemmas.md` A(iv)). Apply R2 with `X = R^2` to each `z_l`, and let `eta` be the smallest of the resulting
     radii. Then for every `x` with `|x - x_0| < eta`, each `B(z_l, rho)` contains exactly `mult(z_l)` roots of
     `p(.; x)`.
  4. These multiplicities sum to `m`, so they account for all roots of `p(.; x)`. Hence
     `n_+(x) = sum over {l : Re z_l > 0} of mult(z_l) = n_+(x_0)`.
  5. This holds for every `x` in an open disk around `x_0`, in particular for every `x` in that disk ∩ `U`. So `n_+`
     is locally constant on `U` in the subspace topology.
- *Conclusion.* The level sets `{x in U : n_+(x) = k}` (`k = 0..m`) are open in `U`, pairwise disjoint, and cover
  `U`. Since `U` is connected, exactly one is nonempty, so `n_+` is constant on `U`. No root lies on `iR` anywhere on
  `U`, so the number of roots in the open left half-plane, `n_- = m - n_+`, is also constant on `U`. Then
  `p Hurwitz <=> n_- = m`, which has constant truth value on `U`.

**(iii) Conclusion.** By (A.1), `SFP` is the conjunction of the clauses in (i) and (ii). Each has constant truth value
on `U`, so `SFP` does too. ∎

**Corollary B'-comp (component form), PROVED.**
- (a) `S ∩ Z_0 = ∅`.
- (b) `R \ Z_0` is open. Each connected component `V` of `R \ Z_0` is open and connected, and either `V ⊂ S` or
  `V ∩ S = ∅`.

So `S` is a union of connected components of `R \ Z_0`.

*Proof.*
- (a) At a point of `S`, (A.2) gives `a_i > 0`, `o_k < 0` and `D > 0`, and (b0) gives `Delta_(m-1) > 0`. So the
  point is not in `Z_0`.
- (b) `R` is open and `Z_0` is closed, so `R \ Z_0` is open. Connected components of open subsets of `R^2` are open,
  because open disks are connected. On `V`, none of `a_i`, `o_k`, `D`, `Delta_(m-1)` vanishes, so Lemma B' applies. ∎

Remarks.
- Only `Delta_(m-1)` is needed, not `Delta_(m-2)` or a Lienard-Chipart subset. The imaginary-axis crossing at
  `+-iw` with `w != 0` is detected exactly by `Delta_(m-1)`, via (b3), and the crossing at 0 by `c_m = D`.
  `region2d.curves_lc` uses exactly `{D, Delta_(m-1), L}`, in agreement with this.
- The hypotheses on `a_i` and `o_k` could be weakened using (A.2). They are kept as the task states them.
- Lemma B' is about arbitrary connected sets that avoid `Z_0`. It says nothing yet about the cells of a particular
  decomposition, nor about whether every stable point is seen by some cell. Those are Corollary B'' and Lemma B'''.

---

## 2. Corollary B'' (correctness of deciding by open cells)

**Setting.**
- `C ⊂ R[beta, t]` is a finite set of polynomials containing `a_i` (all i in sigma), `o_k` (all k outside sigma),
  `D`, `Delta_(m-1)` and `L`.
- A finite family `{U_j}_(j in J)` of nonempty, connected, open, pairwise disjoint subsets of `R^2` is
  **C-sign-invariant** if no member of `C` vanishes at any point of any `U_j`. By Conn this is equivalent to: each
  member of `C` has a constant nonzero sign on each `U_j`.
- `J_R := { j : U_j ∩ R != ∅ }`.
- **Density hypothesis (DH):** every nonempty open subset of `R_+` meets `U_(j in J_R) U_j`.
  - Recall `R_+ = R ∩ (B+ x R)` and `S ⊂ R_+` by (A.2).
  - (DH) holds in particular if `R_+ \ U_j U_j` has empty interior.
  - It is the only property of the decomposition, beyond sign-invariance, that the conclusions below use.
  - Theorem B''' (iii) (§3.5) proves (DH) for the decomposition built by `region2d.py`.
  - (DH) also holds for the 2-cells of **any** CAD of `R^2` (PROVED). The other cells are finitely many points,
    vertical open segments or lines, and graphs of continuous functions over open intervals, and each of these is
    nowhere dense.
    - For a graph `Gamma` of a continuous `xi` on an open interval `I'`: a limit of points `(beta_k, xi(beta_k))` with
      limit `beta in I'` is `(beta, xi(beta))`. So any open box inside `cl(Gamma)` could be shrunk to one over a
      subinterval of `I'`, and would then lie in `Gamma` itself. That is impossible, since a box contains two points
      with the same `beta`.
    - A finite union of nowhere dense sets is nowhere dense, so the union of the 2-cells is dense in `R^2`. Every
      nonempty open subset of `R_+` then meets a 2-cell, which by (1) below lies in `R`.

**Corollary B'' (PROVED).**

**(1) Each cell is inside `R` or disjoint from it.** For every `j`, either `U_j ⊂ R` or `U_j ∩ R = ∅`.

**(2) Constant verdict on a cell.** For `j in J_R`, `SFP` is constant on `U_j`. So the verdict on `U_j` equals the
exact verdict at any single point of `U_j`.

**(3) Undecided samples.** Let `j in J_R` and `x in U_j`, and suppose some coefficient `c_l(x) = 0`
(`1 <= l <= m`) or some Hurwitz minor `Delta_k(x) = 0` (`1 <= k <= m`). This covers every point at which
`ctln.hurwitz_stable` returns `None` (*undecided*) on `p(.; x)` or on any positive multiple of it. The converse is
false: if an earlier minor is negative, `hurwitz_stable` returns `False` before it reaches a vanishing one. Then:

- **(a)** The vanishing quantity is some `c_l` with `1 <= l <= m-1`, or some `Delta_k` with `1 <= k <= m-2`. Indeed
  `c_0 = 1`, while `c_m = D`, `Delta_(m-1)` and `Delta_m = D Delta_(m-1)` (b2) do not vanish on `U_j`. So an
  undecided sample does not contradict sign-invariance: the polynomials `c_l` (`l <= m-1`) and `Delta_k`
  (`k <= m-2`) are not in `C`, and their zero sets may cross the interior of a cell.
- **(b)** `SFP` is false at `x`, and therefore, by (2), on all of `U_j`.
- **(c)** Every other point `y in U_j` has the cell's verdict. So re-sampling inside the same cell is legitimate, and
  any decided re-sample must return False. A re-sample returning True would contradict (b) together with (2), and
  would therefore indicate a fault upstream: a wrong projection, or a failure of the trust base in §4.
- **(d)** Equivalently: on a cell where `SFP` holds at one point, every `c_l > 0` and every `Delta_k > 0` at
  *every* point of the cell. The zero sets of the `c_l` and of the lower minors only ever cross non-stable cells.
  This is why they need not be in `C`.

**(4) Lower-dimensional pieces never need sampling.** Assume (DH). If `S != ∅`, then `S` meets some `U_j` with
`j in J_R`, and then `U_j ⊂ S`. Equivalently:

    no cell U_j (j in J_R) has a stable sample   ==>   S = ∅  (never stable at any legal parameter).

**(5) Exact description of `S` from the stable cells.** Assume (DH). Put

    S_c := U { U_j : j in J_R, U_j ⊂ S }   (the union of the stable cells).

Then `S_c ⊂ S ⊂ cl(S_c)`, and more precisely

    S = ( R ∩ cl(S_c) ) \ Z_0.

So the stable cells determine `S` completely. A point of `R` on a lower-dimensional piece is stable iff it lies off
`Z_0` and in the closure of the stable cells.

*Proof.*

- (1) `L in C` has a constant nonzero sign on `U_j`. If `U_j` meets `R`, that sign is `+`, so `U_j ⊂ {L > 0} = R ⊔ R'`
  (0.1). `R` and `R'` are disjoint open sets and `U_j` is connected, so `U_j` lies in one of them. It meets `R`, so
  `U_j ⊂ R`.
- (2) By (1), `U_j ⊂ R`. It is connected, and `a_i`, `o_k`, `D`, `Delta_(m-1)` (all in `C`) do not vanish on it, so
  Lemma B' applies.
- (3)
  - (a) is stated with its proof.
  - (b) If some `c_l(x) = 0`, then `p(.; x)` is not Hurwitz by (b1). If some `Delta_k(x) = 0`, it is not Hurwitz by
    (b0), direction `=>`. By (A.1), `SFP` is false at `x`, and by (2) on all of `U_j`.

    *That the hypothesis covers `hurwitz_stable = None`.* Read off `ctln.hurwitz_stable` (`ctln.py`, lines 199-230).
    With leading coefficient `> 0`, it returns `None` in exactly two situations. Either some coefficient is 0, or
    all coefficients are `> 0` and the first non-positive minor in the order `Delta_1, ..., Delta_m` is 0.

    A positive multiple `lambda p` (`fastregion.stability_coeffs` multiplies by the denominator of `t`) has
    coefficients `lambda c_l` and minors `lambda^k Delta_k`, so zero patterns and signs are the same.
  - (c) follows from (2).
  - (d) If `SFP` holds at one point of `U_j`, it holds on all of `U_j` by (2). Then `p` is Hurwitz at every point of
    `U_j`, so every `c_l > 0` (b1) and every `Delta_k > 0` (b0).
- (4) `S` is open (Lemma O') and nonempty, and `S ⊂ R_+` (A.2). By (DH) it meets some `U_j`, `j in J_R`, at a point
  `y`. `SFP(y)` is true, so by (2) `SFP` holds on all of `U_j`.
- (5)
  - `S_c ⊂ S` by definition.
  - `S ⊂ cl(S_c)`: let `x in S` and let `O` be an open neighbourhood of `x`. Then `O ∩ S` is a nonempty open subset of
    `R_+`. By (DH) it meets some `U_j`, `j in J_R`, at a point of `S`. By (2) that `U_j ⊂ S_c`. So `O` meets `S_c`.
  - `S ⊂ (R ∩ cl(S_c)) \ Z_0`: `S ⊂ R`, `S ⊂ cl(S_c)`, and `S ∩ Z_0 = ∅` (B'-comp (a)).
  - Conversely, let `x in R ∩ cl(S_c)` with `x` not in `Z_0`. `R \ Z_0` is open (0.1, and `Z_0` is closed), so
    some open disk `Q ∋ x` lies in `R \ Z_0`. `Q` is connected and avoids the zeros of `a_i`, `o_k`, `D` and
    `Delta_(m-1)`, so by Lemma B' `SFP` is constant on `Q`. Since `x in cl(S_c)`, `Q` meets `S_c ⊂ S`, so the
    constant value is True. Hence `x in S`. ∎

Remarks.
- (2) and (3) need no density hypothesis. (4) and (5) need (DH) and openness (Lemma O'). Without openness, a stable
  set lying entirely on a curve would be invisible to every full-dimensional cell. Lemma O' rules that out.
- (3) is a statement about **cells**. It does not touch CLAUDE.md §3. At an arbitrary sample point a vanishing minor
  says nothing about the neighbourhood, so a pointwise sweep must still call it *undecided*. Inside a sign-invariant
  cell, (2) has already settled the neighbourhood question. `region2d.py` currently re-samples and then escalates,
  which is conservative. Whether to use (3) to close such cells as "not stable" is a policy decision, not a
  mathematical one. The point is moot for the existing sweeps: all three n=6 shard tallies and the n=5 tally report
  `escalate: 0`.
- **Consistency check implied by (3)(c)** (a hardening suggestion; nothing was changed). A cell whose first sample is
  undecided must never produce a True verdict at a re-sample. If it ever did, that would be a contradiction, and the
  run should stop, exactly as for a counterexample.

---

## 3. Lemma B''' (the projection set suffices: sectors are sign-invariant cells)

### 3.1 Setting

- **Curves.** `F ⊂ Q[beta, t]` is a finite set of polynomials, each with `deg_t f >= 1`.
- **Vertical polynomials.** `V ⊂ Q[beta]` is a finite set of nonzero polynomials.
- `lc_t f` is the coefficient of `t^(deg_t f)` in `f`. It is a nonzero element of `Q[beta]`.
- `res_t(f, g) := det Syl_t(f, g)`. This is the Sylvester determinant with respect to the formal degrees
  `p = deg_t f`, `q = deg_t g`, whose entries are the `t`-coefficients of `f` and `g`. Only its **zero set** is used
  below, so a global sign convention is irrelevant. (This matters: sympy 1.14 returns `-det Syl` when `p < q` and
  `pq` is odd; see the Appendix.)
- `disc_t f` (for `d := deg_t f >= 2`) is defined by `res_t(f, ∂_t f) = (-1)^(d(d-1)/2) lc_t(f) disc_t(f)` in
  `Q[beta]`. The quotient exists: the first column of `Syl(f, ∂_t f)` has only two nonzero entries, `lc_t f` (row 0)
  and `d lc_t f` (row `d-1`), so the Laplace expansion along it is divisible by `lc_t f`. This is sympy's convention;
  SANITY-checked as an identity of polynomials in the Appendix.
- **Projection set.**

      P := V ∪ { lc_t f : f in F } ∪ { disc_t f : f in F, deg_t f >= 2 } ∪ { res_t(f, g) : f != g in F }.

- **Hypothesis (P≠0):** no member of `P` is the zero polynomial.
- `I` is a nonempty open interval of `R` containing no real root of any member of `P`. For example, `I` can be the
  open interval between two consecutive real roots of `prod P`.

**Sylvester kernel lemma (PROVED).** Let `A` be a commutative ring and `ev : A -> C` a ring homomorphism (below:
evaluation at a real `beta`). Let `f = sum_(k=0..p) f_k t^(p-k)` and `g = sum_(k=0..q) g_k t^(q-k)` lie in `A[t]`
with `p, q >= 1`. If `ev(f)` and `ev(g)` have a common root `z in C`, then `ev(det Syl(f, g)) = 0`.
**No assumption on the leading coefficients is needed.**

*Proof.*
1. `ev` commutes with determinants, so `ev(det Syl(f,g)) = det ev(Syl(f,g))`.
2. Row `i` (`0 <= i < q`) of `Syl` carries `f_0..f_p` in columns `i..i+p`. Row `q+i` (`0 <= i < p`) carries
   `g_0..g_q` in columns `i..i+q`.
3. Let `v := (z^(p+q-1), ..., z, 1)^T != 0`. Row `i` times `v` is `z^(q-1-i) ev(f)(z) = 0`, and row `q+i` times `v`
   is `z^(p-1-i) ev(g)(z) = 0`.
4. So `ev(Syl)` has a nonzero kernel vector over `C`, and its determinant is 0. ∎

### 3.2 Statement

**Lemma B''' (PROVED).** Assume (P≠0), and let `I` be as above. Then:

- **(a) Fixed degree, simple roots.** For `f in F` and `beta in I`, `f(beta, .)` has degree exactly `d_f := deg_t f`
  and `d_f` distinct complex roots.
- **(b) Delineability of each curve.** For `f in F`, the number `r_f` of distinct real roots of `f(beta, .)` is
  constant on `I`. The ordered real roots `xi_(f,1)(beta) < ... < xi_(f,r_f)(beta)` are continuous functions on `I`.
- **(c) No crossings between curves.** For `f != g` in `F` and `beta in I`, `f(beta, .)` and `g(beta, .)` have no
  common complex root.
- **(d) Sections.** For `beta in I`, the real `t` with `prod_(f in F) f(beta, t) = 0` are
  `xi_1(beta) < ... < xi_r(beta)`.
  - `r = sum_f r_f` is constant on `I`.
  - Each `xi_j` is continuous on `I`.
  - For each `j` there is one `f_(j) in F` such that `xi_j(beta)` is a root of `f_(j)(beta, .)` for every `beta in I`.
- **(e) Sectors.** Put `xi_0 := -inf` and `xi_(r+1) := +inf`, and for `j = 0..r` let

      S_j := { (beta, t) : beta in I,  xi_j(beta) < t < xi_(j+1)(beta) }.

  - Each `S_j` is nonempty, open and connected.
  - `I x R` is the disjoint union of the `S_j` and the graphs `Gamma_j` of the `xi_j` (`1 <= j <= r`).
  - Every member of `F ∪ V` has a constant nonzero sign on each `S_j`.
  - More generally, let `h in R[beta, t]` satisfy the **covering condition**
    `Z(h) ∩ (I x R) ⊂ U_(f in F) Z(f) ∪ U_(v in V) (Z(v) x R)`. Then `h` has a constant nonzero sign on each `S_j`.
- **(f) Position relative to `R`.** Suppose `I ⊂ (2, inf)` and `F` contains a nonzero rational multiple of `L`. Then:
  - there is an index `j_L`, constant on `I`, with `xi_(j_L)(beta) = 1/(beta - 2)` for all `beta in I`;
  - `S_j ⊂ R` for `j >= j_L`, and `S_j ∩ R = ∅` for `j < j_L`;
  - `R ∩ (I x R) = U_(j >= j_L) S_j  ∪  U_(j > j_L) Gamma_j`.
- **(g) Slices.** For each `b in I`, `S_j ∩ ({b} x R) = {b} x (xi_j(b), xi_(j+1)(b))`. Under (f), the real roots of
  `prod_f f(b, .)` that are `> t_min = 1/(b-2)` cut `(t_min, inf)` into open intervals. In increasing order, these
  are the slices at `b` of `S_(j_L), S_(j_L+1), ..., S_r`.

### 3.3 Proof

- **(a)**
  - `lc_t f in P` has no root in `I`, so `deg f(beta, .) = d_f`.
  - If `d_f = 1` there is one root, which is simple.
  - Let `d_f >= 2`, and suppose `z` is a repeated root of `f(beta, .)`. Then `f(beta, z) = (∂_t f)(beta, z) = 0`.
    The Sylvester kernel lemma, with formal degrees `d_f` and `d_f - 1 >= 1`, gives `res_t(f, ∂_t f)(beta) = 0`. So
    `lc_t f(beta) disc_t f(beta) = 0`. Both factors are in `P` and nonzero at `beta`, a contradiction.
- **(b)**
  1. For `beta in I`, let `q(.; beta) := f(beta, .)/lc_t f(beta)`. It is monic of degree `d_f`, with coefficients that
     are continuous on `I` (the denominator does not vanish there).
  2. Fix `beta_0 in I`, and let `z_1..z_(d_f)` be the roots of `q(.; beta_0)`, all simple by (a).
  3. Let `rho > 0` be smaller than half of every `|z_l - z_l'|` (`l != l'`), and smaller than `|Im z_l|` for every
     non-real `z_l`.
  4. By R2 (with `X = I`) there is `eta > 0` such that for `|beta - beta_0| < eta` each disk `B(z_l, rho)` contains
     exactly one root `w_l(beta)` of `q(.; beta)`. These are all `d_f` roots.
  5. If `z_l` is real, the disk is symmetric under complex conjugation. `q` is real, so `conj(w_l(beta))` is also a
     root in that disk. By uniqueness `w_l(beta) = conj(w_l(beta))`, i.e. `w_l(beta)` is real.
  6. If `z_l` is not real, the disk does not meet `R`, because `rho < |Im z_l|`. So `w_l(beta)` is not real.
  7. Hence `r_f(beta) = r_f(beta_0)` for `|beta - beta_0| < eta`. So `r_f` is locally constant on the interval `I`.
     Its level sets are therefore open and partition the connected set `I`, so `r_f` is constant.
  8. *Continuity.* The real disks `(z_l - rho, z_l + rho)` (`z_l` real) are disjoint and ordered like their centres.
     So `xi_(f,i)(beta)` lies within `rho` of `xi_(f,i)(beta_0)` for `|beta - beta_0| < eta`. Since `rho` can be taken
     arbitrarily small, each `xi_(f,i)` is continuous at `beta_0`.
- **(c)** A common root at `beta in I` would give `res_t(f, g)(beta) = 0` by the Sylvester kernel lemma. But
  `res_t(f, g) in P`, so this is impossible.
- **(d)**
  1. The functions `xi_(f,i)` (`f in F`, `1 <= i <= r_f`) are continuous on `I` by (b).
  2. Two of them are never equal at any point: within one `f` the roots are distinct, and across different `f` use (c).
  3. For any two of them, `u - v` is continuous and nonvanishing on `I`, so its sign is constant (Conn). The relative
     order of all of them is therefore the same at every `beta in I`.
  4. Listing them in that order gives `xi_1 < ... < xi_r` with constant labels `f_(j)`.
- **(e)**
  - *Open.* `S_j = { (beta,t) in I x R : t - xi_j(beta) > 0, xi_(j+1)(beta) - t > 0 }`. This is the preimage of an
    open set under a continuous map on the open set `I x R`. The conditions involving `+-inf` are void.
  - *Connected, nonempty.* `S_j` is the image of the connected set `I x (0, 1)` (or `I x (0, inf)`) under a continuous
    map:
    - `(beta, s) -> (beta, xi_j + s (xi_(j+1) - xi_j))` when `1 <= j < r`;
    - `(beta, u) -> (beta, xi_r(beta) + u)` when `j = r >= 1`;
    - `(beta, u) -> (beta, xi_1(beta) - u)` when `j = 0 < r`;
    - `S_0 = I x R` when `r = 0`.
  - *Partition.* For `beta in I` and real `t`, either `t` equals some `xi_j(beta)`, or `t` lies strictly between two
    consecutive ones (counting `+-inf`).
  - *Signs.*
    - If `f in F` vanished at `(beta, t) in S_j`, then `t` would be a real root of `f(beta, .)`, i.e. some
      `xi_i(beta)`. No `xi_i(beta)` lies strictly between consecutive `xi_j(beta)` and `xi_(j+1)(beta)`.
    - `v in V ⊂ P` has no root in `I`.
    - Under the covering condition, `h` can vanish on `S_j` only where some `f` or `v` does, and neither does.

    So all of these are nonvanishing on the connected `S_j`, and Conn gives constant signs.
- **(f)**
  - For `beta > 2`, `L(beta, .) = (beta-2) t - 1` has exactly one root, `1/(beta-2)`. It is a real root of the member
    `c L in F`, so by (d) it is `xi_(j_L)(beta)` for an index `j_L` that does not depend on `beta in I`.
  - Points of `S_j` with `j >= j_L` satisfy `t > xi_j(beta) >= 1/(beta-2)` and `beta > 2`, so they are in `R`.
  - Points with `j < j_L` satisfy `t < xi_(j+1)(beta) <= 1/(beta-2)`, so they are not in `R`.
  - The last identity follows from the partition in (e).
- **(g)** Immediate from the definitions and (f). ∎

### 3.4 When (P≠0) and the covering condition hold

**(P≠0) is automatic for the construction stated in the task (PROVED; Gauss's lemma CITED).** Take for `F` the
irreducible factors over `Q` with `deg_t >= 1` of `D`, `Delta_(m-1)` and `L`, made pairwise non-associate; distinct
irreducibles are coprime. Put the factors with `deg_t = 0` into `V`, together with the nonconstant `a_i` and `o_k`.

- `lc_t f != 0` by definition of degree.
- Members of `V` are nonzero:
  - factors of nonzero polynomials are nonzero;
  - the `a_i` and `o_k` are nonzero by Lemma N.
- `res_t(f, g) != 0` for coprime `f, g` with `deg_t >= 1`:
  1. Let `K = Q(beta)`. If `det Syl(f, g) = 0` in `Q[beta]`, the `p + q` polynomials `t^i f` (`i < q`) and `t^j g`
     (`j < p`) are linearly dependent over `K`, because they are the rows of `Syl` in the monomial basis.
  2. So `u f + v g = 0` for some `u, v in K[t]`, not both zero, with `deg u < q` and `deg v < p`.
  3. If `gcd_(K[t])(f, g) = 1`, then `f | v`. Since `deg v < p`, this forces `v = 0`, and then `u = 0`. That is a
     contradiction.
  4. So `f` and `g` have a common factor of positive `t`-degree in `K[t]`. By Gauss's lemma its primitive part in
     `Q[beta][t]` divides both `f` and `g` in `Q[beta, t]`. That contradicts coprimality.
- `disc_t f != 0` for irreducible `f` with `deg_t f >= 2`:
  1. Irreducible with `deg_t >= 1` implies primitive in `Q[beta][t]`, and hence, by Gauss's lemma, irreducible in
     `K[t]`.
  2. `∂_t f != 0` has smaller degree, so `gcd_(K[t])(f, ∂_t f) = 1`.
  3. By the previous item `res_t(f, ∂_t f) != 0`, so `disc_t f != 0`.

Gauss's lemma: e.g. S. Lang, *Algebra*, rev. 3rd ed., Springer 2002, Ch. IV §2 (section number from memory).

**The covering condition holds for `h in {D, Delta_(m-1), L}`** whenever `h = c prod_(f in F) f^(e_f) prod_(v in V) v^(e_v)`
with `c in Q \ {0}`: `Z(h)` is the union of the `Z(f)` and the `Z(v) x R` over the factors that occur. For the task's
construction this is unique factorisation. For `region2d.py` see §4.

### 3.5 Theorem B''' (the decomposition is sign-invariant and dense)

**Theorem B''' (PROVED, given the hypotheses).** Let `F`, `V`, `P` satisfy (P≠0), and suppose:
- every nonconstant `a_i` and `o_k` is in `V`;
- `F` contains a nonzero rational multiple of `L`;
- `D`, `Delta_(m-1)` and `L` satisfy the covering condition on every strip `I x R`.

Let `B+` be the positivity set, and let `{I_alpha}` be the connected components of `B+ \ {real roots of members of P}`.
These are finitely many open intervals inside `(2, inf)`. Let `j_L(alpha)` be the index of the `L`-section over
`I_alpha`. Then the family

    { S^alpha_j : alpha,  j_L(alpha) <= j <= r(alpha) }

- (i) consists of subsets of `R`;
- (ii) is `C`-sign-invariant for `C = {a_i, o_k, D, Delta_(m-1), L}`;
- (iii) satisfies (DH).

So Corollary B'' (1)-(5) holds for it, with `J_R` = the whole family.

*Proof.*
- (i) This is B''' (f).
- (ii)
  - The sectors are nonempty, open and connected by B''' (e).
  - They are pairwise disjoint: sectors over different `I_alpha` have disjoint `beta`-ranges, and sectors over one
    `I_alpha` are disjoint by definition.
  - `a_i` and `o_k` are either constant (and then equal to 1 by Lemma N) or in `V ⊂ P`, so they do not vanish on any
    `I_alpha`.
  - `D`, `Delta_(m-1)` and `L` have constant nonzero sign on each sector by B''' (e) and the covering condition.
- (iii)
  1. Let `Pi := prod_(p in P) p(beta) * prod_(f in F) f(beta, t)`. This is a nonzero polynomial by (P≠0).
  2. Let `x = (beta, t) in R_+` with `Pi(x) != 0`. Then `beta in B+`, and `beta` is not a root of any member of `P`,
     so `beta in I_alpha` for some `alpha`.
  3. `t` is not a root of any `f(beta, .)`, so it is not a section value, and it lies in some `S^alpha_j`. Since
     `t > 1/(beta - 2)`, B''' (f) gives `j >= j_L(alpha)`.
  4. So `R_+ \ Z(Pi)` is contained in the union of the family.
  5. A nonempty open `O ⊂ R_+` is not contained in `Z(Pi)` (Lemma Z), so it meets the union. ∎

### 3.6 Relation to the general theory (CITED, for context only; nothing above depends on it)

- G. E. Collins, *Quantifier elimination for real closed fields by cylindrical algebraic decomposition*, in: Automata
  Theory and Formal Languages, 2nd GI Conference, LNCS 33, Springer 1975, pp. 134-183.
  - Collins's projection also contains the lower coefficients and the principal subresultant coefficients. That is
    needed to stay delineable over base cells where leading coefficients vanish, and to build sign-invariant cells
    of every dimension.
- S. McCallum, *An improved projection operation for cylindrical algebraic decomposition of three-dimensional
  space*, J. Symbolic Comput. 5 (1988) 141-161.
- S. McCallum, *An improved projection operation for cylindrical algebraic decomposition*, in: B. F. Caviness,
  J. R. Johnson (eds.), Quantifier Elimination and Cylindrical Algebraic Decomposition, Springer 1998, pp. 242-268.
  - McCallum's projection for a squarefree, pairwise coprime basis is coefficients + discriminants + resultants,
    under a well-orientedness hypothesis.
- C. W. Brown, *Improved projection for cylindrical algebraic decomposition*, J. Symbolic Comput. 32 (2001) 447-465.
  - Brown replaces "all coefficients" by "leading coefficients" plus point checks. The set `P` used here is that
    projection in the plane.
- S. Basu, R. Pollack, M.-F. Roy, *Algorithms in Real Algebraic Geometry*, 2nd ed., Springer 2006, Ch. 4 (resultants,
  subresultants) and Ch. 5 (CAD).

All bibliographic details are from memory; verify them before citing externally.

The elementary proof above covers exactly what `region2d.py` needs:
- **full-dimensional** sectors over **open** `beta`-intervals on which every leading coefficient is nonzero;
- plus density.

It makes **no claim** about:
- sections (the graphs `Gamma_j`) or vertical lines, i.e. cells of dimension 0 and 1;
- sign-invariance of polynomials outside `F ∪ V ∪ C`, for example `Delta_k` with `k <= m-2` and the `c_l` with
  `l <= m-1`. Those may vanish inside sectors, and B'' (3) explains why that is harmless.

None of this is needed, because B'' (4)-(5) uses only full-dimensional cells.

---

## 4. What `region2d.py` guarantees (code as of commit `b5beb30`)

**Scope.** This section covers `region2d.region2d(A, sigma)` for one `(graph, support)`. It does not cover which
graphs are passed to it: that selection is made by `run_region_sweep.py` / `run_region2d_sweep.py` and is outside
this note.

### 4.1 Hypotheses of Theorem B''', and where the code meets them

| hypothesis | where in `region2d.py` | how it is ensured |
|---|---|---|
| inputs `a_i, o_k, d, s, phi, omega` exact | `fastregion.beta_data` | Lemma A degree bounds (all `<= m` in `beta`); integer interpolation at `m+1` points; one extra check point |
| `D = c_m`, `Delta_(m-1)` as defined in §0 | `curves_lc` (l. 52) | `c[k] = phi_k + T omega_k`, Hurwitz entry `c[2i-j]`, `c[0] = 1`. Matches `lemmas.md` §0. **TRUSTED:** sympy Berkowitz determinant. SANITY T1 |
| `F` = curves | `projection` (l. 101) + `_coprime` (l. 74) | Curves of `deg_t >= 2` (in practice only from `Delta_(m-1)`) are split by `sp.factor_list` (l. 109). `D` and `L` (`deg_t <= 1`) are kept whole. `_coprime` replaces any pair with a nontrivial gcd `g` by `(f_i/g, f_j/g, g)` until all pairs are coprime |
| `V` = vertical polynomials | `projection`, `_coprime` | factors / gcd pieces with `deg_t = 0` go to `polys` (l. 89, 112); the nonconstant `a_i`, `o_k` are `pos_polys`; both enter `intp` (l. 291) |
| (P≠0) | l. 124, l. 130 | `lc` nonzero by definition. `disc` nonzero is **asserted** (l. 124). `res` nonzero is **checked**: it raises `ArithmeticError` (l. 130). `V` nonzero by Lemma N and factorisation |
| covering condition for `D`, `Delta_(m-1)`, `L` | `projection`, `_coprime` | **By construction**, with exact `factor_list`, `gcd`, `cancel`: each of `D`, `Delta_(m-1)`, `L` equals a nonzero constant times a product of powers of members of `F` and `V`. Factorisation is an identity, and each `_coprime` step is an exact division that drops only nonzero constants. A curve such as `D` may carry vertical content `g(beta)`; its zeros are roots of `lc_t D = s`, which is in `P`. **TRUSTED:** sympy. SANITY COV |
| `F` contains `c L` | `curves_lc`, `_coprime` | `L = (beta-2) t - 1` is irreducible, so every gcd with `L` is 1 or `L` up to a constant, and some member of `F` is always `c L` |
| `I_alpha` = components of `B+ \ roots(P)` | `positivity_intervals` (l. 256); main loop l. 294-329 | `B+` = the cells of `(2, inf)` between consecutive roots of the nonconstant `a_i`, `o_k` whose sample passes (A.2; signs are constant between roots). Inside each `B+`-cell, joint isolation of `intp` (l. 297), **replayed**. Beta-cells lie between consecutive isolating intervals; the sample `b0` is the midpoint of a rational inner interval, so it lies in the true `I_alpha` |
| slices at `b0` (B''' (g)) | l. 330-342 | Curves are specialised at `b0` (each has full `t`-degree there, since `lc(b0) != 0`) and isolated on `[t_min, inf)`, **replayed**. The `L`-root `t_min` is dropped as a degenerate interval (l. 341); it is the only possible root at `t_min`, since `res(f, L)(b0) != 0`. `fastregion._cells` builds the open intervals between consecutive isolating intervals. The k-th t-cell is the slice of `S_(j_L + k)` |
| exact verdict at a sample | `fastregion.decide_point`, `ctln.hurwitz_stable` | (A.1) with exact `Fraction`s. `D0 = 0` (`positive = None`) cannot occur in a sector. Cross-checked against `ctln.analyze_support` when `crosscheck=True` (the sweep default) |

**Two details.**
- **The `L`-root is a degenerate interval.** Suppose sympy's intervals lie in `[inf, sup] = [t_min, inf)`. This is
  what its `inf`/`sup` arguments are for; I did not verify it in sympy's source. Then a non-degenerate interval
  `[t_min, c]` for `L` would fail the replay's strict sign-change test, because `L(t_min) = 0`. So the interval
  filtered at l. 341 is the degenerate one.
  - If that assumption ever failed, the interval would survive the filter. `fastregion._cells` would then keep
    refining an interval whose root is exactly `t_min`, which is a hang, not a wrong verdict.
  - SANITY: the interval was degenerate in 1060 of 1060 specialisations, and the n=5 and n=6 sweeps completed.
- **Constant-in-`t` "curves".** A curve that is constant in `t` (for example `Delta_0 = 1` when `m = 1`) can remain in
  the `curves` list. It is harmless: its `lc` is itself, so it is in `P` if nonconstant, and it has no `t`-roots.

### 4.2 Conclusions for the three verdicts

Assume the trust base below and that every runtime check passed. Then the family of (beta-cell, t-cell) pairs that
`region2d.region2d` samples is exactly the family of Theorem B''' for `F` = the `curves` of positive `t`-degree and
`V` = the vertical polynomials together with the nonconstant `a_i`, `o_k`. Each sector is sampled exactly once
(re-samples stay in the same sector). Hence:

- **`"never"` (positivity set empty).** `S = ∅` by (A.2) alone.
- **`"never"` (cells sampled, none stable, none undecided).** `S = ∅` by Corollary B'' (4).
- **`"stable_somewhere"`.**
  - `S != ∅` is certified by the exact pointwise verdict at a rational sample. This uses no CAD argument at all
    (cf. `notes/lemma_F_reduction.md` Lemma O, consequence (O1)).
  - The CAD adds that each reported stable sector is entirely stable (B'' (2)).
  - It also adds the exact description `S = (R ∩ cl(S_c)) \ Z_0`, where `S_c` is the union of the stable sectors
    (B'' (5)).
- **`"escalate"`.**
  - A cell whose three samples are all undecided is in fact **not stable**, by B'' (3)(b). So an `escalate` verdict
    with no stable cell is mathematically `never`.
  - The code stays conservative. Changing that is a policy decision (§3 of CLAUDE.md is about points, not cells).
  - The existing tallies have `escalate = 0` at n=5 and n=6, so no current result depends on this.

### 4.3 Trust base (TRUSTED, not re-verified at runtime by `region2d.py`)

| item | role | sanity |
|---|---|---|
| T1 | sympy `Matrix.det(method="berkowitz")`, for `Delta_(m-1)` in `Z[beta, t]` | SANITY, 35 cases |
| T2 | sympy `factor_list`; only the identity `Delta_(m-1) = c prod g_i^(e_i)` matters | SANITY COV, 35 cases |
| T3 | sympy `gcd` / `cancel` in `_coprime`; only exactness of the division matters | SANITY COV, 35 cases |
| T4 | sympy `resultant` and `discriminant`. **Delineability rests on these**, and only their zero sets matter | SANITY, 35 cases |
| T5 | sympy `dup_sqf_part` inside the Descartes replay, and the replay code itself, which I read but did not audit line by line | not checked here |
| T6 | Python `Fraction` and the `fastregion` integer routines (`idet`, `icharpoly`, `cramer_ones`, `interp`) | used project-wide and covered by `test_ctln.py`; not re-audited here |

T1-T4 were re-verified independently on 35 cases, by Bareiss determinants at enough integer points to force identity
of polynomials (Appendix). That is evidence, not a proof for all 34,071 n=6 records (28,350 of them with a nonempty CAD).

**Hardening suggestions (not implemented; this note modifies no code).**
- (h1) Re-derive `Delta_(m-1)` by `fastregion.hurwitz_det_poly` at `m` values of `t` and compare.
- (h2) Multiply factorisations and gcd splits back and compare.
- (h3) Recompute every resultant and discriminant by integer Sylvester determinants at `deg-bound + 1` points, as in
  `check_projection.py`. Compare up to one global sign.
- (h4) Enforce the B'' (3)(c) consistency check: an undecided first sample followed by a True re-sample must stop the
  run.

### 4.4 Findings from this audit (no soundness defect found)

- **F1, reporting caveat.** The `stable_cells` tuples `(beta_lo, beta_hi, beta_sample, t_lo, t_hi)` are not
  rectangles.
  - `beta_lo` and `beta_hi` are rational *inner* bounds of the true beta-cell, whose endpoints are algebraic and listed
    in `beta_critical`.
  - `t_lo` and `t_hi` are bounds of the t-slice **at `beta = beta_sample` only**. The sector's `t`-range moves with
    `beta`.
  - Example: `notes/n6_map.md` line 20 reports `['2', '17/8', '33/16', '16', None]`. Read as a box, that contains the
    illegal point `(beta, t) = (2.01, 17)`, where `L < 0`.
  - The correct reading: `16 = 1/(33/16 - 2)` is just `t_min` at the sample. So this stable sector is the whole legal
    region over the beta-cell, `{ 2 < beta < rho_1, t > 1/(beta-2) }`. Here `rho_1` is the algebraic number isolated
    in `[17/8, 15/7]` by the degree-14 polynomial in `beta_critical`.
  - The second tuple on the same line, `['15/7', '1029/475', '7164/3325', '29', None]`, is different. There `29` is a
    rational inner bound, at the sample only, for a section of a curve above `t_min ≈ 6.47`.
  - `aggregate.py` l. 169 prints these tuples verbatim. No soundness consequence; documentation only.
- **F2, sympy resultant sign.** sympy 1.14's `resultant(f, g)` equals `-det Syl(f, g)` when `deg f < deg g` and
  `deg f * deg g` is odd.
  - Example: `resultant(t - 1, t^3 - 8, t)` returns `7`, while `Res = g(1) = -7`.
  - Only zero sets are used, so `region2d.py` is unaffected.
  - Any future use of a resultant's **sign** (e.g. a Sturm-Habicht or signed-subresultant argument) must not rely on
    sympy's sign. Suggested ANOMALIES entry; I have not edited that file.
- **F3, the zero-polynomial case.** Lemma N shows that `a_i`, `o_k`, `D` and every `Delta_k` are never the zero
  polynomial, in two variables and also for fixed `t >= 0`. The zero-polynomial branch of `lemmas.md` Lemma B(d) is
  therefore vacuous, and `region2d.py` never needs to handle it.
- **F4, escalation is conservative.** See 4.2.

---

## 5. Not proved here / open

- **Correctness of sympy's exact algebra (T1-T5) is TRUSTED.** It was sanity-checked on 35 cases, not proved for the
  full sweep. Implementing (h1)-(h3) would remove T1-T4 from the trust base at modest cost.
- **The Descartes-bisection replay (T5) was read but not audited line by line.** In particular, I did not prove
  termination of `_count_open01` or the exactness of `count_roots_closed`'s change of variable. Its role is to
  confirm sympy's isolation; a bug there could only make the replay *fail to catch* an isolation error, not create a
  wrong verdict by itself.
- **The Routh-Hurwitz theorem (b0) is CITED** (`lemmas.md` B.1). Everything else in §0-§3 is proved here or in the
  cited project notes.
- **Cells of dimension 0 and 1 are not analysed** (no sign-invariance claim on sections or vertical lines). By B'' (5)
  this is not needed for any statement about `S`.
- **The Lean formalisation of B'/B''/B''' has not been attempted.** For a Lean port:
  - (R)/R2 would come from Mathlib's continuity of roots, if available; I have not checked what Mathlib provides.
  - The Sylvester kernel lemma, Conn, Z and the sector construction are elementary.

---

## Appendix: sanity checks run (tests, not proofs)

All scripts are in `/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverB/`. They
use exact integer / `Fraction` arithmetic only, and import the project modules read-only. JSON outputs sit next to
the scripts.

**Cases (35).**
- Geneson's graph, full support.
- The 4 `stable_somewhere` records of the n=6 region2d sweep.
- The 6 records with the most beta-cells (22-23).
- 14 random n=6 records with a nonempty CAD.
- 10 random proper supports (sizes 2-5) of random n=6 graphs.

**`check_projection.py`** prints `ALL PROJECTION CHECKS PASSED on 35 cases`.
- **N:** `a_i(0) = o_k(0) = D(0,0) = 1` and `Delta_(m-1)(0,0) > 0`.
- **H:** region2d's `Delta_(m-1)` equals the Bareiss determinant of the Hurwitz matrix built from
  `c_l = phi_l + t omega_l`, on an integer grid of size `(deg bound + 1)` in each variable. That is an identity of
  polynomials.
- **RES:** each of the 186 projection resultants equals the Bareiss determinant of the Sylvester matrix (formal
  degrees) at `deg-bound + 1` integer `beta` points, up to one global sign per pair.
  - The sign is `-1` exactly for `(deg_t f, deg_t g) in {(1,3), (1,5)}`: 14 pairs, all with `p < q` and `pq` odd.
  - It is `+1` in all 172 other pairs (finding F2).
- **DSC:** `res_t(f, ∂_t f) = (-1)^(d(d-1)/2) lc_t f disc_t f` for every curve of `t`-degree `>= 2`.
- **COV:** each of `D`, `Delta_(m-1)` and `L` reduces to a nonzero constant after exact division by the curves and
  the projection polynomials.

**`sector_constancy.py`** prints `ALL SECTOR CHECKS PASSED on 35 cases`.
- Runs `region2d.region2d` on each case. In every beta-cell it takes 4 random rational `beta` and recomputes the
  t-sections exactly.
- **DEL:** the number of sections above `t_min`, and the curve owning each section, are identical to those at the
  sample `b0`. This is delineability evidence.
- **SEC:** at a random rational `t` strictly inside each slice, the exact verdict equals the cell verdict. 1012
  points, 0 undecided, 0 disagreements.
- **LROOT:** the `L`-root is a degenerate isolating interval in 1060 of 1060 specialisations.

**`undecided_inside_cells.py`** prints `B''(3) SANITY PASSED`.
- `c_l = phi_l + t omega_l` is linear in `t`. So at each beta-cell sample `b0` and each `l = 1..m-1`, the point
  `t* = -phi_l(b0)/omega_l(b0)` lies exactly on `{c_l = 0}`.
- 77 such points were legal and strictly inside a t-cell. All 77 are undecided by `ctln.hurwitz_stable`, as they
  must be, and none lies in a stable cell, as B'' (3)(d) requires.

Re-run with:

    cd /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverB
    python3 check_projection.py; python3 sector_constancy.py; python3 undecided_inside_cells.py

They take a few minutes in total.
