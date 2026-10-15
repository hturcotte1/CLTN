# Lean formalisation: Remark F, Lemma R (forward), Lemma A.2, the e_2 edge bound, Lemma D

Source mathematics: notes/lemmas.md section 7 (Remark F), notes/lemma_M_refuted.md section 3.1 (Lemma R),
ANOMALIES.md entry "2026-10-14 Lemma R". Conventions: CLAUDE.md section 2 and notes/lemmas.md section 0.
Nothing in the repository was modified; everything lives under this directory.

## Toolchain and Mathlib (what worked)

- elan 4.2.4 (pre-installed). Installed `leanprover/lean4:v4.35.0` (the toolchain of Mathlib's `v4.35.0` tag
  and of Mathlib master on 2026-10-15), set as elan default.
- Lake project `CtlnCert/` requires Mathlib at tag `v4.35.0` (commit 56a005855ea17e0b0958f73ec801ae140b0a1679).
- `lake update` fetched Mathlib and its dependencies and its post-update hook downloaded the prebuilt cache from
  https://cache.mathlib.org: 9031 .olean archives, about 95 s. A following `lake exe cache get` reported
  "No files to download". Mathlib was NOT built from source. `.lake` is about 7.3 GB.

## Files

| file | content |
|---|---|
| `CtlnCert/lakefile.toml`, `lean-toolchain`, `lake-manifest.json` | project, pinned toolchain and Mathlib |
| `CtlnCert/CtlnCert/Basic.lean` | `J`, `K N beta t = 1 + beta N + t J`, `Eig` (eigenvector over C), `cplx`, `PosStable`, `HurwitzStable`; `eig_iff_isRoot_charpoly`, `posStable_iff_charpoly` |
| `CtlnCert/CtlnCert/RayIdentity.lean` | Remark F: `affine_identity`; charpoly scaling `charpoly_affine_eval`, `charpoly_affine`, `charpoly_ray`; `eig_affine`, `cplx_affine`, `eig_ray_decomp`; `posStable_ray`, `posStable_convex` |
| `CtlnCert/CtlnCert/IVTStep.lean` | `fR`, `ivt_step` (Mathlib IVT), `curve_point_of_root`; explicit root `lamStar`, `explicit_root_abstract` (abstract sqrt), `explicit_root` (Real.sqrt), `root_unique` |
| `CtlnCert/CtlnCert/Bridge.lean` | `W` exactly as CLAUDE.md section 2 (reads `A j i`), `Nmat A = J - 1 - Aᵀ`; `closed_form`, `IminusW_eq`, `legal_iff`, `hurwitz_iff_posStable` |
| `CtlnCert/CtlnCert/LemmaR.lean` | `lemmaR_forward` (K level), `lemmaR_forward_ctln` (eps, delta level), `curve_shift_legal` |
| `CtlnCert/CtlnCert/LemmaA2.lean` | Lemma A.2 (`det_pos_of_posStable`), Cramer/adjugate consequences, spectral half of the e_2 bound (`sum_minors_two_pos`) |
| `CtlnCert/CtlnCert/EdgeBound.lean` | `e2`, pair counts `nB nN nO arcs`, `e2_closed_form`, `bidir_gt_nonadj`, `arcs_lower_bound`, `edge_bound_posStable`, `edge_bound_ctln` |
| `CtlnCert/CtlnCert/Domination.lean` | `Dominates` (Curto Def 1.1), `comparison` (D.1), `lemmaD`, `lemmaD_full_support` (D.3), `lemmaD_induced` (D.4) |
| `CtlnCert/CtlnCert.lean` | root module importing all of the above |
| `CtlnCert/Audit.lean` | `#print axioms` for all 57 theorems (26 batch 1 + 31 batch 2) |
| `CtlnCert/build.log`, `CtlnCert/audit.log` | output of the last clean build and of the axiom audit |

## Commands (re-run)

```
export PATH=$HOME/.elan/bin:$PATH
cd <this directory>/CtlnCert
lake exe cache get          # no-op once the cache is present
lake build                  # builds CtlnCert (default target); "Build completed successfully"
lake env lean Audit.lean    # every theorem: [propext, Classical.choice, Quot.sound] (affine_identity, K_apply, minor_pair: [propext, Quot.sound]); no sorryAx
grep -rn "sorry\|admit\|^axiom" CtlnCert CtlnCert.lean   # empty
```

The only build warnings are Lake's advisory that Mathlib "is designed for use with the module system".

## Statements proved (no sorry)

- `affine_identity` (any commutative ring, any l): `K N (l*beta) (l*t) = (1 - l) • 1 + l • K N beta t`.
- `charpoly_affine_eval` (any field, l ≠ 0): `eval x (charpoly (c • 1 + l • M)) = l^m * eval ((x - c)/l) (charpoly M)`;
  `charpoly_affine` the same as a polynomial identity (infinite field); `charpoly_ray` its specialisation to K.
- `eig_ray_decomp` (l ≠ 0): every complex eigenvalue mu of K(l beta, l t) is `(1 - l) + l nu` for an eigenvalue nu of K(beta, t).
- `posStable_ray`: K(beta, t) positive stable and 0 < l ≤ 1 implies K(l beta, l t) positive stable.
- `ivt_step`: t > 0, t(beta - 2) > 1 implies some l in (2/beta, 1) with `l t (l beta - 2) - 1 = 0`;
  `curve_point_of_root`: then l t > 0, l beta > 2, (l t)(l beta - 2) = 1, l beta = 2 + 1/(l t).
- `explicit_root`: `lamStar = (t + sqrt(t^2 + beta t))/(beta t)` is such a root and lies in (2/beta, 1);
  `explicit_root_abstract` the same with any s > 0, s^2 = t^2 + beta t (no sqrt, no IVT);
  `root_unique`: f has only one positive root.
- `closed_form`: `1 - W = (-delta) • 1 + (1 + delta) • J - (eps + delta) • Aᵀ` (with the transpose).
- `IminusW_eq`: `1 - W = eps • K (J - 1 - Aᵀ) (1 + delta/eps) ((1 - eps)/eps)` (eps ≠ 0).
- `legal_iff` (0 < eps < 1): `delta(1-eps)/eps > 1 ↔ t (beta - 2) > 1`.
- `hurwitz_iff_posStable` (eps > 0): `-1 + W` Hurwitz stable ↔ `K` positive stable.
- `lemmaR_forward`: K(beta, t) positive stable, t > 0, t(beta - 2) > 1 implies some t' in (0, t) with
  2 + 1/t' < beta and K(2 + 1/t', t') positive stable.
- `lemmaR_forward_ctln`: if `-1 + W A eps delta` is Hurwitz stable at a legal point (0 < eps < 1, q > 1), then
  `-1 + W A eps' delta'` is Hurwitz stable for some eps < eps' < 1, delta' > 0 with q(eps', delta') = 1 exactly
  (eps' = 1/(1 + t'), delta' = 1/t').

## Batch 2 (2026-10-15): Lemma A.2, the e_2 edge bound, Lemma D -- statements proved (no sorry, no extra hypothesis)

Sources: notes/lemmas.md section 1 (Lemma A.2) and section 5 (Lemma D, D.1, D.3, D.4);
notes/lemma_M_refuted.md section 5 (e_2 closed form and the edge bound).

**Lemma A.2 (`LemmaA2.lean`).**
- `eval_zero_pos_of_no_nonneg_root`: real polynomial f, leadingCoeff f > 0, no root in [0, inf) => f(0) > 0 (IVT).
- `negK_root_re_neg`: K positive stable => every complex root of charpoly(-K) (= det(x I + K)) has Re < 0.
- `det_pos_of_posStable`: **PosStable K -> 0 < det K.** Proof differs from notes/lemmas.md (conjugate pairing):
  x |-> det(x I + K) is monic, has no root in [0, inf) (a root x >= 0 would make -x an eigenvalue of K), so it is
  positive at 0. Same statement.
- `adjugate_mulVec_eq`: K x = c -> adj(K) c = det K • x (any commutative ring, no invertibility).
- `pos_iff_adjugate_pos`: PosStable K, K x = 1 -> ((forall i, 0 < x i) <-> forall i, 0 < (adj K 1) i).
- `solution_eq_adjugate`: PosStable K, K x = 1 -> x = (det K)^-1 • adj(K) 1.
- `ctln_pos_iff_adjugate_pos`: IsAdj A, eps > 0, theta > 0, -I + W Hurwitz stable, (I - W) x = theta 1 ->
  (x > 0 <-> adj(K(J - I - A^T, 1 + delta/eps, (1-eps)/eps)) 1 > 0). This is "permission is t-free" in Lean.
- `iterate_derivative_rootSet_re_neg`: every root of every iterated derivative of charpoly(-K) has Re < 0
  (Mathlib's Gauss-Lucas `rootSet_derivative_subset_convexHull_rootSet`, iterated).
- `sum_minors_two_pos`: PosStable K, card >= 2 -> 0 < sum over 2-subsets s of det K[s,s].
  Uses Mathlib's `Matrix.charpoly_coeff_eq_sum_minors` (coefficient x^(m-2) of charpoly(-K) = sum of 2x2 principal
  minors of K), then the (m-2)-th derivative of charpoly(-K) is a real quadratic with positive leading coefficient
  and no root in [0, inf), so its value at 0, (m-2)! * e_2, is positive.

**e_2 edge bound (`EdgeBound.lean`).** Nodes carry a `LinearOrder`; unordered pairs are `pairs n = {(i,j) : i < j}`.
- `e2 M := sum over i<j of (M i i M j j - M i j M j i)`; `e2_eq_sum_minors`: equals Mathlib's powersetCard-2 sum.
- `e2_closed_form` (any nontrivial commutative ring R with decidable equality; A 0/1 entries, zero diagonal):
  `e2 (K (J - 1 - A^T) beta t) = B(1+2t) + Nn(1-beta)(1+beta+2t) + O(1 - t(beta-2))`.
- `pair_count`: B + Nn + O = #pairs; `two_mul_card_pairs`: 2 #pairs = m(m-1); `arcs_eq`: #arcs = 2B + O.
- `bidir_gt_nonadj`: beta > 2, t > 0, t(beta-2) > 1, e2 > 0 -> Nn < B.
- `arcs_lower_bound`: same hypotheses -> m(m-1) + 2 <= 2 #arcs (i.e. #arcs >= m(m-1)/2 + 1).
- `e2_pos_of_posStable`: PosStable K, m >= 2 -> 0 < e2 K. **The spectral step is proved, not assumed.**
- `edge_bound_posStable`: PosStable K(J - 1 - A^T, beta, t), t > 0, t(beta-2) > 1, m >= 2 -> Nn < B and
  m(m-1) + 2 <= 2 #arcs.
- `edge_bound_ctln`: IsAdj A, 0 < eps < 1, delta(1-eps)/eps > 1, m >= 2, -I + W A eps delta Hurwitz stable ->
  Nn < B and m(m-1) + 2 <= 2 #arcs. (m >= 2 is necessary: for m <= 1, B = Nn = 0.)

**Lemma D (`Domination.lean`)**, W from Bridge.lean (W i j reads A j i).
- `Dominates A k j := k != j and A j k = 1 and not A k j = 1 and forall i, i != j -> i != k -> A i j = 1 -> A i k = 1`
  (Curto 2510.05098 Def 1.1 as quoted in notes/lemmas.md section 5; "k does not send to j" is written `A k j != 1`,
  which for 0/1 entries is A k j = 0 and is exactly what W tests).
- `input A eps delta theta x l := theta l + (W x) l`; `IsFixedPoint ... x := forall l, x l = max (input ... l) 0`.
- `comparison` (D.1): Dominates A k j, eps + delta >= 0, theta j <= theta k, x >= 0 ->
  y_j - y_k <= (1 - eps) x_j - (1 + delta) x_k.
- `lemmaD`: Dominates A k j, eps > 0, delta >= 0, theta j <= theta k, x a fixed point of x = [W x + theta]_+ -> x j = 0.
- `lemmaD_full_support` (D.3): same hypotheses, no x > 0 with (I - W) x = theta (any theta vector with
  theta j <= theta k); `lemmaD_full_support_uniform`: with theta • 1, any real theta.
- `dominated_not_complete`: a graph with a dominated node is not complete bidirected.
- `W_submatrix`: W(A[f,f]) = W(A)[f,f] for injective f; `lemmaD_induced` (D.4): domination inside G[sigma] rules out
  x > 0 on sigma with (I - W_sigma) x = theta 1.
- As in notes/lemmas.md: eps < 1, legality, theta > 0 and 0/1 entries are NOT used.

## Not formalised

- The converse (<=) of Lemma R needs openness of positive stability (continuity of eigenvalues in the
  entries). Mathlib has no ready-made statement for that, so it is not formalised. Only the elementary half,
  `curve_shift_legal` (the point (2 + 1/t' + eta, t') is legal for eta > 0), is proved.
- Condition 1 is covered only through Lemma A.2's consequence (`ctln_pos_iff_adjugate_pos`: positivity of x is the
  sign of adj(K) 1) and the domination obstructions; condition 2 (off-support, o_k < 0 of Lemma A (v)) and the
  rank-one identities of Lemma A (ii)-(v) are not formalised yet. The full certificate lemma (three conditions ->
  stable fixed point of the threshold-linear ODE) is still open in Lean.
