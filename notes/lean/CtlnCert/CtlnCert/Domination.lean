/-
  Lemma D of notes/lemmas.md section 5 (graphical domination kills every fixed point through the
  dominated node), with W exactly as in Bridge.lean / CLAUDE.md section 2 (W i j reads A j i).

  Definition (Curto, arXiv:2510.05098, Def. 1.1, as quoted in notes/lemmas.md section 5):
  node k graphically dominates node j (k ≠ j) iff
    (i)  for every i not in {j, k}, i -> j implies i -> k, and
    (ii) j -> k and k ↛ j.
  With A i j = 1 iff i -> j:  A j k = 1, A k j ≠ 1 (for a 0/1 matrix: A k j = 0), and A i j = 1 -> A i k = 1 for i ∉ {j, k}.

  * `comparison`        D.1: x ≥ 0, theta j ≤ theta k, eps + delta ≥ 0
                        =>  y_j(x) - y_k(x) ≤ (1 - eps) x_j - (1 + delta) x_k,
                        where y_l(x) = theta_l + sum_i W l i x_i
  * `lemmaD`            Lemma D: eps > 0, delta ≥ 0, theta j ≤ theta k: every fixed point x of
                        x = [W x + theta]_+ has x_j = 0 (equivalently, not x_j > 0)
  * `lemmaD_full_support`  Corollary D.3: no x > 0 with (I - W) x = theta (any theta vector with
                        theta j ≤ theta k; in particular uniform theta of any sign)
  * `lemmaD_full_support_uniform`  the same with theta • 1
  * `dominated_not_complete`  the graph is not complete bidirected (k ↛ j)
  * `W_submatrix`, `lemmaD_induced`  Corollary D.4: domination inside an induced subgraph G[sigma]
                        (condition (i) quantified over sigma only) rules out a positive solution of
                        condition 1 on sigma.
  Not used: eps < 1, legality q > 1, theta > 0, 0/1 entries of A.
-/
import CtlnCert.Bridge

open Matrix

namespace CtlnCert

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Curto's graphical domination (arXiv:2510.05098, Def. 1.1): `k` dominates `j`. -/
def Dominates (A : Matrix n n ℝ) (k j : n) : Prop :=
  k ≠ j ∧ A j k = 1 ∧ ¬ A k j = 1 ∧ ∀ i, i ≠ j → i ≠ k → A i j = 1 → A i k = 1

/-- Input to node `l` at state `x`: `y_l(x) = theta_l + sum_i W l i x_i`. -/
noncomputable def input (A : Matrix n n ℝ) (eps delta : ℝ) (theta x : n → ℝ) (l : n) : ℝ :=
  theta l + (W A eps delta *ᵥ x) l

/-- Fixed point of `dx/dt = -x + [W x + theta]_+`. -/
def IsFixedPoint (A : Matrix n n ℝ) (eps delta : ℝ) (theta x : n → ℝ) : Prop :=
  ∀ l, x l = max (input A eps delta theta x l) 0

/-- **D.1 (comparison inequality).** -/
theorem comparison (A : Matrix n n ℝ) {k j : n} (hd : Dominates A k j) {eps delta : ℝ}
    (hed : 0 ≤ eps + delta) (theta x : n → ℝ) (hth : theta j ≤ theta k) (hx : ∀ i, 0 ≤ x i) :
    input A eps delta theta x j - input A eps delta theta x k
      ≤ (1 - eps) * x j - (1 + delta) * x k := by
  obtain ⟨hkj, hjk, hkj0, hi⟩ := hd
  have hjk' : j ≠ k := Ne.symm hkj
  have hA01 : ¬ A k j = 1 := hkj0
  unfold input
  simp only [mulVec, dotProduct]
  have hsum : ∑ i, W A eps delta j i * x i - ∑ i, W A eps delta k i * x i
      ≤ ∑ i, ((if i = j then (1 - eps) * x j else 0) + (if i = k then -(1 + delta) * x k else 0)) := by
    rw [← Finset.sum_sub_distrib]
    apply Finset.sum_le_sum
    intro i _
    by_cases hij : i = j
    · subst hij
      simp [W, hjk', hkj, hjk]
      ring_nf; rfl
    · by_cases hik : i = k
      · subst hik
        simp [W, hjk', hA01, hij]
        ring_nf; rfl
      · simp only [hij, hik, ↓reduceIte, add_zero]
        simp only [W, of_apply, Ne.symm hij, Ne.symm hik, ↓reduceIte]
        by_cases hij1 : A i j = 1
        · simp [hij1, hi i hij hik hij1]
        · simp only [hij1, ↓reduceIte]
          have : (-1 - delta) - (if A i k = 1 then -1 + eps else -1 - delta) ≤ 0 := by
            split_ifs <;> linarith
          rw [← sub_mul]
          exact mul_nonpos_of_nonpos_of_nonneg this (hx i)
  rw [Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.sum_ite_eq'] at hsum
  simp only [Finset.mem_univ, ↓reduceIte] at hsum
  linarith

/-- **Lemma D** (notes/lemmas.md section 5). If `k` dominates `j`, `eps > 0`, `delta ≥ 0` and
`theta j ≤ theta k`, then every fixed point `x` of `x = [W x + theta]_+` has `x j = 0`. -/
theorem lemmaD (A : Matrix n n ℝ) {k j : n} (hd : Dominates A k j) {eps delta : ℝ}
    (he : 0 < eps) (hdel : 0 ≤ delta) (theta x : n → ℝ) (hth : theta j ≤ theta k)
    (hfp : IsFixedPoint A eps delta theta x) : x j = 0 := by
  have hx : ∀ i, 0 ≤ x i := fun i => by rw [hfp i]; exact le_max_right _ _
  by_contra hne
  have hxj : 0 < x j := lt_of_le_of_ne (hx j) (Ne.symm hne)
  have hyj : input A eps delta theta x j = x j := by
    have h := hfp j
    rcases le_total (input A eps delta theta x j) 0 with h0 | h0
    · rw [max_eq_right h0] at h; linarith
    · rw [max_eq_left h0] at h; exact h.symm
  have hed : 0 ≤ eps + delta := by linarith
  have hc := comparison A hd hed theta x hth hx
  have hxk := hx k
  have hdk : 0 ≤ delta * x k := mul_nonneg hdel hxk
  have hej : 0 < eps * x j := mul_pos he hxj
  -- y_k ≥ y_j - (1 - eps) x_j + (1 + delta) x_k = eps x_j + (1 + delta) x_k > 0
  have hyk_pos : 0 < input A eps delta theta x k := by nlinarith
  have hyk : x k = input A eps delta theta x k := by
    rw [hfp k, max_eq_left hyk_pos.le]
  -- then x_k ≥ eps x_j + (1 + delta) x_k, i.e. 0 ≥ eps x_j + delta x_k > 0
  nlinarith

/-- **Corollary D.3 (full-support prune).** If `k` dominates `j`, `eps > 0`, `delta ≥ 0`, then
there is no `x > 0` (entrywise) with `(I - W) x = theta`, for any input vector with
`theta j ≤ theta k`. -/
theorem lemmaD_full_support (A : Matrix n n ℝ) {k j : n} (hd : Dominates A k j) {eps delta : ℝ}
    (he : 0 < eps) (hdel : 0 ≤ delta) (theta x : n → ℝ) (hth : theta j ≤ theta k)
    (hpos : ∀ i, 0 < x i) (hsol : (1 - W A eps delta) *ᵥ x = theta) : False := by
  have hyx : ∀ l, input A eps delta theta x l = x l := by
    intro l
    have h := congrFun hsol l
    rw [sub_mulVec, one_mulVec, Pi.sub_apply] at h
    unfold input; linarith
  have hfp : IsFixedPoint A eps delta theta x := by
    intro l; rw [hyx l, max_eq_left (hpos l).le]
  have := lemmaD A hd he hdel theta x hth hfp
  exact (hpos j).ne' this

/-- Corollary D.3 with uniform input `theta • 1` (CLAUDE.md section 2, condition 1). -/
theorem lemmaD_full_support_uniform (A : Matrix n n ℝ) {k j : n} (hd : Dominates A k j)
    {eps delta : ℝ} (he : 0 < eps) (hdel : 0 ≤ delta) (theta : ℝ) (x : n → ℝ)
    (hpos : ∀ i, 0 < x i) (hsol : (1 - W A eps delta) *ᵥ x = fun _ => theta) : False :=
  lemmaD_full_support A hd he hdel (fun _ => theta) x le_rfl hpos hsol

omit [Fintype n] [DecidableEq n] in
/-- A graph with a dominated node is not complete bidirected. -/
theorem dominated_not_complete (A : Matrix n n ℝ) {k j : n} (hd : Dominates A k j) :
    ¬ (∀ a b, a ≠ b → A a b = 1) := by
  intro h
  exact hd.2.2.1 (h k j hd.1)

omit [Fintype n] in
/-- The weight matrix of an induced subgraph is the corresponding submatrix of `W`. -/
theorem W_submatrix {m : Type*} [DecidableEq m] (A : Matrix n n ℝ) (f : m → n)
    (hf : Function.Injective f) (eps delta : ℝ) :
    W (A.submatrix f f) eps delta = (W A eps delta).submatrix f f := by
  ext a b
  by_cases h : a = b
  · subst h; simp [W]
  · simp only [W, of_apply, submatrix_apply, h, hf.ne h, ↓reduceIte]
    congr 1

omit [Fintype n] in
/-- **Corollary D.4 (induced-subgraph prune).** If, inside the induced subgraph `G[sigma]`, the
node `k` dominates `j`, then condition 1 of CLAUDE.md section 2 fails for `sigma`: there is no
`x > 0` on `sigma` with `(I - W_sigma) x = theta 1`. -/
theorem lemmaD_induced (A : Matrix n n ℝ) (sigma : Finset n) {k j : sigma}
    (hd : Dominates (A.submatrix (Subtype.val : sigma → n) Subtype.val) k j) {eps delta : ℝ}
    (he : 0 < eps) (hdel : 0 ≤ delta) (theta : ℝ) (x : sigma → ℝ) (hpos : ∀ i, 0 < x i)
    (hsol : (1 - (W A eps delta).submatrix (Subtype.val : sigma → n) Subtype.val) *ᵥ x
      = fun _ => theta) : False := by
  rw [← W_submatrix A _ Subtype.val_injective] at hsol
  exact lemmaD_full_support_uniform _ hd he hdel theta x hpos hsol

end CtlnCert
