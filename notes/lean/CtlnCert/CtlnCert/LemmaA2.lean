/-
  Lemma A.2 of notes/lemmas.md (sign of D = det K) and the spectral half of the e_2 bound
  (notes/lemma_M_refuted.md section 5, "Stability needs e_2 > 0").

  * `eval_zero_pos_of_no_nonneg_root`  a real polynomial with positive leading coefficient and no
                                        root in [0, inf) is positive at 0 (intermediate value theorem)
  * `negK_root_re_neg`                  K positive stable => every complex root of charpoly(-K) has Re < 0
  * `det_pos_of_posStable`              Lemma A.2: K positive stable => det K > 0
  * `adjugate_mulVec_eq`                K x = c  =>  adj(K) c = det K * x   (Cramer, no invertibility)
  * `pos_iff_adjugate_pos`              K positive stable, K x = 1  =>  (x > 0  <->  adj(K) 1 > 0)
  * `ctln_pos_iff_adjugate_pos`         CTLN form: -I + W Hurwitz, (I - W) x = theta 1, theta > 0, eps > 0
                                        =>  (x > 0  <->  adj(K(N, beta, t)) 1 > 0)
  * `sum_minors_two_pos`                K positive stable, card >= 2  =>  the sum of the 2x2 principal
                                        minors of K is > 0  (Mathlib's charpoly_coeff_eq_sum_minors +
                                        Gauss-Lucas iterated m-2 times + the IVT lemma above)

  The proof of Lemma A.2 here is NOT the conjugate-pairing argument of notes/lemmas.md; it is the
  equivalent continuity argument: x |-> det(x I + K) is a monic real polynomial with no root in
  [0, inf) (a root x >= 0 would make -x <= 0 an eigenvalue of K), so it is positive at x = 0.
-/
import CtlnCert.Bridge

open Matrix Polynomial Filter

namespace CtlnCert

/-- A real polynomial with positive leading coefficient and no root in `[0, ∞)` is positive at 0. -/
theorem eval_zero_pos_of_no_nonneg_root (f : ℝ[X]) (hlc : 0 < f.leadingCoeff)
    (hroot : ∀ x : ℝ, 0 ≤ x → f.eval x ≠ 0) : 0 < f.eval 0 := by
  by_cases hdeg : 0 < f.degree
  · have hev : ∀ᶠ x in atTop, 0 < f.eval x :=
      (Polynomial.tendsto_atTop_of_leadingCoeff_nonneg f hdeg hlc.le).eventually
        (eventually_gt_atTop 0)
    obtain ⟨a, ha⟩ := eventually_atTop.mp hev
    set R := max a 0 with hR
    have hRpos : 0 < f.eval R := ha R (le_max_left _ _)
    have hR0 : 0 ≤ R := le_max_right _ _
    by_contra h0
    push Not at h0
    have hmem : (0 : ℝ) ∈ Set.Icc (f.eval 0) (f.eval R) := ⟨h0, hRpos.le⟩
    obtain ⟨c, hc, hfc⟩ :=
      intermediate_value_Icc hR0 (f.continuous.continuousOn) hmem
    exact hroot c hc.1 hfc
  · push Not at hdeg
    rw [eq_C_of_degree_le_zero hdeg] at hlc ⊢
    simpa using hlc

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- If `K` is positive stable, every complex root of the characteristic polynomial of `-K`
(that is, every root of `x ↦ det(x I + K)`) has negative real part. -/
theorem negK_root_re_neg (K : Matrix n n ℝ) (h : PosStable K) (z : ℂ)
    (hz : ((-K).charpoly.map Complex.ofRealHom).IsRoot z) : z.re < 0 := by
  have hc : (cplx (-K)).charpoly = (-K).charpoly.map Complex.ofRealHom := by
    rw [← charpoly_map]; rfl
  rw [← hc, ← eig_iff_isRoot_charpoly] at hz
  have hK : -K = (0 : ℝ) • (1 : Matrix n n ℝ) + (-1 : ℝ) • K := by
    rw [zero_smul, zero_add, neg_one_smul]
  have hne : ((-1 : ℝ) : ℂ) ≠ 0 := by norm_num
  rw [hK, cplx_affine, eig_affine _ _ _ _ hne] at hz
  have := h _ hz
  have hre : ((z - ((0 : ℝ) : ℂ)) / ((-1 : ℝ) : ℂ)).re = -z.re := by
    simp [div_neg]
  rw [hre] at this
  linarith

/-- Real roots of `x ↦ det(x I + K)` are negative when `K` is positive stable. -/
theorem negK_real_root_neg (K : Matrix n n ℝ) (h : PosStable K) (x : ℝ)
    (hx : (-K).charpoly.eval x = 0) : x < 0 := by
  have hz : ((-K).charpoly.map Complex.ofRealHom).IsRoot (x : ℂ) := by
    rw [IsRoot]
    show eval (Complex.ofRealHom x) _ = 0
    rw [eval_map_apply, hx, map_zero]
  have := negK_root_re_neg K h x hz
  simpa using this

/-- **Lemma A.2** (notes/lemmas.md): a positive stable real matrix has positive determinant. -/
theorem det_pos_of_posStable (K : Matrix n n ℝ) (h : PosStable K) : 0 < K.det := by
  have h0 : (-K).charpoly.eval 0 = K.det := by
    rw [eval_charpoly]; simp
  rw [← h0]
  apply eval_zero_pos_of_no_nonneg_root
  · rw [(charpoly_monic _).leadingCoeff]; exact one_pos
  · intro x hx hroot
    have := negK_real_root_neg K h x hroot
    linarith

/-- Cramer in adjugate form (no invertibility assumed): `K x = c` implies `adj(K) c = det K • x`. -/
theorem adjugate_mulVec_eq {R : Type*} [CommRing R] (K : Matrix n n R) (x c : n → R)
    (hx : K *ᵥ x = c) : adjugate K *ᵥ c = K.det • x := by
  rw [← hx, mulVec_mulVec, adjugate_mul, smul_mulVec, one_mulVec]

/-- Consequence of Lemma A.2: for positive stable `K` and `K x = 1`, `x > 0` entrywise iff
`adj(K) 1 > 0` entrywise. -/
theorem pos_iff_adjugate_pos (K : Matrix n n ℝ) (h : PosStable K) (x : n → ℝ)
    (hx : K *ᵥ x = fun _ => 1) :
    (∀ i, 0 < x i) ↔ ∀ i, 0 < (adjugate K *ᵥ fun _ => (1 : ℝ)) i := by
  have hd := det_pos_of_posStable K h
  have hadj := adjugate_mulVec_eq K x _ hx
  rw [hadj]
  simp only [Pi.smul_apply, smul_eq_mul]
  constructor
  · intro hp i; exact mul_pos hd (hp i)
  · intro hp i; exact pos_of_mul_pos_right (hp i) hd.le

/-- The solution `x` of `K x = 1` for positive stable `K` is `adj(K) 1 / det K`. -/
theorem solution_eq_adjugate (K : Matrix n n ℝ) (h : PosStable K) (x : n → ℝ)
    (hx : K *ᵥ x = fun _ => 1) :
    x = (K.det)⁻¹ • (adjugate K *ᵥ fun _ => (1 : ℝ)) := by
  have hd := det_pos_of_posStable K h
  rw [adjugate_mulVec_eq K x _ hx, smul_smul, inv_mul_cancel₀ hd.ne', one_smul]

/-- CTLN form of the consequence of Lemma A.2 (conditions 1 and 3 of CLAUDE.md section 2):
if `-I + W` is Hurwitz stable, `eps > 0`, `theta > 0` and `(I - W) x = theta 1`, then
`x > 0` entrywise iff `adj(K(N, beta, t)) 1 > 0` entrywise, with `N = J - I - Aᵀ`,
`beta = 1 + delta/eps`, `t = (1 - eps)/eps`. -/
theorem ctln_pos_iff_adjugate_pos (A : Matrix n n ℝ) (hA : IsAdj A) {eps delta theta : ℝ}
    (he : 0 < eps) (hth : 0 < theta) (hs : HurwitzStable (-1 + W A eps delta))
    (x : n → ℝ) (hx : (1 - W A eps delta) *ᵥ x = fun _ => theta) :
    (∀ i, 0 < x i) ↔
      ∀ i, 0 < (adjugate (K (Nmat A) (1 + delta / eps) ((1 - eps) / eps)) *ᵥ
        fun _ => (1 : ℝ)) i := by
  set Km := K (Nmat A) (1 + delta / eps) ((1 - eps) / eps)
  have hK : PosStable Km := (hurwitz_iff_posStable A hA delta he).mp hs
  have hd := det_pos_of_posStable Km hK
  have hKx : Km *ᵥ x = fun _ => theta / eps := by
    rw [IminusW_eq A hA eps delta he.ne', smul_mulVec] at hx
    funext i
    have := congrFun hx i
    simp only [Pi.smul_apply, smul_eq_mul] at this
    field_simp
    linarith
  have hadj := adjugate_mulVec_eq Km x _ hKx
  have hc : (fun _ => theta / eps : n → ℝ) = (theta / eps) • (fun _ => (1 : ℝ)) := by
    funext i; simp
  rw [hc, mulVec_smul] at hadj
  have hte : 0 < theta / eps := div_pos hth he
  constructor
  · intro hp i
    have := congrFun hadj i
    simp only [Pi.smul_apply, smul_eq_mul] at this
    have h2 : 0 < Km.det * x i := mul_pos hd (hp i)
    rw [← this] at h2
    exact pos_of_mul_pos_right h2 hte.le
  · intro hp i
    have := congrFun hadj i
    simp only [Pi.smul_apply, smul_eq_mul] at this
    have h2 : 0 < theta / eps * (adjugate Km *ᵥ fun _ => (1 : ℝ)) i := mul_pos hte (hp i)
    rw [this] at h2
    exact pos_of_mul_pos_right h2 hd.le

/-- Iterated Gauss-Lucas: every complex root of every iterated derivative of
`charpoly(-K)` has negative real part, when `K` is positive stable. -/
theorem iterate_derivative_rootSet_re_neg (K : Matrix n n ℝ) (h : PosStable K) (k : ℕ) :
    (derivative^[k] ((-K).charpoly.map Complex.ofRealHom)).rootSet ℂ ⊆ {z : ℂ | z.re < 0} := by
  induction k with
  | zero =>
    intro z hz
    rw [Function.iterate_zero_apply, mem_rootSet] at hz
    have hr : ((-K).charpoly.map Complex.ofRealHom).IsRoot z := by
      rw [IsRoot]; simpa [coe_aeval_eq_eval] using hz.2
    exact negK_root_re_neg K h z hr
  | succ k ih =>
    rw [Function.iterate_succ_apply']
    set P := derivative^[k] ((-K).charpoly.map Complex.ofRealHom)
    by_cases hdeg : 0 < P.degree
    · exact (rootSet_derivative_subset_convexHull_rootSet hdeg).trans
        (convexHull_min ih (convex_halfSpace_re_lt 0))
    · push Not at hdeg
      rw [eq_C_of_degree_le_zero hdeg, derivative_C, rootSet_zero]
      exact Set.empty_subset _

/-- Spectral half of the e_2 bound: for a positive stable real `K` with at least two rows, the sum
of the 2x2 principal minors of `K` is positive.  (It is the coefficient of `x^(m-2)` in
`det(x I + K)`, and that polynomial's `(m-2)`-th derivative is a real quadratic with positive
leading coefficient and, by Gauss-Lucas, no root in `[0, ∞)`.) -/
theorem sum_minors_two_pos (K : Matrix n n ℝ) (h : PosStable K) (hm : 2 ≤ Fintype.card n) :
    0 < ∑ s ∈ Finset.univ.powersetCard 2,
      (K.submatrix (Subtype.val : s → n) (Subtype.val : s → n)).det := by
  set m := Fintype.card n with hmdef
  set p := (-K).charpoly with hp
  -- the coefficient of x^(m-2) in charpoly(-K) is the sum of the 2x2 principal minors of K
  have hcoef : p.coeff (m - 2) = ∑ s ∈ Finset.univ.powersetCard 2,
      (K.submatrix (Subtype.val : s → n) (Subtype.val : s → n)).det := by
    rw [hp, charpoly_coeff_eq_sum_minors (-K) 2 hm]
    rw [neg_one_sq, one_mul]
    refine Finset.sum_congr rfl fun s hs => ?_
    have hcard : Fintype.card s = 2 := by
      rw [Fintype.card_coe]; exact (Finset.mem_powersetCard.mp hs).2
    rw [submatrix_neg, Pi.neg_apply, Pi.neg_apply, det_neg, hcard, neg_one_sq, one_mul]
  rw [← hcoef]
  set f := derivative^[m - 2] p with hf
  have hpdeg : p.natDegree = m := charpoly_natDegree_eq_dim _
  have hpmon : p.coeff m = 1 := by
    have := (charpoly_monic (-K)).coeff_natDegree
    rwa [← hp, hpdeg] at this
  -- f(0) = (m-2)! * coeff_{m-2} p
  have hf0 : f.eval 0 = ((m - 2).descFactorial (m - 2) : ℝ) * p.coeff (m - 2) := by
    rw [← coeff_zero_eq_eval_zero, hf, coeff_iterate_derivative, zero_add, nsmul_eq_mul]
  -- f has degree <= 2 and coeff 2 > 0
  have hfdeg : f.natDegree ≤ 2 := by
    have := natDegree_iterate_derivative p (m - 2)
    rw [hpdeg, ← hf] at this
    omega
  have hf2 : f.coeff 2 = ((m).descFactorial (m - 2) : ℝ) := by
    rw [hf, coeff_iterate_derivative, nsmul_eq_mul]
    have : 2 + (m - 2) = m := by omega
    rw [this, hpmon, mul_one]
  have hdf : (0 : ℝ) < ((m).descFactorial (m - 2) : ℝ) := by
    have : 0 < (m).descFactorial (m - 2) := Nat.descFactorial_pos.mpr (by omega)
    exact_mod_cast this
  have hf2ne : f.coeff 2 ≠ 0 := by rw [hf2]; exact hdf.ne'
  have hfnat : f.natDegree = 2 := le_antisymm hfdeg (le_natDegree_of_ne_zero hf2ne)
  have hlc : 0 < f.leadingCoeff := by
    rw [leadingCoeff, hfnat, hf2]; exact hdf
  have hfne : f ≠ 0 := by
    intro h0; rw [h0] at hfnat; simp at hfnat
  have hpos : 0 < f.eval 0 := by
    apply eval_zero_pos_of_no_nonneg_root f hlc
    intro x hx hroot
    have hmem : (x : ℂ) ∈ (derivative^[m - 2] ((-K).charpoly.map Complex.ofRealHom)).rootSet ℂ := by
      rw [iterate_derivative_map, mem_rootSet]
      refine ⟨?_, ?_⟩
      · exact (Polynomial.map_ne_zero_iff Complex.ofRealHom.injective).mpr hfne
      · rw [coe_aeval_eq_eval]
        show eval (Complex.ofRealHom x) _ = 0
        rw [eval_map_apply, ← hp, ← hf, hroot, map_zero]
    have := iterate_derivative_rootSet_re_neg K h (m - 2) hmem
    simp only [Set.mem_ofPred_eq, Complex.ofReal_re] at this
    linarith
  rw [hf0] at hpos
  have hfac : (0 : ℝ) < ((m - 2).descFactorial (m - 2) : ℝ) := by
    have : 0 < (m - 2).descFactorial (m - 2) := Nat.descFactorial_pos.mpr le_rfl
    exact_mod_cast this
  exact pos_of_mul_pos_right hpos hfac.le

end CtlnCert
