/-
  Remark F (notes/lemmas.md section 7) and its spectral consequences.

  * `affine_identity`      K(l beta, l t) = (1 - l) I + l K(beta, t)       (any commutative ring)
  * `charpoly_affine_eval` det(x I - (c I + l M)) = l^m det(((x - c)/l) I - M)   (any field, l != 0)
  * `charpoly_affine`      the same as an identity of polynomials (infinite field)
  * `charpoly_ray`         specialisation to K(l beta, l t)
  * `eig_affine`           mu eigenvalue of c I + l M  <->  (mu - c)/l eigenvalue of M
  * `eig_ray_decomp`       every eigenvalue of K(l beta, l t) is (1 - l) + l nu, nu an eigenvalue of K(beta, t)
  * `posStable_ray`        K(beta, t) positive stable and 0 < l <= 1  ==>  K(l beta, l t) positive stable
-/
import CtlnCert.Basic

open Matrix Polynomial

namespace CtlnCert

variable {n : Type*} [Fintype n] [DecidableEq n]

omit [Fintype n] in
/-- **Remark F, affine identity.**  Valid for every `l`, over any commutative ring. -/
theorem affine_identity {R : Type*} [CommRing R] (N : Matrix n n R) (beta t l : R) :
    K N (l * beta) (l * t) = (1 - l) • (1 : Matrix n n R) + l • K N beta t := by
  unfold K
  rw [mul_smul, mul_smul, smul_add, smul_add, sub_smul, one_smul]
  abel

/-- Characteristic polynomial of `c I + l M`, evaluation form: for `l ≠ 0` and every `x`,
`charpoly(c I + l M)(x) = l^m * charpoly(M)((x - c)/l)`. -/
theorem charpoly_affine_eval {F : Type*} [Field F] (M : Matrix n n F) (c l x : F) (hl : l ≠ 0) :
    (c • (1 : Matrix n n F) + l • M).charpoly.eval x
      = l ^ Fintype.card n * M.charpoly.eval ((x - c) / l) := by
  rw [eval_charpoly, eval_charpoly, ← det_smul]
  congr 1
  ext i j
  by_cases hij : i = j
  · subst hij
    simp [scalar_apply, mul_sub]
    field_simp
    ring
  · simp [scalar_apply, one_apply, hij]

/-- Characteristic polynomial of `c I + l M` as a polynomial identity (infinite field, `l ≠ 0`):
`charpoly(c I + l M) = C(l^m) * charpoly(M) ∘ (C l⁻¹ * (X - C c))`. -/
theorem charpoly_affine {F : Type*} [Field F] [Infinite F] (M : Matrix n n F) (c l : F)
    (hl : l ≠ 0) :
    (c • (1 : Matrix n n F) + l • M).charpoly
      = C (l ^ Fintype.card n) * M.charpoly.comp (C l⁻¹ * (X - C c)) := by
  apply Polynomial.funext
  intro x
  rw [charpoly_affine_eval M c l x hl]
  simp [eval_comp, div_eq_inv_mul]

/-- The ray specialisation: for `l ≠ 0`,
`charpoly(K(l beta, l t))(x) = l^m * charpoly(K(beta, t))((x - (1 - l))/l)`. -/
theorem charpoly_ray {F : Type*} [Field F] (N : Matrix n n F) (beta t l x : F) (hl : l ≠ 0) :
    (K N (l * beta) (l * t)).charpoly.eval x
      = l ^ Fintype.card n * (K N beta t).charpoly.eval ((x - (1 - l)) / l) := by
  rw [affine_identity]
  exact charpoly_affine_eval _ _ _ _ hl

/-- Eigenvalues of `c I + l M` (`l ≠ 0`), eigenvector form. -/
theorem eig_affine (M : Matrix n n ℂ) (c l mu : ℂ) (hl : l ≠ 0) :
    Eig (c • (1 : Matrix n n ℂ) + l • M) mu ↔ Eig M ((mu - c) / l) := by
  constructor
  · rintro ⟨v, hv, h⟩
    refine ⟨v, hv, ?_⟩
    rw [add_mulVec, smul_mulVec, smul_mulVec, one_mulVec] at h
    have h2 : l • (M *ᵥ v) = (mu - c) • v := by
      rw [sub_smul, ← h]; abel
    calc M *ᵥ v = l⁻¹ • (l • (M *ᵥ v)) := by rw [smul_smul, inv_mul_cancel₀ hl, one_smul]
      _ = ((mu - c) / l) • v := by rw [h2, smul_smul, div_eq_inv_mul]
  · rintro ⟨v, hv, h⟩
    refine ⟨v, hv, ?_⟩
    rw [add_mulVec, smul_mulVec, smul_mulVec, one_mulVec, h, smul_smul, mul_div_cancel₀ _ hl,
      sub_smul]
    abel

omit [Fintype n] in
/-- Complexification commutes with the affine combination. -/
theorem cplx_affine (M : Matrix n n ℝ) (c l : ℝ) :
    cplx (c • (1 : Matrix n n ℝ) + l • M) = (c : ℂ) • (1 : Matrix n n ℂ) + (l : ℂ) • cplx M := by
  ext i j
  by_cases hij : i = j
  · subst hij; simp [cplx]
  · simp [cplx, one_apply, hij]

/-- **Spectral consequence of Remark F.**  For `l ≠ 0`, every eigenvalue `mu` of `K(l beta, l t)`
has the form `(1 - l) + l nu` with `nu` an eigenvalue of `K(beta, t)` (eigenvalues over ℂ). -/
theorem eig_ray_decomp (N : Matrix n n ℝ) (beta t l : ℝ) (hl : l ≠ 0) (mu : ℂ)
    (hmu : Eig (cplx (K N (l * beta) (l * t))) mu) :
    ∃ nu : ℂ, Eig (cplx (K N beta t)) nu ∧ mu = (1 - l : ℂ) + (l : ℂ) * nu := by
  rw [affine_identity, cplx_affine] at hmu
  have hl' : (l : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hl
  refine ⟨(mu - ((1 - l : ℝ) : ℂ)) / l, (eig_affine _ _ _ _ hl').mp hmu, ?_⟩
  push_cast
  field_simp
  ring

/-- **Stability propagates down rays.**  If `K(beta, t)` is positive stable and `0 < l ≤ 1`,
then `K(l beta, l t)` is positive stable. -/
theorem posStable_ray (N : Matrix n n ℝ) (beta t l : ℝ) (h : PosStable (K N beta t))
    (hl0 : 0 < l) (hl1 : l ≤ 1) : PosStable (K N (l * beta) (l * t)) := by
  intro mu hmu
  obtain ⟨nu, hnu, rfl⟩ := eig_ray_decomp N beta t l hl0.ne' mu hmu
  have hre : 0 < nu.re := h nu hnu
  have : ((1 - l : ℂ) + (l : ℂ) * nu).re = (1 - l) + l * nu.re := by
    simp [Complex.add_re, Complex.mul_re]
  rw [this]
  nlinarith

/-- The same statement for an arbitrary real matrix in place of `K(beta, t)`:
`(1 - l) I + l M` is positive stable whenever `M` is and `0 < l ≤ 1`. -/
theorem posStable_convex (M : Matrix n n ℝ) (h : PosStable M) (l : ℝ) (hl0 : 0 < l)
    (hl1 : l ≤ 1) : PosStable ((1 - l) • (1 : Matrix n n ℝ) + l • M) := by
  intro mu hmu
  rw [cplx_affine] at hmu
  have hl' : (l : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hl0.ne'
  have hnu := h _ ((eig_affine _ _ _ _ hl').mp hmu)
  set nu := (mu - ((1 - l : ℝ) : ℂ)) / (l : ℂ) with hnudef
  have hmu_eq : mu = ((1 - l : ℝ) : ℂ) + (l : ℂ) * nu := by
    rw [hnudef]; field_simp; ring
  rw [hmu_eq]
  have : (((1 - l : ℝ) : ℂ) + (l : ℂ) * nu).re = (1 - l) + l * nu.re := by
    simp [Complex.add_re, Complex.mul_re]
  rw [this]
  nlinarith

end CtlnCert
