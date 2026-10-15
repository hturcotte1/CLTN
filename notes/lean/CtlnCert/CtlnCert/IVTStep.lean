/-
  Lemma R, IVT step (notes/lemma_M_refuted.md section 3.1).

  For t > 0 with t (beta - 2) > 1 (this forces beta > 2), put
      f(l) := l t (l beta - 2) - 1.
  Then f(2/beta) = -1 < 0 and f(1) = t(beta - 2) - 1 > 0, so there is l* in (2/beta, 1)
  with f(l*) = 0, i.e. (l* t)(l* beta - 2) = 1 and l* beta > 2.

  Two proofs are given:
  * `ivt_step`      : Mathlib's intermediate value theorem on the continuous polynomial f;
  * `explicit_root` : the explicit root l* = (t + sqrt(t^2 + beta t)) / (beta t), verified
                      algebraically, with 2/beta < l* < 1 proved from real arithmetic and
                      sqrt facts only; plus `root_unique`: it is the only root in (0, inf).
-/
import Mathlib

namespace CtlnCert

/-- The function of the IVT step. -/
def fR (beta t l : ℝ) : ℝ := l * t * (l * beta - 2) - 1

lemma beta_gt_two {beta t : ℝ} (ht : 0 < t) (hleg : t * (beta - 2) > 1) : 2 < beta := by
  by_contra h
  push Not at h
  nlinarith

lemma fR_at_two_over_beta {beta t : ℝ} (hb : beta ≠ 0) : fR beta t (2 / beta) = -1 := by
  unfold fR
  field_simp
  ring

lemma fR_at_one (beta t : ℝ) : fR beta t 1 = t * (beta - 2) - 1 := by
  unfold fR; ring

/-- **IVT step of Lemma R** (Mathlib IVT). -/
theorem ivt_step {beta t : ℝ} (ht : 0 < t) (hleg : t * (beta - 2) > 1) :
    ∃ l ∈ Set.Ioo (2 / beta) 1, fR beta t l = 0 := by
  have hb : 2 < beta := beta_gt_two ht hleg
  have hb0 : beta ≠ 0 := by linarith
  have hlt : 2 / beta < 1 := by rw [div_lt_one (by linarith)]; exact hb
  have hcont : ContinuousOn (fR beta t) (Set.Icc (2 / beta) 1) := by
    unfold fR; fun_prop
  have h0 : (0 : ℝ) ∈ Set.Ioo (fR beta t (2 / beta)) (fR beta t 1) := by
    rw [fR_at_two_over_beta hb0, fR_at_one]
    constructor <;> linarith
  obtain ⟨l, hl, hfl⟩ := intermediate_value_Ioo hlt.le hcont h0
  exact ⟨l, hl, hfl⟩

/-- Consequences of a root in `(2/beta, 1)`: the scaled point lies on the legality curve. -/
theorem curve_point_of_root {beta t l : ℝ} (ht : 0 < t) (hleg : t * (beta - 2) > 1)
    (hl : l ∈ Set.Ioo (2 / beta) 1) (hf : fR beta t l = 0) :
    0 < l * t ∧ 2 < l * beta ∧ (l * t) * (l * beta - 2) = 1 ∧ l * beta = 2 + 1 / (l * t) := by
  have hb : 2 < beta := beta_gt_two ht hleg
  have hbpos : 0 < beta := by linarith
  obtain ⟨hl1, _⟩ := hl
  have hlb : 2 < l * beta := by
    have := (div_lt_iff₀ hbpos).mp hl1
    linarith
  have hlpos : 0 < l := lt_trans (div_pos two_pos hbpos) hl1
  have hlt : 0 < l * t := mul_pos hlpos ht
  have hcurve : (l * t) * (l * beta - 2) = 1 := by unfold fR at hf; linarith
  refine ⟨hlt, hlb, hcurve, ?_⟩
  field_simp
  linarith

/-! ### Explicit root (no IVT): the quadratic formula -/

/-- The explicit root `l* = (t + s)/(beta t)` where `s = sqrt(t^2 + beta t)`. -/
noncomputable def lamStar (beta t : ℝ) : ℝ := (t + Real.sqrt (t ^ 2 + beta * t)) / (beta * t)

/-- Algebraic core, with an abstract square root `s` (`s^2 = t^2 + beta t`, `s > 0`):
`(t + s)/(beta t)` is a root of `f`, and lies in `(2/beta, 1)`. No IVT and no `Real.sqrt` used. -/
theorem explicit_root_abstract {beta t s : ℝ} (ht : 0 < t) (hleg : t * (beta - 2) > 1)
    (hs : s ^ 2 = t ^ 2 + beta * t) (hs0 : 0 < s) :
    fR beta t ((t + s) / (beta * t)) = 0 ∧
      2 / beta < (t + s) / (beta * t) ∧ (t + s) / (beta * t) < 1 := by
  have hb : 2 < beta := beta_gt_two ht hleg
  have hbpos : 0 < beta := by linarith
  have hbt : 0 < beta * t := mul_pos hbpos ht
  refine ⟨?_, ?_, ?_⟩
  · unfold fR
    field_simp
    nlinarith [hs]
  · -- 2/beta < (t+s)/(beta t)  <=>  2 t < t + s  <=>  t < s
    have hts : t < s := by
      by_contra h; push Not at h
      nlinarith [mul_pos hbpos ht]
    rw [div_lt_div_iff₀ hbpos hbt]
    nlinarith
  · -- (t+s)/(beta t) < 1  <=>  s < (beta - 1) t  <=>  t^2 + beta t < (beta-1)^2 t^2 (both sides > 0)
    rw [div_lt_one hbt]
    have hpos : 0 < (beta - 1) * t := mul_pos (by linarith) ht
    have hsq : s ^ 2 < ((beta - 1) * t) ^ 2 := by
      rw [hs]
      -- (beta-1)^2 t^2 - t^2 - beta t = t (t beta (beta - 2) - beta) = beta t (t(beta-2) - 1) > 0
      nlinarith [mul_pos hbt (by linarith : (0:ℝ) < t * (beta - 2) - 1)]
    have : s < (beta - 1) * t := by
      by_contra h; push Not at h
      nlinarith
    linarith

/-- **Explicit root** with `Real.sqrt`: `l* = (t + sqrt(t^2 + beta t))/(beta t)` is a root of `f`
in `(2/beta, 1)`. -/
theorem explicit_root {beta t : ℝ} (ht : 0 < t) (hleg : t * (beta - 2) > 1) :
    fR beta t (lamStar beta t) = 0 ∧ lamStar beta t ∈ Set.Ioo (2 / beta) 1 := by
  have hb : 2 < beta := beta_gt_two ht hleg
  have hpos : 0 < t ^ 2 + beta * t := by nlinarith
  have hs := Real.sq_sqrt hpos.le
  have hs0 := Real.sqrt_pos.mpr hpos
  obtain ⟨h1, h2, h3⟩ := explicit_root_abstract ht hleg hs hs0
  exact ⟨h1, h2, h3⟩

/-- Uniqueness: `f` has exactly one positive root (beta, t > 0), so the IVT root and the
explicit root coincide. -/
theorem root_unique {beta t l₁ l₂ : ℝ} (hb : 0 < beta) (ht : 0 < t)
    (h1 : 0 < l₁) (h2 : 0 < l₂) (hf1 : fR beta t l₁ = 0) (hf2 : fR beta t l₂ = 0) : l₁ = l₂ := by
  unfold fR at hf1 hf2
  -- beta t (l1 - l2)(l1 + l2) = 2 t (l1 - l2); if l1 != l2 then beta (l1 + l2) = 2,
  -- and then 1 = l1 t (l1 beta - 2) = -l1 l2 beta t < 0, contradiction.
  by_contra hne
  have hsub : l₁ - l₂ ≠ 0 := sub_ne_zero.mpr hne
  have key : (l₁ - l₂) * (t * (beta * (l₁ + l₂) - 2)) = 0 := by linear_combination hf1 - hf2
  have hsum : t * (beta * (l₁ + l₂) - 2) = 0 := (mul_eq_zero.mp key).resolve_left hsub
  have hsum' : beta * (l₁ + l₂) = 2 := by
    have := (mul_eq_zero.mp hsum).resolve_left ht.ne'
    linarith
  have : l₁ * t * (l₁ * beta - 2) = -(l₁ * l₂ * beta * t) := by
    have : l₁ * beta - 2 = -(l₂ * beta) := by linarith
    rw [this]; ring
  nlinarith [mul_pos (mul_pos (mul_pos h1 h2) hb) ht]

end CtlnCert
