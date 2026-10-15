/-
  Lemma R (notes/lemma_M_refuted.md section 3.1), forward direction, fully formal.

  * `lemmaR_forward`       K-level: K(beta, t) positive stable with t > 0, t(beta - 2) > 1
                           ==> exists t' in (0, t) with K(2 + 1/t', t') positive stable and
                               2 + 1/t' < beta.
  * `lemmaR_forward_ctln`  CTLN level, CLAUDE.md section 2 conventions:
                           -I + W(eps, delta) Hurwitz stable at a legal point (0 < eps < 1,
                           delta > 0, q > 1) ==> -I + W(eps', delta') Hurwitz stable at a point
                           with q(eps', delta') = 1 exactly, eps < eps' < 1, delta' > 0.
  * `curve_shift_legal`    the elementary half of the converse: any point (2 + 1/t' + eta, t'),
                           eta > 0, t' > 0, is legal.  The converse itself also needs openness of
                           positive stability (continuity of eigenvalues), which is NOT formalised here.
-/
import CtlnCert.Bridge
import CtlnCert.IVTStep

open Matrix

namespace CtlnCert

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Lemma R, forward direction (K level).** -/
theorem lemmaR_forward (N : Matrix n n ℝ) {beta t : ℝ} (ht : 0 < t) (hleg : t * (beta - 2) > 1)
    (hs : PosStable (K N beta t)) :
    ∃ t' : ℝ, 0 < t' ∧ t' < t ∧ 2 + 1 / t' < beta ∧ PosStable (K N (2 + 1 / t') t') := by
  obtain ⟨l, hl, hf⟩ := ivt_step ht hleg
  obtain ⟨hlt, hlb, _, hcurve⟩ := curve_point_of_root ht hleg hl hf
  have hb : 2 < beta := beta_gt_two ht hleg
  have hlpos : 0 < l := lt_trans (div_pos two_pos (by linarith)) hl.1
  have hstab := posStable_ray N beta t l hs hlpos hl.2.le
  refine ⟨l * t, hlt, ?_, ?_, ?_⟩
  · nlinarith [hl.2]
  · rw [← hcurve]; nlinarith [hl.2]
  · rw [← hcurve]; exact hstab

/-- **Lemma R, forward direction, in CTLN coordinates.**  Stability of the full support at some
legal `(eps, delta)` forces stability at a point `(eps', delta')` on the legality boundary
`q = delta (1 - eps)/eps = 1`, with `eps' > eps`. -/
theorem lemmaR_forward_ctln (A : Matrix n n ℝ) (hA : IsAdj A) {eps delta : ℝ}
    (he0 : 0 < eps) (he1 : eps < 1) (hq : delta * (1 - eps) / eps > 1)
    (hs : HurwitzStable (-1 + W A eps delta)) :
    ∃ eps' delta' : ℝ, eps < eps' ∧ eps' < 1 ∧ 0 < delta' ∧
      delta' * (1 - eps') / eps' = 1 ∧ HurwitzStable (-1 + W A eps' delta') := by
  set t := (1 - eps) / eps with htdef
  set beta := 1 + delta / eps with hbdef
  have ht : 0 < t := div_pos (by linarith) he0
  have hleg : t * (beta - 2) > 1 := (legal_iff he0 he1).mp hq
  have hK : PosStable (K (Nmat A) beta t) := (hurwitz_iff_posStable A hA delta he0).mp hs
  obtain ⟨t', ht'0, ht't, _, hK'⟩ := lemmaR_forward (Nmat A) ht hleg hK
  refine ⟨1 / (1 + t'), 1 / t', ?_, ?_, ?_, ?_, ?_⟩
  · -- eps = 1/(1 + t) < 1/(1 + t')
    have heps : eps = 1 / (1 + t) := by
      rw [htdef]; field_simp; ring
    rw [heps]
    exact one_div_lt_one_div_of_lt (by linarith) (by linarith)
  · rw [div_lt_one (by linarith)]; linarith
  · exact one_div_pos.mpr ht'0
  · field_simp
    ring
  · have he' : 0 < 1 / (1 + t') := one_div_pos.mpr (by linarith)
    rw [hurwitz_iff_posStable A hA (1 / t') he']
    have h1 : 1 + 1 / t' / (1 / (1 + t')) = 2 + 1 / t' := by field_simp; ring
    have h2 : (1 - 1 / (1 + t')) / (1 / (1 + t')) = t' := by field_simp; ring
    rw [h1, h2]
    exact hK'

/-- Elementary half of the converse of Lemma R: shifting a curve point to the right in beta
lands in the legal region. -/
theorem curve_shift_legal {t' eta : ℝ} (ht' : 0 < t') (heta : 0 < eta) :
    t' * ((2 + 1 / t' + eta) - 2) > 1 := by
  have : t' * ((2 + 1 / t' + eta) - 2) = 1 + t' * eta := by field_simp; ring
  rw [this]
  nlinarith [mul_pos ht' heta]

end CtlnCert
