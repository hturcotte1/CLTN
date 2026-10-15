/-
  Bridge to the project conventions (CLAUDE.md section 2, notes/lemmas.md section 0).

  A[i][j] = 1 iff there is an edge i -> j (0/1 entries, zero diagonal).
  W[i][j] = 0 if i = j;  -1 + eps if A[j][i] = 1;  -1 - delta otherwise.   (note the transpose)

  * `closed_form`      I - W = -delta I + (1 + delta) J - (eps + delta) A^T
                       (CLAUDE.md writes A_sigma; with W reading A[j][i] it is A^T, as already logged
                        in ANOMALIES.md 2026-10-14; this file proves the transposed form)
  * `IminusW_eq`       I - W = eps * K(N, beta, t) with N = J - I - A^T,
                       beta = 1 + delta/eps, t = (1 - eps)/eps
  * `legal_iff`        for 0 < eps < 1, delta > 0:  q > 1  <->  t (beta - 2) > 1,
                       q = delta (1 - eps)/eps
  * `hurwitz_iff_posStable`  -I + W Hurwitz stable  <->  K(N, beta, t) positive stable  (eps > 0)
-/
import CtlnCert.RayIdentity

open Matrix

namespace CtlnCert

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- The CTLN weight matrix, exactly as in CLAUDE.md section 2 (reads `A j i`). -/
noncomputable def W (A : Matrix n n ℝ) (eps delta : ℝ) : Matrix n n ℝ :=
  Matrix.of fun i j => if i = j then 0 else if A j i = 1 then -1 + eps else -1 - delta

/-- The receiver-first non-edge matrix `N = J - I - A^T`. -/
def Nmat (A : Matrix n n ℝ) : Matrix n n ℝ := J ℝ - 1 - Aᵀ

/-- Hypotheses on an adjacency matrix: 0/1 entries and zero diagonal. -/
def IsAdj (A : Matrix n n ℝ) : Prop := (∀ i j, A i j = 0 ∨ A i j = 1) ∧ ∀ i, A i i = 0

omit [Fintype n] in
/-- CLAUDE.md section 2 closed form (with the transpose made explicit). -/
theorem closed_form (A : Matrix n n ℝ) (hA : IsAdj A) (eps delta : ℝ) :
    1 - W A eps delta = (-delta) • (1 : Matrix n n ℝ) + (1 + delta) • J ℝ - (eps + delta) • Aᵀ := by
  ext i j
  by_cases hij : i = j
  · subst hij
    simp [W, J, hA.2 i]
  · rcases hA.1 j i with h | h
    · simp [W, J, hij, h]
      ring
    · simp [W, J, hij, h]
      ring

omit [Fintype n] in
/-- The (beta, t) reparametrisation: `I - W = eps * K(N, 1 + delta/eps, (1 - eps)/eps)`. -/
theorem IminusW_eq (A : Matrix n n ℝ) (hA : IsAdj A) (eps delta : ℝ) (he : eps ≠ 0) :
    1 - W A eps delta = eps • K (Nmat A) (1 + delta / eps) ((1 - eps) / eps) := by
  rw [closed_form A hA]
  ext i j
  by_cases hij : i = j
  · subst hij
    simp [K, Nmat, J, hA.2 i]
    field_simp
    ring
  · simp [K, Nmat, J, hij]
    field_simp
    ring

/-- Legality in (beta, t) coordinates: `q > 1 <-> t (beta - 2) > 1`. -/
theorem legal_iff {eps delta : ℝ} (he0 : 0 < eps) (he1 : eps < 1) :
    delta * (1 - eps) / eps > 1 ↔ ((1 - eps) / eps) * ((1 + delta / eps) - 2) > 1 := by
  have h1 : 0 < 1 - eps := by linarith
  have key : ((1 - eps) / eps) * ((1 + delta / eps) - 2) - 1
      = (delta * (1 - eps) / eps - 1) / eps := by
    field_simp; ring
  constructor
  · intro h
    have : 0 < (delta * (1 - eps) / eps - 1) / eps := div_pos (by linarith) he0
    linarith
  · intro h
    have : 0 < (delta * (1 - eps) / eps - 1) / eps := by linarith
    have := (div_pos_iff.mp this)
    rcases this with ⟨h2, _⟩ | ⟨_, h3⟩
    · linarith
    · linarith

/-- Stability condition 3 of CLAUDE.md section 2 in (beta, t) coordinates:
`-I + W` Hurwitz stable iff `K(N, beta, t)` positive stable, for `eps > 0`. -/
theorem hurwitz_iff_posStable (A : Matrix n n ℝ) (hA : IsAdj A) {eps : ℝ} (delta : ℝ)
    (he : 0 < eps) :
    HurwitzStable (-1 + W A eps delta)
      ↔ PosStable (K (Nmat A) (1 + delta / eps) ((1 - eps) / eps)) := by
  set Km := K (Nmat A) (1 + delta / eps) ((1 - eps) / eps)
  have hW : -1 + W A eps delta = (0 : ℝ) • (1 : Matrix n n ℝ) + (-eps) • Km := by
    have := IminusW_eq A hA eps delta he.ne'
    rw [zero_smul, zero_add, neg_smul, ← this]
    abel
  have hne : ((-eps : ℝ) : ℂ) ≠ 0 := by
    rw [Complex.ofReal_ne_zero]; linarith
  have heq : ∀ mu : ℂ, Eig (cplx (-1 + W A eps delta)) mu ↔ Eig (cplx Km) (-mu / (eps : ℂ)) := by
    intro mu
    rw [hW, cplx_affine, eig_affine _ _ _ _ hne]
    have : (mu - ((0 : ℝ) : ℂ)) / ((-eps : ℝ) : ℂ) = -mu / (eps : ℂ) := by
      push_cast; rw [sub_zero, div_neg, neg_div]
    rw [this]
  have hepsC : (eps : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr he.ne'
  constructor
  · intro h nu hnu
    have hmu : Eig (cplx (-1 + W A eps delta)) (-(eps : ℂ) * nu) := by
      rw [heq]
      have : -(-(eps : ℂ) * nu) / (eps : ℂ) = nu := by field_simp
      rw [this]; exact hnu
    have := h _ hmu
    have hre : (-(eps : ℂ) * nu).re = -eps * nu.re := by simp [Complex.mul_re]
    rw [hre] at this
    nlinarith
  · intro h mu hmu
    have := h _ ((heq mu).mp hmu)
    have hre : (-mu / (eps : ℂ)).re = -mu.re / eps := by
      simp [Complex.div_ofReal_re]
    rw [hre] at this
    have : 0 < -mu.re := by
      by_contra hc; push Not at hc
      have : -mu.re / eps ≤ 0 := div_nonpos_of_nonpos_of_nonneg hc he.le
      linarith
    linarith

end CtlnCert
