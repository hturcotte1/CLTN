/-
  Basic definitions for the CTLN ray-reduction formalisation.

  Conventions follow notes/lemmas.md section 0 of the project:
    K(beta, t) := I + beta * N + t * J,   J = all-ones matrix.
  For a real matrix M, "positive stable" means every complex eigenvalue of M has
  strictly positive real part.  Eigenvalues are defined by eigenvectors over the
  complex numbers, and `eig_iff_isRoot_charpoly` shows this is the same as being a
  root of the characteristic polynomial (the object the project's Hurwitz test uses).
-/
import Mathlib

open Matrix

namespace CtlnCert

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- The all-ones matrix `J`. -/
def J (R : Type*) [One R] : Matrix n n R := Matrix.of fun _ _ => 1

/-- `K(beta, t) = I + beta N + t J`. -/
def K {R : Type*} [CommRing R] (N : Matrix n n R) (beta t : R) : Matrix n n R :=
  1 + beta • N + t • J R

/-- `mu` is an eigenvalue of the complex matrix `M`: there is a nonzero `v` with `M v = mu v`. -/
def Eig (M : Matrix n n ℂ) (mu : ℂ) : Prop :=
  ∃ v : n → ℂ, v ≠ 0 ∧ M *ᵥ v = mu • v

/-- Complexification of a real matrix. -/
noncomputable def cplx (M : Matrix n n ℝ) : Matrix n n ℂ := M.map (fun x : ℝ => (x : ℂ))

/-- Positive stability of a real matrix: every complex eigenvalue has positive real part. -/
def PosStable (M : Matrix n n ℝ) : Prop :=
  ∀ mu : ℂ, Eig (cplx M) mu → 0 < mu.re

/-- Hurwitz stability of a real matrix: every complex eigenvalue has negative real part. -/
def HurwitzStable (M : Matrix n n ℝ) : Prop :=
  ∀ mu : ℂ, Eig (cplx M) mu → mu.re < 0

/-- Eigenvalues (eigenvector definition) are exactly the roots of the characteristic polynomial. -/
theorem eig_iff_isRoot_charpoly (M : Matrix n n ℂ) (mu : ℂ) :
    Eig M mu ↔ M.charpoly.IsRoot mu := by
  rw [Polynomial.IsRoot, eval_charpoly, ← exists_mulVec_eq_zero_iff]
  constructor
  · rintro ⟨v, hv, h⟩
    refine ⟨v, hv, ?_⟩
    rw [sub_mulVec, h, scalar_apply, sub_eq_zero]
    ext i
    simp [mulVec_diagonal]
  · rintro ⟨v, hv, h⟩
    refine ⟨v, hv, ?_⟩
    rw [sub_mulVec, sub_eq_zero] at h
    rw [← h, scalar_apply]
    ext i
    simp [mulVec_diagonal]

/-- Positive stability, stated through the characteristic polynomial of the real matrix. -/
theorem posStable_iff_charpoly (M : Matrix n n ℝ) :
    PosStable M ↔ ∀ z : ℂ, (M.charpoly.map Complex.ofRealHom).IsRoot z → 0 < z.re := by
  have hc : (cplx M).charpoly = M.charpoly.map Complex.ofRealHom := by
    rw [← charpoly_map]; rfl
  unfold PosStable
  simp only [eig_iff_isRoot_charpoly, hc]

end CtlnCert
