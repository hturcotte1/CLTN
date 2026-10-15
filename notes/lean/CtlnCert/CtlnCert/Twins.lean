/-
  Remark J.4' of notes/lemma_J_equivalence.md (exact twin eigenvalue), full-support form.

  u ≠ v are *in-twins* (relative to the full support) if `A i u = A i v` for every `i ∉ {u, v}`
  (same in-neighbours among the other nodes).  With W from Bridge.lean (W i j reads A j i):

  * `twin_left_eig`          generic: rows u, v of M agree off columns {u, v} and
                             M u u - M v u = lam = M v v - M u v  =>  xiᵀ M = lam xiᵀ, xi = e_u - e_v
  * `twin_W_nonadj`          non-adjacent in-twins (A u v = A v u = 0):  xiᵀ W = (1 + delta) xiᵀ
  * `twin_W_bidir`           bidirected in-twins   (A u v = A v u = 1):  xiᵀ W = (1 - eps) xiᵀ
  * `twin_negIW_nonadj`      xiᵀ (-1 + W) = delta xiᵀ
  * `twin_negIW_bidir`       xiᵀ (-1 + W) = (-eps) xiᵀ
  * `charpoly_root_of_left_eig`   a left eigenvector gives a root of charpoly (det(lam I - M) = 0;
                             equivalently, via Matrix.charpoly_transpose, an eigenvalue of Mᵀ)
  * `eig_of_left_eig`        ... hence an eigenvalue (Eig, right eigenvector over C) of cplx M
  * `twin_charpoly_root_nonadj`, `twin_charpoly_root_bidir`   delta resp. -eps is a root of
                             charpoly(-1 + W)
  * `twin_not_hurwitz`       **non-adjacent in-twins => ¬ HurwitzStable (-1 + W A eps delta)** for every
                             eps and every delta > 0
  * `twin_not_hurwitz_induced`  the same inside an induced subgraph G[sigma] (sigma = range of an
                             injective f): twins relative to sigma make -1 + W_sigma not Hurwitz stable
  Not used: IsAdj, eps > 0, eps < 1, legality, theta.  Only `delta > 0` is used, for the instability.
-/
import CtlnCert.Domination

open Matrix Polynomial

namespace CtlnCert

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- The twin vector `xi = e_u - e_v`. -/
noncomputable def twinVec (u v : n) : n → ℝ := Pi.single u 1 - Pi.single v 1

omit [Fintype n] in
theorem twinVec_ne_zero {u v : n} (huv : u ≠ v) : twinVec u v ≠ 0 := by
  intro h
  have := congrFun h u
  simp [twinVec, huv] at this

/-- Generic twin-row computation: `xiᵀ M = lam xiᵀ` for `xi = e_u - e_v`. -/
theorem twin_left_eig (M : Matrix n n ℝ) {u v : n} (huv : u ≠ v) (lam : ℝ)
    (hrow : ∀ j, j ≠ u → j ≠ v → M u j = M v j)
    (hu : M u u - M v u = lam) (hv : M u v - M v v = -lam) :
    twinVec u v ᵥ* M = lam • twinVec u v := by
  ext j
  rw [twinVec, sub_vecMul, single_one_vecMul, single_one_vecMul]
  simp only [Pi.sub_apply, row_apply, Pi.smul_apply, smul_eq_mul, Pi.single_apply]
  by_cases hju : j = u
  · subst hju; simp [huv, hu]
  · by_cases hjv : j = v
    · subst hjv; simp [hju, hv]
    · simp [hju, hjv, hrow j hju hjv]

/-- In-twin hypothesis (relative to the full support). -/
def InTwins (A : Matrix n n ℝ) (u v : n) : Prop :=
  u ≠ v ∧ ∀ i, i ≠ u → i ≠ v → A i u = A i v

omit [Fintype n] in
theorem W_twin_rows (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v) (eps delta : ℝ) :
    ∀ j, j ≠ u → j ≠ v → W A eps delta u j = W A eps delta v j := by
  intro j hju hjv
  simp [W, Ne.symm hju, Ne.symm hjv, h.2 j hju hjv]

/-- **Non-adjacent in-twins:** `xiᵀ W = (1 + delta) xiᵀ`. -/
theorem twin_W_nonadj (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v)
    (huv : A u v = 0) (hvu : A v u = 0) (eps delta : ℝ) :
    twinVec u v ᵥ* W A eps delta = (1 + delta) • twinVec u v := by
  apply twin_left_eig _ h.1 _ (W_twin_rows A h eps delta)
  · simp [W, Ne.symm h.1, huv]; ring
  · simp [W, h.1, hvu]; ring

/-- **Bidirected in-twins:** `xiᵀ W = (1 - eps) xiᵀ`. -/
theorem twin_W_bidir (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v)
    (huv : A u v = 1) (hvu : A v u = 1) (eps delta : ℝ) :
    twinVec u v ᵥ* W A eps delta = (1 - eps) • twinVec u v := by
  apply twin_left_eig _ h.1 _ (W_twin_rows A h eps delta)
  · simp [W, Ne.symm h.1, huv]; ring
  · simp [W, h.1, hvu]; ring

theorem vecMul_negOne_add (M : Matrix n n ℝ) (x : n → ℝ) :
    x ᵥ* (-1 + M) = -x + x ᵥ* M := by
  rw [vecMul_add, vecMul_neg, vecMul_one]

/-- Non-adjacent in-twins: `xiᵀ (-1 + W) = delta xiᵀ`. -/
theorem twin_negIW_nonadj (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v)
    (huv : A u v = 0) (hvu : A v u = 0) (eps delta : ℝ) :
    twinVec u v ᵥ* (-1 + W A eps delta) = delta • twinVec u v := by
  rw [vecMul_negOne_add, twin_W_nonadj A h huv hvu, add_smul, one_smul]
  abel

/-- Bidirected in-twins: `xiᵀ (-1 + W) = (-eps) xiᵀ`. -/
theorem twin_negIW_bidir (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v)
    (huv : A u v = 1) (hvu : A v u = 1) (eps delta : ℝ) :
    twinVec u v ᵥ* (-1 + W A eps delta) = (-eps) • twinVec u v := by
  rw [vecMul_negOne_add, twin_W_bidir A h huv hvu, sub_smul, one_smul, neg_smul]
  abel

/-- A (nonzero) left eigenvector gives a root of the characteristic polynomial. -/
theorem charpoly_root_of_left_eig (M : Matrix n n ℝ) (x : n → ℝ) (hx : x ≠ 0) (lam : ℝ)
    (h : x ᵥ* M = lam • x) : M.charpoly.IsRoot lam := by
  rw [IsRoot, eval_charpoly, ← exists_vecMul_eq_zero_iff]
  refine ⟨x, hx, ?_⟩
  rw [vecMul_sub, h, scalar_apply]
  ext i
  simp

/-- The same, as the transpose statement: `lam` is a root of `charpoly Mᵀ` (Matrix.charpoly_transpose). -/
theorem charpoly_transpose_root_of_left_eig (M : Matrix n n ℝ) (x : n → ℝ) (hx : x ≠ 0) (lam : ℝ)
    (h : x ᵥ* M = lam • x) : Mᵀ.charpoly.IsRoot lam := by
  rw [charpoly_transpose]; exact charpoly_root_of_left_eig M x hx lam h

/-- A real left eigenvector gives a complex eigenvalue (right eigenvector, `Eig`) of `cplx M`. -/
theorem eig_of_left_eig (M : Matrix n n ℝ) (x : n → ℝ) (hx : x ≠ 0) (lam : ℝ)
    (h : x ᵥ* M = lam • x) : Eig (cplx M) (lam : ℂ) := by
  have hc : (cplx M).charpoly = M.charpoly.map Complex.ofRealHom := by
    rw [← charpoly_map]; rfl
  rw [eig_iff_isRoot_charpoly, hc, IsRoot, ← Complex.ofRealHom_eq_coe, eval_map_apply,
    (charpoly_root_of_left_eig M x hx lam h).eq_zero, map_zero]

/-- Non-adjacent in-twins: `delta` is a root of `charpoly(-1 + W)`. -/
theorem twin_charpoly_root_nonadj (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v)
    (huv : A u v = 0) (hvu : A v u = 0) (eps delta : ℝ) :
    (-1 + W A eps delta).charpoly.IsRoot delta :=
  charpoly_root_of_left_eig _ _ (twinVec_ne_zero h.1) _ (twin_negIW_nonadj A h huv hvu eps delta)

/-- Bidirected in-twins: `-eps` is a root of `charpoly(-1 + W)`. -/
theorem twin_charpoly_root_bidir (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v)
    (huv : A u v = 1) (hvu : A v u = 1) (eps delta : ℝ) :
    (-1 + W A eps delta).charpoly.IsRoot (-eps) :=
  charpoly_root_of_left_eig _ _ (twinVec_ne_zero h.1) _ (twin_negIW_bidir A h huv hvu eps delta)

/-- **Remark J.4' (CTLN corollary).**  Two non-adjacent in-twins (relative to the full support) make
`-1 + W` not Hurwitz stable, for every `eps` and every `delta > 0`: `delta` is an eigenvalue. -/
theorem twin_not_hurwitz (A : Matrix n n ℝ) {u v : n} (h : InTwins A u v)
    (huv : A u v = 0) (hvu : A v u = 0) (eps : ℝ) {delta : ℝ} (hd : 0 < delta) :
    ¬ HurwitzStable (-1 + W A eps delta) := by
  intro hs
  have := hs _ (eig_of_left_eig _ _ (twinVec_ne_zero h.1) _ (twin_negIW_nonadj A h huv hvu eps delta))
  simp only [Complex.ofReal_re] at this
  linarith

omit [Fintype n] in
/-- Induced form: if `f : m → n` is injective and `u ≠ v` in `m` are non-adjacent in-twins of the induced
graph `A[f, f]` (twins relative to `sigma = range f`), then `-1 + W_sigma` is not Hurwitz stable. -/
theorem twin_not_hurwitz_induced {m : Type*} [Fintype m] [DecidableEq m] (A : Matrix n n ℝ)
    (f : m → n) (hf : Function.Injective f) {u v : m} (h : InTwins (A.submatrix f f) u v)
    (huv : A (f u) (f v) = 0) (hvu : A (f v) (f u) = 0) (eps : ℝ) {delta : ℝ} (hd : 0 < delta) :
    ¬ HurwitzStable (-1 + (W A eps delta).submatrix f f) := by
  rw [← W_submatrix A f hf]
  exact twin_not_hurwitz _ h huv hvu eps hd

end CtlnCert
