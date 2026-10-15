/-
  Stodola's necessary condition, and its CTLN corollary (the (b1) "coefficients never all
  positive" prune of n7/curve_test.py, `coeff_never`).

  * `PosCoeffs f`                   every coefficient c_0, ..., c_(natDegree f) of f is > 0
  * `PosCoeffs.mul`                 positive-coefficient polynomials are closed under products
  * `posCoeffs_linear`, `posCoeffs_quadratic`   X + a (a > 0) and X^2 + b X + c (b, c > 0)
  * `conj_root`                     a real polynomial's complex roots come in conjugate pairs
  * `stodola`                       **monic real p, all complex roots with Re < 0  =>  every
                                    coefficient c_k (k <= deg p) is > 0**
  * `coeff_pos_of_hurwitz`          M Hurwitz stable  =>  every coefficient of charpoly M is > 0
  * `coeff_pos_of_posStable`        K positive stable =>  every coefficient of charpoly(-K) =
                                    det(mu I + K) is > 0 (k = 0 is det K > 0, Lemma A.2)
  * `ctln_coeff_pos`                CTLN form: -I + W Hurwitz stable, eps > 0  =>  all coefficients of
                                    charpoly(-I + W) and of det(mu I + K(N, 1 + delta/eps, (1-eps)/eps))
                                    are > 0

  Proof route: strong induction on the degree.  A monic p of positive degree has a complex root z
  (Re z < 0).  If z is real, p = (X + a) q with a = -z > 0.  Otherwise z and conj z are distinct roots,
  (X - z)(X - conj z) = X^2 + b X + c with b = -2 Re z > 0, c = |z|^2 > 0 divides p over C, hence over R
  (`Polynomial.map_dvd_map'`).  In both cases q is monic of lower degree, its roots are roots of p, and
  the product of two positive-coefficient polynomials has positive coefficients.
-/
import CtlnCert.LemmaA2

open Matrix Polynomial ComplexConjugate

namespace CtlnCert

/-- All coefficients up to (and including) the degree are strictly positive. -/
def PosCoeffs (f : ℝ[X]) : Prop := ∀ k ≤ f.natDegree, 0 < f.coeff k

theorem PosCoeffs.nonneg {f : ℝ[X]} (hf : PosCoeffs f) (k : ℕ) : 0 ≤ f.coeff k := by
  by_cases hk : k ≤ f.natDegree
  · exact (hf k hk).le
  · rw [coeff_eq_zero_of_natDegree_lt (by omega)]

theorem PosCoeffs.ne_zero {f : ℝ[X]} (hf : PosCoeffs f) : f ≠ 0 := by
  intro h
  have := hf 0 (Nat.zero_le _)
  simp [h] at this

/-- Products of positive-coefficient polynomials have positive coefficients. -/
theorem PosCoeffs.mul {f g : ℝ[X]} (hf : PosCoeffs f) (hg : PosCoeffs g) : PosCoeffs (f * g) := by
  intro k hk
  rw [natDegree_mul hf.ne_zero hg.ne_zero] at hk
  rw [coeff_mul]
  apply Finset.sum_pos' (fun x _ => mul_nonneg (hf.nonneg _) (hg.nonneg _))
  refine ⟨(min k f.natDegree, k - min k f.natDegree), ?_, ?_⟩
  · rw [Finset.HasAntidiagonal.mem_antidiagonal]; omega
  · exact mul_pos (hf _ (min_le_right _ _)) (hg _ (by omega))

theorem posCoeffs_one : PosCoeffs (1 : ℝ[X]) := by
  intro k hk
  rw [natDegree_one] at hk
  have : k = 0 := by omega
  subst this; simp

theorem posCoeffs_linear {a : ℝ} (ha : 0 < a) : PosCoeffs (X + C a) := by
  intro k hk
  rw [natDegree_X_add_C] at hk
  interval_cases k <;> simp [ha]

/-- The quadratic `X^2 + b X + c` written in Mathlib's `C 1 * X^2 + C b * X + C c` normal form. -/
theorem posCoeffs_quadratic {b c : ℝ} (hb : 0 < b) (hc : 0 < c) :
    PosCoeffs (C 1 * X ^ 2 + C b * X + C c) := by
  intro k hk
  rw [natDegree_quadratic one_ne_zero] at hk
  interval_cases k <;> simp [coeff_X, coeff_X_pow, hb, hc]

/-- For a real polynomial, `p(conj w) = conj (p(w))` after mapping to `ℂ`. -/
theorem eval_map_conj (p : ℝ[X]) (w : ℂ) :
    (p.map Complex.ofRealHom).eval (conj w) = conj ((p.map Complex.ofRealHom).eval w) := by
  simp [eval_map, eval₂_eq_sum_range, map_sum, map_mul, map_pow]

theorem conj_root (p : ℝ[X]) {z : ℂ} (hz : (p.map Complex.ofRealHom).IsRoot z) :
    (p.map Complex.ofRealHom).IsRoot (conj z) := by
  rw [IsRoot, eval_map_conj, hz.eq_zero, map_zero]

/-- Root-hereditary step: roots of a factor are roots of the product. -/
theorem roots_of_factor {F q p : ℝ[X]} (hpq : p = F * q)
    (h : ∀ z : ℂ, (p.map Complex.ofRealHom).IsRoot z → z.re < 0) :
    ∀ z : ℂ, (q.map Complex.ofRealHom).IsRoot z → z.re < 0 := by
  intro z hz
  apply h z
  rw [hpq, Polynomial.map_mul]
  exact root_mul_left_of_isRoot _ hz

/-- Inductive form of Stodola's condition. -/
theorem stodola_aux (d : ℕ) : ∀ p : ℝ[X], p.natDegree = d → p.Monic →
    (∀ z : ℂ, (p.map Complex.ofRealHom).IsRoot z → z.re < 0) → PosCoeffs p := by
  induction d using Nat.strong_induction_on with
  | _ d ih =>
  intro p hd hp hroots
  -- splitting off a monic positive-coefficient factor of positive degree
  have key : ∀ F q : ℝ[X], F.Monic → 1 ≤ F.natDegree → PosCoeffs F → p = F * q →
      PosCoeffs p := by
    intro F q hF hFd hFp hpq
    have hq : q.Monic := hF.of_mul_monic_left (hpq ▸ hp)
    have hdeg : p.natDegree = F.natDegree + q.natDegree := by
      rw [hpq, hF.natDegree_mul hq]
    have := ih q.natDegree (by omega) q rfl hq (roots_of_factor hpq hroots)
    rw [hpq]; exact hFp.mul this
  rcases Nat.eq_zero_or_pos d with h0 | hpos
  · rw [hp.natDegree_eq_zero.mp (hd.trans h0)]
    exact posCoeffs_one
  · have hdegC : 0 < (p.map Complex.ofRealHom).degree := by
      rw [degree_map, degree_eq_natDegree hp.ne_zero]
      exact_mod_cast hd ▸ hpos
    obtain ⟨z, hz⟩ := Complex.exists_root hdegC
    have hzre := hroots z hz
    by_cases him : z.im = 0
    · -- real root r = Re z < 0: p = (X + a) q with a = -r > 0
      set r := z.re
      have hzr : z = (r : ℂ) := Complex.ext (by simp [r]) (by simp [him])
      have hpr : p.IsRoot r := by
        rw [IsRoot]
        have := hz.eq_zero
        rw [hzr, ← Complex.ofRealHom_eq_coe, eval_map_apply] at this
        simpa using this
      have hfac := (mul_divByMonic_eq_iff_isRoot (p := p) (a := r)).mpr hpr
      have hXa : X - C r = X + C (-r) := by rw [C_neg, sub_eq_add_neg]
      refine key (X + C (-r)) (p /ₘ (X - C r)) (monic_X_add_C _) (by rw [natDegree_X_add_C])
        (posCoeffs_linear (by linarith)) ?_
      rw [← hXa, hfac]
    · -- non-real root: the conjugate pair gives a real quadratic factor
      set b := -(2 * z.re)
      set c := Complex.normSq z
      have hb : 0 < b := by simp only [b]; linarith
      have hc : 0 < c := Complex.normSq_pos.mpr (fun h => him (by simp [h]))
      set Q : ℝ[X] := C 1 * X ^ 2 + C b * X + C c with hQ
      have hQmon : Q.Monic := by
        rw [hQ, C_1, one_mul, add_assoc]
        exact monic_X_pow_add (lt_of_le_of_lt degree_linear_le (by norm_num))
      have hQmap : Q.map Complex.ofRealHom = (X - C z) * (X - C (conj z)) := by
        have h1 : (X - C z) * (X - C (conj z)) = X ^ 2 - C (z + conj z) * X + C (z * conj z) := by
          rw [C_add, C_mul]; ring
        rw [h1, Complex.add_conj, Complex.mul_conj, hQ]
        simp only [Polynomial.map_add, Polynomial.map_mul, map_C, map_X, Polynomial.map_pow,
          Complex.ofRealHom_eq_coe]
        simp only [b, c, Complex.ofReal_neg, Complex.ofReal_one, C_neg, C_1, one_mul]
        ring
      have hdvdC : Q.map Complex.ofRealHom ∣ p.map Complex.ofRealHom := by
        rw [hQmap]
        set P := p.map Complex.ofRealHom
        have h1 := (mul_divByMonic_eq_iff_isRoot (p := P) (a := z)).mpr hz
        set R1 := P /ₘ (X - C z)
        have hcz : P.IsRoot (conj z) := conj_root p hz
        have hne : conj z - z ≠ 0 := by
          intro h
          have : (conj z - z).im = 0 := by rw [h]; rfl
          simp at this
          exact him (by linarith)
        have hR1 : R1.IsRoot (conj z) := by
          rw [IsRoot] at hcz ⊢
          rw [← h1, eval_mul, eval_sub, eval_X, eval_C] at hcz
          exact (mul_eq_zero.mp hcz).resolve_left hne
        have h2 := (mul_divByMonic_eq_iff_isRoot (p := R1) (a := conj z)).mpr hR1
        refine ⟨R1 /ₘ (X - C (conj z)), ?_⟩
        rw [mul_assoc, h2, h1]
      have hdvd : Q ∣ p := (map_dvd_map' Complex.ofRealHom).mp hdvdC
      obtain ⟨q, hq⟩ := hdvd
      refine key Q q hQmon ?_ (posCoeffs_quadratic hb hc) hq
      rw [hQ, natDegree_quadratic one_ne_zero]; norm_num

/-- **Stodola's necessary condition.**  A monic real polynomial all of whose complex roots have
negative real part has every coefficient `c_0, ..., c_(deg p)` strictly positive. -/
theorem stodola (p : ℝ[X]) (hp : p.Monic)
    (h : ∀ z : ℂ, (p.map Complex.ofRealHom).IsRoot z → z.re < 0) :
    ∀ k ≤ p.natDegree, 0 < p.coeff k :=
  stodola_aux p.natDegree p rfl hp h

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- A Hurwitz stable real matrix has a characteristic polynomial with all coefficients positive. -/
theorem coeff_pos_of_hurwitz (M : Matrix n n ℝ) (h : HurwitzStable M) :
    ∀ k ≤ Fintype.card n, 0 < M.charpoly.coeff k := by
  have hc : (cplx M).charpoly = M.charpoly.map Complex.ofRealHom := by
    rw [← charpoly_map]; rfl
  intro k hk
  apply stodola M.charpoly (charpoly_monic M) _ k (by rwa [charpoly_natDegree_eq_dim])
  intro z hz
  rw [← hc, ← eig_iff_isRoot_charpoly] at hz
  exact h z hz

/-- Positive stable `K`: every coefficient of `charpoly(-K)`, i.e. of `mu ↦ det(mu I + K)`, is
positive.  The case `k = 0` is `det K > 0` (Lemma A.2). -/
theorem coeff_pos_of_posStable (K : Matrix n n ℝ) (h : PosStable K) :
    ∀ k ≤ Fintype.card n, 0 < (-K).charpoly.coeff k := by
  intro k hk
  exact stodola (-K).charpoly (charpoly_monic _) (negK_root_re_neg K h) k
    (by rwa [charpoly_natDegree_eq_dim])

/-- The k = 0 case recovers Lemma A.2 (consistency check): `(charpoly (-K)).coeff 0 = det K`. -/
theorem coeff_zero_negK (K : Matrix n n ℝ) : (-K).charpoly.coeff 0 = K.det := by
  rw [coeff_zero_eq_eval_zero, eval_charpoly]; simp

/-- **CTLN corollary.**  If `-I + W` is Hurwitz stable (`eps > 0`), then every coefficient of
`charpoly(-I + W)` and every coefficient of `det(mu I + K(N, 1 + delta/eps, (1 - eps)/eps))` is
positive.  Contrapositive: if for some `k` a coefficient is `<= 0`, the full support is not stable
(the (b1) test `coeff_never` of n7/curve_test.py). -/
theorem ctln_coeff_pos (A : Matrix n n ℝ) (hA : IsAdj A) {eps : ℝ} (delta : ℝ) (he : 0 < eps)
    (hs : HurwitzStable (-1 + W A eps delta)) :
    (∀ k ≤ Fintype.card n, 0 < (-1 + W A eps delta).charpoly.coeff k) ∧
    (∀ k ≤ Fintype.card n,
      0 < (-(K (Nmat A) (1 + delta / eps) ((1 - eps) / eps))).charpoly.coeff k) :=
  ⟨coeff_pos_of_hurwitz _ hs,
   coeff_pos_of_posStable _ ((hurwitz_iff_posStable A hA delta he).mp hs)⟩

end CtlnCert
