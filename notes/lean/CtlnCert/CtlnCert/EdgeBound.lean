/-
  The e_2 edge bound (notes/lemma_M_refuted.md section 5; ANOMALIES.md 2026-10-14 audit entry).

  Setting: A is a 0/1 adjacency matrix with zero diagonal, N = J - I - Aᵀ (N i j = 1 iff i ≠ j and
  j does not send to i), K = I + beta N + t J.  Unordered pairs {i, j} are indexed by i < j.
    B  = #bidirected pairs (A i j = A j i = 1)
    Nn = #non-adjacent pairs (A i j = A j i = 0)
    O  = #one-way pairs.

  * `e2`                       sum over i < j of (M i i M j j - M i j M j i)
  * `e2_eq_sum_minors`         e2 M = Mathlib's sum of the 2x2 principal minors (powersetCard 2)
  * `e2_closed_form`           over any nontrivial commutative ring:
                               e2 K = B(1+2t) + Nn(1-beta)(1+beta+2t) + O(1 - t(beta-2))
  * `pair_count`               B + Nn + O = #pairs
  * `two_mul_card_pairs`       2 #pairs = m(m-1)
  * `arcs_eq`                  #arcs = 2B + O
  * `bidir_gt_nonadj`          beta > 2, t > 0, t(beta-2) > 1, e2 K > 0  =>  Nn < B
  * `arcs_lower_bound`         same hypotheses  =>  2 #arcs >= m(m-1) + 2
  * `e2_pos_of_posStable`      K positive stable, m >= 2  =>  e2 K > 0   (via LemmaA2.sum_minors_two_pos)
  * `edge_bound_posStable`     K(N, beta, t) positive stable at a legal (beta, t), m >= 2
                               =>  Nn < B  and  2 #arcs >= m(m-1) + 2
  * `edge_bound_ctln`          the same in (eps, delta) language: -I + W Hurwitz stable,
                               0 < eps < 1, q = delta(1-eps)/eps > 1, m >= 2
-/
import CtlnCert.LemmaA2
import CtlnCert.IVTStep

open Matrix

namespace CtlnCert

variable {n : Type*} [Fintype n] [LinearOrder n]

/-- Ordered representatives `(i, j)`, `i < j`, of the unordered pairs of distinct nodes. -/
def pairs (n : Type*) [Fintype n] [LinearOrder n] : Finset (n × n) :=
  Finset.univ.filter (fun p : n × n => p.1 < p.2)

/-- `e_2(M)`: the sum of the 2x2 principal minors, indexed by pairs `i < j`. -/
def e2 {R : Type*} [CommRing R] (M : Matrix n n R) : R :=
  ∑ p ∈ pairs n, (M p.1 p.1 * M p.2 p.2 - M p.1 p.2 * M p.2 p.1)

section counts

variable {R : Type*} [CommRing R] [DecidableEq R]

/-- Number of bidirected unordered pairs. -/
def nB (A : Matrix n n R) : ℕ :=
  ((pairs n).filter (fun p => A p.1 p.2 = 1 ∧ A p.2 p.1 = 1)).card

/-- Number of non-adjacent unordered pairs. -/
def nN (A : Matrix n n R) : ℕ :=
  ((pairs n).filter (fun p => A p.1 p.2 = 0 ∧ A p.2 p.1 = 0)).card

/-- Number of one-way unordered pairs. -/
def nO (A : Matrix n n R) : ℕ :=
  ((pairs n).filter
    (fun p => (A p.1 p.2 = 1 ∧ A p.2 p.1 = 0) ∨ (A p.1 p.2 = 0 ∧ A p.2 p.1 = 1))).card

/-- Number of arcs (ordered pairs `(i, j)` with `A i j = 1`). -/
def arcs (A : Matrix n n R) : ℕ :=
  (Finset.univ.filter (fun p : n × n => A p.1 p.2 = 1)).card

end counts

/-! ### Combinatorics of pairs -/

/-- Splitting a sum over `n × n` into pairs `i < j` (both orders) and the diagonal. -/
theorem sum_prod_split {M : Type*} [AddCommMonoid M] (f : n × n → M) :
    ∑ p, f p = ∑ p ∈ pairs n, (f p + f p.swap) + ∑ i, f (i, i) := by
  rw [Finset.sum_add_distrib]
  rw [← Finset.sum_filter_add_sum_filter_not Finset.univ (fun p : n × n => p.1 < p.2)]
  have h2 : Finset.univ.filter (fun p : n × n => ¬ p.1 < p.2)
      = Finset.univ.filter (fun p : n × n => p.2 < p.1)
        ∪ Finset.univ.filter (fun p : n × n => p.1 = p.2) := by
    ext p
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_union, not_lt]
    constructor
    · intro h; rcases lt_or_eq_of_le h with h | h
      · exact Or.inl h
      · exact Or.inr h.symm
    · rintro (h | h)
      · exact h.le
      · exact h.ge
  have hdisj : Disjoint (Finset.univ.filter (fun p : n × n => p.2 < p.1))
      (Finset.univ.filter (fun p : n × n => p.1 = p.2)) := by
    rw [Finset.disjoint_filter]
    intro p _ h h'
    rw [h'] at h; exact lt_irrefl _ h
  rw [h2, Finset.sum_union hdisj, ← add_assoc]
  congr 1
  · congr 1
    apply Finset.sum_nbij' Prod.swap Prod.swap
    · intro p hp; simpa [pairs] using hp
    · intro p hp; simpa [pairs] using hp
    · intro p _; rfl
    · intro p _; rfl
    · intro p _; rfl
  · apply Finset.sum_nbij' (fun p => p.1) (fun i => (i, i))
    · intro p _; exact Finset.mem_univ _
    · intro i _; simp
    · intro p hp
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hp
      ext
      · rfl
      · exact hp
    · intro i _; rfl
    · intro p hp
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hp
      obtain ⟨a, b⟩ := p
      simp only at hp
      subst hp; rfl

/-- `2 * #pairs = m (m - 1)`. -/
theorem two_mul_card_pairs : 2 * (pairs n).card = Fintype.card n * (Fintype.card n - 1) := by
  have h := sum_prod_split (n := n) (fun _ => (1 : ℕ))
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_prod, smul_eq_mul, mul_one] at h
  -- h : m * m = #pairs * 2 + m
  have : Fintype.card n * Fintype.card n = 2 * (pairs n).card + Fintype.card n := by
    rw [h]; ring
  rcases Nat.eq_zero_or_pos (Fintype.card n) with h0 | h0
  · rw [h0] at this ⊢; omega
  · have e : Fintype.card n * (Fintype.card n - 1) = Fintype.card n * Fintype.card n - Fintype.card n := by
      rw [Nat.mul_sub, mul_one]
    rw [e]; omega

omit [Fintype n] in
/-- The 2x2 principal minor on `{i, j}`. -/
theorem det_submatrix_pair {R : Type*} [CommRing R] (M : Matrix n n R) {i j : n} (hij : i ≠ j) :
    (M.submatrix (Subtype.val : ({i, j} : Finset n) → n) Subtype.val).det
      = M i i * M j j - M i j * M j i := by
  have hc : Fintype.card ({i, j} : Finset n) = 2 := by
    rw [Fintype.card_coe, Finset.card_pair hij]
  let e : Fin 2 ≃ ({i, j} : Finset n) := (Fintype.equivFinOfCardEq hc).symm
  rw [← det_submatrix_equiv_self e, det_fin_two]
  simp only [submatrix_apply]
  have h0 := (e 0).2
  have h1 := (e 1).2
  have hne : ((e 0 : ({i, j} : Finset n)) : n) ≠ e 1 := by
    intro h
    have := e.injective (Subtype.ext h)
    exact absurd this (by decide)
  simp only [Finset.mem_insert, Finset.mem_singleton] at h0 h1
  rcases h0 with h0 | h0 <;> rcases h1 with h1 | h1 <;> rw [h0, h1] at hne ⊢ <;>
    first | exact absurd rfl hne | ring

/-- `e2` is Mathlib's sum of the 2x2 principal minors. -/
theorem e2_eq_sum_minors {R : Type*} [CommRing R] (M : Matrix n n R) :
    e2 M = ∑ s ∈ Finset.univ.powersetCard 2,
      (M.submatrix (Subtype.val : s → n) (Subtype.val : s → n)).det := by
  unfold e2
  apply Finset.sum_bij (fun p _ => ({p.1, p.2} : Finset n))
  · intro p hp
    simp only [pairs, Finset.mem_filter, Finset.mem_univ, true_and] at hp
    rw [Finset.mem_powersetCard]
    exact ⟨Finset.subset_univ _, Finset.card_pair hp.ne⟩
  · intro p hp q hq hpq
    simp only [pairs, Finset.mem_filter, Finset.mem_univ, true_and] at hp hq
    have a1 : p.1 ∈ ({q.1, q.2} : Finset n) := by rw [← hpq]; simp
    have a2 : p.2 ∈ ({q.1, q.2} : Finset n) := by rw [← hpq]; simp
    simp only [Finset.mem_insert, Finset.mem_singleton] at a1 a2
    obtain ⟨p1, p2⟩ := p
    obtain ⟨q1, q2⟩ := q
    simp only at hp hq a1 a2
    rcases a1 with a1 | a1 <;> rcases a2 with a2 | a2 <;> subst a1 <;> subst a2
    · exact absurd hp (lt_irrefl _)
    · rfl
    · exact absurd (hp.trans hq) (lt_irrefl _)
    · exact absurd hp (lt_irrefl _)
  · intro s hs
    rw [Finset.mem_powersetCard] at hs
    obtain ⟨x, y, hxy, rfl⟩ := Finset.card_eq_two.mp hs.2
    rcases lt_or_gt_of_ne hxy with h | h
    · exact ⟨(x, y), by simp [pairs, h], rfl⟩
    · exact ⟨(y, x), by simp [pairs, h], Finset.pair_comm y x⟩
  · intro p hp
    simp only [pairs, Finset.mem_filter, Finset.mem_univ, true_and] at hp
    rw [det_submatrix_pair M hp.ne]

/-! ### The closed form -/

section closed

variable {R : Type*} [CommRing R] [DecidableEq R] [Nontrivial R]

/-- 0/1 entries and zero diagonal, over a general ring. -/
def IsAdjR (A : Matrix n n R) : Prop := (∀ i j, A i j = 0 ∨ A i j = 1) ∧ ∀ i, A i i = 0

omit [Fintype n] [DecidableEq R] [Nontrivial R] in
theorem K_apply (A : Matrix n n R) (beta t : R) (i j : n) :
    K (J R - 1 - Aᵀ) beta t i j = (if i = j then 1 else 0) + beta * ((1 : R) - (if i = j then 1 else 0) - A j i) + t := by
  simp [K, J, Matrix.one_apply, transpose_apply]

omit [Fintype n] in
/-- Per-pair minor, by type. -/
theorem minor_pair (A : Matrix n n R) (hA : IsAdjR A) (beta t : R) {i j : n} (hij : i ≠ j) :
    let Km := K (J R - 1 - Aᵀ) beta t
    Km i i * Km j j - Km i j * Km j i
      = (if A i j = 1 ∧ A j i = 1 then 1 else 0) * (1 + 2 * t)
        + (if A i j = 0 ∧ A j i = 0 then 1 else 0) * ((1 - beta) * (1 + beta + 2 * t))
        + (if (A i j = 1 ∧ A j i = 0) ∨ (A i j = 0 ∧ A j i = 1) then 1 else 0)
            * (1 - t * (beta - 2)) := by
  intro Km
  have h01 : (0 : R) ≠ 1 := zero_ne_one
  have h10 : (1 : R) ≠ 0 := one_ne_zero
  rcases hA.1 i j with h1 | h1 <;> rcases hA.1 j i with h2 | h2 <;>
    simp [Km, K_apply, hij, Ne.symm hij, hA.2 i, hA.2 j, h1, h2, h01, h10] <;> ring

/-- **e_2 closed form** (notes/lemma_M_refuted.md section 5), over any nontrivial commutative ring. -/
theorem e2_closed_form (A : Matrix n n R) (hA : IsAdjR A) (beta t : R) :
    e2 (K (J R - 1 - Aᵀ) beta t)
      = (nB A : R) * (1 + 2 * t) + (nN A : R) * ((1 - beta) * (1 + beta + 2 * t))
        + (nO A : R) * (1 - t * (beta - 2)) := by
  unfold e2
  rw [Finset.sum_congr rfl (fun p hp => by
    simp only [pairs, Finset.mem_filter, Finset.mem_univ, true_and] at hp
    exact minor_pair A hA beta t hp.ne)]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.sum_mul,
    ← Finset.sum_mul]
  simp only [Finset.sum_boole, nB, nN, nO]

/-- Every pair is of exactly one type: `B + Nn + O = #pairs`. -/
theorem pair_count (A : Matrix n n R) (hA : IsAdjR A) :
    nB A + nN A + nO A = (pairs n).card := by
  have h01 : (0 : R) ≠ 1 := zero_ne_one
  have h10 : (1 : R) ≠ 0 := one_ne_zero
  simp only [nB, nN, nO, Finset.card_filter]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib, Finset.card_eq_sum_ones]
  refine Finset.sum_congr rfl fun p _ => ?_
  rcases hA.1 p.1 p.2 with h1 | h1 <;> rcases hA.1 p.2 p.1 with h2 | h2 <;> simp [h1, h2, h01, h10]

/-- `#arcs = 2 B + O`. -/
theorem arcs_eq (A : Matrix n n R) (hA : IsAdjR A) : arcs A = 2 * nB A + nO A := by
  have h01 : (0 : R) ≠ 1 := zero_ne_one
  have h10 : (1 : R) ≠ 0 := one_ne_zero
  unfold arcs
  rw [Finset.card_filter, sum_prod_split]
  have hdiag : ∑ i : n, (if A (i, i).1 (i, i).2 = 1 then 1 else 0) = 0 := by
    apply Finset.sum_eq_zero; intro i _; simp [hA.2 i, h01]
  rw [hdiag, add_zero]
  simp only [nB, nO, Finset.card_filter, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun p _ => ?_
  simp only [Prod.fst_swap, Prod.snd_swap]
  rcases hA.1 p.1 p.2 with h1 | h1 <;> rcases hA.1 p.2 p.1 with h2 | h2 <;> simp [h1, h2, h01, h10]

end closed

/-! ### The arithmetic consequence (over ℝ) -/

/-- At a legal point with `beta > 2`, `e_2 > 0` forces more bidirected than non-adjacent pairs. -/
theorem bidir_gt_nonadj (A : Matrix n n ℝ) (hA : IsAdjR A) {beta t : ℝ} (hb : 2 < beta)
    (ht : 0 < t) (hleg : t * (beta - 2) > 1) (he : 0 < e2 (K (J ℝ - 1 - Aᵀ) beta t)) :
    nN A < nB A := by
  rw [e2_closed_form A hA] at he
  by_contra hc
  push Not at hc
  have hBN : (nB A : ℝ) ≤ (nN A : ℝ) := by exact_mod_cast hc
  have hB0 : (0 : ℝ) ≤ nB A := Nat.cast_nonneg _
  have hO0 : (0 : ℝ) ≤ nO A := Nat.cast_nonneg _
  have hc1 : 0 ≤ (beta - 1) * (beta + 1 + 2 * t) := by
    apply mul_nonneg <;> linarith
  have hc2 : 1 + 2 * t ≤ (beta - 1) * (beta + 1 + 2 * t) := by nlinarith
  have hO : (nO A : ℝ) * (1 - t * (beta - 2)) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos hO0 (by linarith)
  have h1 : (nB A : ℝ) * (1 + 2 * t) ≤ (nB A : ℝ) * ((beta - 1) * (beta + 1 + 2 * t)) :=
    mul_le_mul_of_nonneg_left hc2 hB0
  have h2 : (nB A : ℝ) * ((beta - 1) * (beta + 1 + 2 * t))
      ≤ (nN A : ℝ) * ((beta - 1) * (beta + 1 + 2 * t)) :=
    mul_le_mul_of_nonneg_right hBN hc1
  have h3 : (nN A : ℝ) * ((1 - beta) * (1 + beta + 2 * t))
      = -((nN A : ℝ) * ((beta - 1) * (beta + 1 + 2 * t))) := by ring
  linarith

/-- Same hypotheses: `2 #arcs ≥ m(m-1) + 2`, i.e. `#arcs ≥ m(m-1)/2 + 1`. -/
theorem arcs_lower_bound (A : Matrix n n ℝ) (hA : IsAdjR A) {beta t : ℝ} (hb : 2 < beta)
    (ht : 0 < t) (hleg : t * (beta - 2) > 1) (he : 0 < e2 (K (J ℝ - 1 - Aᵀ) beta t)) :
    Fintype.card n * (Fintype.card n - 1) + 2 ≤ 2 * arcs A := by
  have h := bidir_gt_nonadj A hA hb ht hleg he
  have h1 := pair_count A hA
  have h2 := arcs_eq A hA
  have h3 := two_mul_card_pairs (n := n)
  omega

/-- Spectral input: `e_2(K) > 0` for positive stable `K` with `m ≥ 2` (Mathlib's
`charpoly_coeff_eq_sum_minors`, iterated Gauss-Lucas, IVT; see `LemmaA2.sum_minors_two_pos`). -/
theorem e2_pos_of_posStable (Km : Matrix n n ℝ) (h : PosStable Km) (hm : 2 ≤ Fintype.card n) :
    0 < e2 Km := by
  rw [e2_eq_sum_minors]
  exact sum_minors_two_pos Km h hm

/-- **Edge bound, K form.** If `K(N, beta, t)` (with `N = J - I - Aᵀ`) is positive stable at a legal
point `t > 0`, `t (beta - 2) > 1`, and `m ≥ 2`, then `Nn < B` and `#arcs ≥ m(m-1)/2 + 1`. -/
theorem edge_bound_posStable (A : Matrix n n ℝ) (hA : IsAdjR A) {beta t : ℝ}
    (ht : 0 < t) (hleg : t * (beta - 2) > 1) (hm : 2 ≤ Fintype.card n)
    (hK : PosStable (K (J ℝ - 1 - Aᵀ) beta t)) :
    nN A < nB A ∧ Fintype.card n * (Fintype.card n - 1) + 2 ≤ 2 * arcs A := by
  have hb := beta_gt_two ht hleg
  have he := e2_pos_of_posStable _ hK hm
  exact ⟨bidir_gt_nonadj A hA hb ht hleg he, arcs_lower_bound A hA hb ht hleg he⟩

/-- **Edge bound, CTLN form.** If the full support is stable (`-I + W` Hurwitz) at a legal point
`0 < eps < 1`, `q = delta (1 - eps)/eps > 1`, on `m ≥ 2` nodes, then the graph has more bidirected
than non-adjacent pairs, and at least `m(m-1)/2 + 1` arcs. -/
theorem edge_bound_ctln (A : Matrix n n ℝ) (hA : IsAdj A) {eps delta : ℝ}
    (he0 : 0 < eps) (he1 : eps < 1) (hq : delta * (1 - eps) / eps > 1)
    (hm : 2 ≤ Fintype.card n) (hs : HurwitzStable (-1 + W A eps delta)) :
    nN A < nB A ∧ Fintype.card n * (Fintype.card n - 1) + 2 ≤ 2 * arcs A := by
  have ht : 0 < (1 - eps) / eps := div_pos (by linarith) he0
  have hleg := (legal_iff he0 he1).mp hq
  have hK := (hurwitz_iff_posStable A hA delta he0).mp hs
  exact edge_bound_posStable A hA ht hleg hm hK

end CtlnCert
