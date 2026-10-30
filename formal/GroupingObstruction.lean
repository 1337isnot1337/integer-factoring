/-
SPDX-License-Identifier: Apache-2.0

A connected unit path obstructs grouping below a quadratic vertex count.
Uses the exact grpSize, grpY, and totU definitions of the pinned upstream.
This is a lower bound for universal grouping, not for factoring algorithms.
-/
module

public import FactoringAux.Grouping

@[expose] public section

namespace Factoring.Improvements

open Finset Factoring.Aux

/-- A connected support must have a crossing edge in any balanced grouping
with more than two groups. One such edge forces this integral lower bound. -/
theorem connected_grouping_lower_bound {g K : ℕ} (hg : 0 < g) (hK : 2 < K)
    (c : Fin g → Fin g → ℕ) (P D : ℕ)
    (hconnected : ∀ χ : Fin g → Fin K,
      (∀ a b, c a b ≠ 0 → χ a = χ b) → ∀ a b, χ a = χ b)
    (χ : Fin g → Fin K)
    (hsize : ∀ γ, (K : ℤ) * grpSize χ γ ≤ 2 * g)
    (hpair : ∀ γ γ', γ ≠ γ' → (D : ℤ) * (K : ℤ) ^ 2 * grpY c χ γ γ' ≤
      D * totU c + g * P) :
    (D : ℤ) * (K : ℤ) ^ 2 ≤ D * totU c + g * P := by
  by_contra h
  have hsame : ∀ a b, c a b ≠ 0 → χ a = χ b := by
    intro a b hab
    by_contra hne
    have hy : 1 ≤ grpY c χ (χ a) (χ b) := by
      calc
        1 ≤ (c a b : ℤ) := by exact_mod_cast (Nat.one_le_iff_ne_zero.mpr hab)
        _ = (c a b : ℤ) * ind2 χ (χ a) (χ b) (a, b) := by simp [ind2]
        _ ≤ grpY c χ (χ a) (χ b) := by
          unfold grpY
          exact Finset.single_le_sum (fun q _ =>
            mul_nonneg (Nat.cast_nonneg (c q.1 q.2))
              (ind2_nonneg χ (χ a) (χ b) q)) (mem_univ (a, b))
    have hmul := mul_le_mul_of_nonneg_left hy (by positivity : (0 : ℤ) ≤ D * (K : ℤ) ^ 2)
    have hbound := hpair (χ a) (χ b) hne
    simp only [mul_one] at hmul
    exact h (le_trans hmul hbound)
  have hconstant := hconnected χ hsame
  let a : Fin g := ⟨0, hg⟩
  have hfull : grpSize χ (χ a) = g := by
    have hone (b : Fin g) : (if χ b = χ a then (1 : ℤ) else 0) = 1 :=
      ite_eq_left (hconstant b a)
    simp only [grpSize, hone, sum_const, card_univ, Fintype.card_fin,
      nsmul_eq_mul, mul_one]
  have hbad := hsize (χ a)
  rw [hfull] at hbad
  have hK' : (2 : ℤ) < K := by exact_mod_cast hK
  have hg' : (0 : ℤ) < g := by exact_mod_cast hg
  nlinarith

/-- Unit-weight path on m+1 vertices. The two directed indicators are disjoint. -/
def pathWeight (m : ℕ) (a b : Fin (m + 1)) : ℕ :=
  (if a.val + 1 = b.val then 1 else 0) + (if b.val + 1 = a.val then 1 else 0)

theorem pathWeight_diag (m : ℕ) (a : Fin (m + 1)) : pathWeight m a a = 0 := by
  simp [pathWeight]

theorem pathWeight_symm (m : ℕ) (a b : Fin (m + 1)) :
    pathWeight m a b = pathWeight m b a := by
  simp only [pathWeight, add_comm]

theorem pathWeight_row_le_two (m : ℕ) (a : Fin (m + 1)) :
    ∑ b, pathWeight m a b ≤ 2 := by
  have hf : (∑ b : Fin (m + 1), if a.val + 1 = b.val then (1 : ℕ) else 0) ≤ 1 := by
    simp only [Finset.sum_boole, Nat.cast_id]
    rw [Finset.card_le_one]
    intro b hb c hc
    simp only [mem_filter, mem_univ, true_and] at hb hc
    apply Fin.ext
    omega
  have hb : (∑ b : Fin (m + 1), if b.val + 1 = a.val then (1 : ℕ) else 0) ≤ 1 := by
    simp only [Finset.sum_boole, Nat.cast_id]
    rw [Finset.card_le_one]
    intro b hb c hc
    simp only [mem_filter, mem_univ, true_and] at hb hc
    apply Fin.ext
    omega
  simpa only [pathWeight, Finset.sum_add_distrib] using Nat.add_le_add hf hb

theorem pathWeight_total (m : ℕ) : totU (pathWeight m) = 2 * m := by
  have hf : (∑ a : Fin (m + 1), ∑ b : Fin (m + 1),
      if a.val + 1 = b.val then (1 : ℤ) else 0) = m := by
    rw [Finset.sum_comm, Fin.sum_univ_succ]
    have hsucc (b : Fin m) : (∑ a : Fin (m + 1),
        if a.val + 1 = b.succ.val then (1 : ℤ) else 0) = 1 := by
      have heq (a : Fin (m + 1)) : a.val + 1 = b.succ.val ↔ a = b.castSucc := by
        simp [Fin.ext_iff]
      simp_rw [heq]
      simp
    have hzero : (∑ a : Fin (m + 1), if a.val + 1 = (0 : Fin (m + 1)).val
        then (1 : ℤ) else 0) = 0 := by simp
    rw [hzero]
    simp only [hsucc, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one,
      zero_add]
  have hr : (∑ a : Fin (m + 1), ∑ b : Fin (m + 1),
      if b.val + 1 = a.val then (1 : ℤ) else 0) = m := by
    rw [Finset.sum_comm]
    exact hf
  simp only [totU, pathWeight, Nat.cast_add, Nat.cast_ite, Nat.cast_one,
    Nat.cast_zero, Fintype.sum_prod_type, Finset.sum_add_distrib]
  rw [hf, hr]
  ring

/-- Every valid grouping of the path satisfies
D K² ≤ 2D(g-1)+2g. Thus universal grouping requires Ω(K²) vertices. -/
theorem path_grouping_lower_bound {m K : ℕ} (hK : 2 < K) (D : ℕ)
    (χ : Fin (m + 1) → Fin K)
    (hsize : ∀ γ, (K : ℤ) * grpSize χ γ ≤ 2 * (m + 1))
    (hpair : ∀ γ γ', γ ≠ γ' → (D : ℤ) * (K : ℤ) ^ 2 * grpY (pathWeight m) χ γ γ' ≤
      D * totU (pathWeight m) + (m + 1) * 2) :
    (D : ℤ) * ((K : ℤ) ^ 2 + 2) ≤ 2 * (D + 1) * (m + 1) := by
  have hc : ∀ ψ : Fin (m + 1) → Fin K,
      (∀ a b, pathWeight m a b ≠ 0 → ψ a = ψ b) → ∀ a b, ψ a = ψ b := by
    intro ψ hedge
    have hzero : ∀ a : Fin (m + 1), ψ a = ψ 0 := by
      intro a
      refine Fin.induction ?_ (fun i hi => ?_) a
      · rfl
      · exact (hedge i.castSucc i.succ (by simp [pathWeight])).symm.trans hi
    intro a b
    exact (hzero a).trans (hzero b).symm
  have h := connected_grouping_lower_bound (Nat.zero_lt_succ m) hK
    (pathWeight m) 2 D hc χ (by simpa using hsize) (by simpa using hpair)
  rw [pathWeight_total] at h
  push_cast at h ⊢
  nlinarith

#print axioms connected_grouping_lower_bound
#print axioms pathWeight_diag
#print axioms pathWeight_symm
#print axioms pathWeight_row_le_two
#print axioms pathWeight_total
#print axioms path_grouping_lower_bound

end Factoring.Improvements
