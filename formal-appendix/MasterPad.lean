/-
SPDX-License-Identifier: Apache-2.0

One explicit pad budget for grouping search, top-table writing, and the
Appendix evaluator. This certifies the parameters, not the row query.
-/
import TopParamsPoly
import SelectedPadCost
import GroupingSearchCost

open Groundwork Levels.Gen Factoring.Aux IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

def masterExponent (x C n : ℕ) : ℕ :=
  let l := selectedLevel x n
  let P := K x l
  let L := topWidth (sOf n) P + 1
  let B := topEntryBits (sOf n) P
  3 * P ^ 9 + (P - 1) * L + B + padExponent P L B C (2 ^ (l / 2))

theorem masterExponent_unary (x C : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : Unary fun a => masterExponent x C (n a) := by
  have hK := selected_K_unary x hx hn
  have hg := Unary.mul (Unary.const 3) (unary_pow_const hK 9)
  have ht := selected_tableExponent_unary x hx hn
  have hb := selected_topEntryBits_unary x hx hn
  have he := selected_padExponent_unary x C hx hn
  exact Unary.add (Unary.add (Unary.add hg ht) hb) he

theorem master_pad_suffices (x C n : ℕ) (pad : List Unit)
    (hpad : 2 ^ masterExponent x C n ≤ pad.length) :
    let l := selectedLevel x n
    let P := K x l
    let L := topWidth (sOf n) P + 1
    let B := topEntryBits (sOf n) P
    P ^ (3 * P ^ 8) ≤ pad.length ∧
      2 ^ ((P - 1) * L) ≤ pad.length ∧
      2 ^ B ≤ pad.length ∧
      2 ^ padExponent P L B C (2 ^ (l / 2)) ≤ pad.length := by
  dsimp only
  have hgroup : 3 * (K x (selectedLevel x n)) ^ 9 ≤ masterExponent x C n := by
    unfold masterExponent; dsimp; omega
  have htable : (K x (selectedLevel x n) - 1) *
      (topWidth (sOf n) (K x (selectedLevel x n)) + 1) ≤ masterExponent x C n := by
    unfold masterExponent; dsimp; omega
  have heval : padExponent (K x (selectedLevel x n))
      (topWidth (sOf n) (K x (selectedLevel x n)) + 1)
      (topEntryBits (sOf n) (K x (selectedLevel x n))) C (2 ^ (selectedLevel x n / 2)) ≤
      masterExponent x C n := by simp [masterExponent]
  have hbits : topEntryBits (sOf n) (K x (selectedLevel x n)) ≤ masterExponent x C n := by
    unfold masterExponent; dsimp; omega
  exact ⟨(grouping_count_le_pad _).trans ((Nat.pow_le_pow_right (by norm_num) hgroup).trans hpad),
    (Nat.pow_le_pow_right (by norm_num) htable).trans hpad,
    (Nat.pow_le_pow_right (by norm_num) hbits).trans hpad,
    (Nat.pow_le_pow_right (by norm_num) heval).trans hpad⟩

theorem masterExponent_cost (ρ : ℝ) (x C n : ℕ) (hx : 1 ≤ x)
    (ha : AppRunningTime.PartA ρ x) (hn : (K x 0) ^ 40 ≤ n) :
    2 ^ (selectedLevel x n / 2) * masterExponent x C n ≤ (20004 + C * 500) * n := by
  let l := selectedLevel x n
  let P := K x l
  let r := 2 ^ (l / 2)
  have hs := selectedLevel_spec x n hx hn
  have hP : 2 ≤ P := K_two x l hx
  have hr : r ≤ P := half_level_le_K x l hx
  have hgroup : r * (3 * P ^ 9) ≤ 3 * n := grouping_pad_cost P n r hP hs.1 hr
  have hwritten := selected_written_exponents P n r hP hs.1 hr
  have htable := hwritten.1
  have hbits := hwritten.2
  have heval := scheduled_pad_cost ρ x l n C hx ha hs.1 (le_of_lt hs.2)
  change r * padExponent P (topWidth (sOf n) P + 1) (topEntryBits (sOf n) P) C r ≤ _ at heval
  change r * (3 * P ^ 9 + (P - 1) * (topWidth (sOf n) P + 1) + topEntryBits (sOf n) P +
    padExponent P (topWidth (sOf n) P + 1) (topEntryBits (sOf n) P) C r) ≤ _
  nlinarith

theorem masterExponent_sqrt_cost (ρ : ℝ) (x C n : ℕ) (hx : 1 ≤ x)
    (ha : AppRunningTime.PartA ρ x) (hn : 1 < n) (hlarge : (K x 0) ^ 40 ≤ n) :
    (masterExponent x C n : ℝ) ≤ ((20004 + C * 500 : ℕ) : ℝ) *
      Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * (n : ℝ) / Real.sqrt (Real.log n) := by
  let l := selectedLevel x n
  let r := 2 ^ (l / 2)
  let D : ℝ := Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2)
  have hs := selectedLevel_spec x n hx hlarge
  have hq : (r : ℝ) * masterExponent x C n ≤ ((20004 + C * 500 : ℕ) : ℝ) * n := by
    exact_mod_cast masterExponent_cost ρ x C n hx ha hlarge
  have hsqrt : Real.sqrt (Real.log n) ≤ D * r := sqrt_log_le_half_level ρ x l n ha (by omega) (le_of_lt hs.2)
  have hD : 0 ≤ D := Real.sqrt_nonneg _
  have hden : 0 < Real.sqrt (Real.log n) := Real.sqrt_pos.2 (Real.log_pos (by exact_mod_cast hn))
  apply (le_div_iff₀ hden).2
  calc
    (masterExponent x C n : ℝ) * Real.sqrt (Real.log n) ≤ masterExponent x C n * (D * r) :=
      mul_le_mul_of_nonneg_left hsqrt (by positivity)
    _ = D * ((r : ℝ) * masterExponent x C n) := by ring
    _ ≤ D * (((20004 + C * 500 : ℕ) : ℝ) * n) := mul_le_mul_of_nonneg_left hq hD
    _ = ((20004 + C * 500 : ℕ) : ℝ) * D * n := by ring

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.masterExponent_unary
#print axioms IntegerFactoring.HalfAppendix.master_pad_suffices
#print axioms IntegerFactoring.HalfAppendix.masterExponent_cost
#print axioms IntegerFactoring.HalfAppendix.masterExponent_sqrt_cost
