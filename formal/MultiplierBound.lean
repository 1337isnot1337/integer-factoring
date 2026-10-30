/-
SPDX-License-Identifier: Apache-2.0

Bounds for the proposed level multipliers. These statements do not construct
or verify the modified factoring evaluator.
-/
module

public import FactoringAux.Moduli
public import Mathlib.Tactic.NormNum
public import Mathlib.Tactic.Linarith

@[expose] public section

namespace Factoring.Improvements

open Factoring.Aux Finset

/-- The redundant choices and scaled rectangles at one level. -/
def levelMultiplier (g t r : ℕ) : ℕ := (t + 1) ^ (g ^ 4 - r) * 8 ^ (t * r)

/-- A bound independent of the number of actual rectangles. -/
theorem levelMultiplier_le (g t r : ℕ) (hr : r ≤ g ^ 4) :
    levelMultiplier g t r ≤ 2 ^ (4 * g ^ 4 * t) := by
  have ht : t + 1 ≤ 2 ^ t := Nat.lt_two_pow_self (n := t)
  have h1 : (t + 1) ^ (g ^ 4 - r) ≤ (2 ^ t) ^ (g ^ 4 - r) :=
    Nat.pow_le_pow_left ht _
  unfold levelMultiplier
  calc
    (t + 1) ^ (g ^ 4 - r) * 8 ^ (t * r) ≤
        (2 ^ t) ^ (g ^ 4 - r) * 8 ^ (t * r) := Nat.mul_le_mul_right _ h1
    _ = 2 ^ (t * (g ^ 4 - r) + 3 * (t * r)) := by
      rw [← pow_mul, show (8 : ℕ) = 2 ^ 3 by norm_num, ← pow_mul, ← pow_add]
    _ ≤ 2 ^ (4 * g ^ 4 * t) := Nat.pow_le_pow_right (by norm_num) (by nlinarith [Nat.sub_add_cancel hr])

/-- The product of level multipliers has the stated uniform binary envelope. -/
theorem multipliers_le {ι : Type*} (s : Finset ι) (g t r : ι → ℕ)
    (hr : ∀ i ∈ s, r i ≤ g i ^ 4) :
    (∏ i ∈ s, levelMultiplier (g i) (t i) (r i)) ≤
      2 ^ (∑ i ∈ s, 4 * g i ^ 4 * t i) := by
  calc
    (∏ i ∈ s, levelMultiplier (g i) (t i) (r i)) ≤
        ∏ i ∈ s, 2 ^ (4 * g i ^ 4 * t i) := by
          exact Finset.prod_le_prod fun i hi => levelMultiplier_le _ _ _ (hr i hi)
    _ = 2 ^ (∑ i ∈ s, 4 * g i ^ 4 * t i) := by
          rw [Finset.prod_pow_eq_pow_sum]

/-- Scaled level values remain positive; no modular invertibility is needed. -/
theorem levelMultiplier_pos (g t r : ℕ) : 0 < levelMultiplier g t r := by
  unfold levelMultiplier
  positivity

theorem multipliers_pos {ι : Type*} (s : Finset ι) (g t r : ι → ℕ) :
    0 < ∏ i ∈ s, levelMultiplier (g i) (t i) (r i) := by
  exact Nat.pos_iff_ne_zero.mpr (Finset.prod_ne_zero_iff.mpr
    fun i _ => Nat.ne_of_gt (levelMultiplier_pos _ _ _))

/-- An integer count multiplied by a positive factor can still be detected
by the source's small-moduli test, even if the factor shares primes with a
tested modulus. -/
theorem scaled_count_detected {M H c : ℕ} (hM : 0 < M) (hH : 0 < H) (hc : 0 < c)
    (hbound : (2 * M + 1) * (H * c) < 4 ^ M) :
    ∃ m, 2 ≤ m ∧ m ≤ 2 * M ∧ ¬m ∣ H * c := by
  by_contra h
  push Not at h
  have hall : ∀ m, 1 ≤ m → m ≤ 2 * M → m ∣ H * c := by
    intro m hm hmM
    by_cases hm2 : 2 ≤ m
    · exact h m hm2 hmM
    · have : m = 1 := by omega
      subst m
      exact one_dvd _
  have hz : H * c = 0 := Factoring.Aux.eq_zero_of_all_dvd hM hall hbound
  have hp : 0 < H * c := Nat.mul_pos hH hc
  omega

/-- The explicit modulus cutoff used by the proposed integration suffices
when the multiplier and count have the stated binary-size bounds. -/
theorem scaled_count_detected_explicit {H c n S : ℕ} (hH : 0 < H)
    (hc : 0 < c) (hHbound : H ≤ 2 ^ S) (hcbound : c < 2 ^ n) :
    ∃ m, 2 ≤ m ∧ m ≤ 2 * (n + S + 4) ∧ ¬m ∣ H * c := by
  let M := n + S + 4
  have hsmall : H * c < 2 ^ (n + S) := by
    rw [pow_add]
    calc H * c ≤ 2 ^ S * c := Nat.mul_le_mul_right _ hHbound
      _ < 2 ^ S * 2 ^ n := Nat.mul_lt_mul_of_pos_left hcbound (by positivity)
      _ = 2 ^ n * 2 ^ S := by rw [mul_comm]
  have hpower : 2 * M + 1 ≤ 2 ^ (n + S + 8) := by
    have h := Nat.lt_two_pow_self (n := n + S)
    rw [pow_add]
    norm_num
    dsimp [M]
    omega
  apply scaled_count_detected (M := M) (by dsimp [M]; omega) hH hc
  calc (2 * M + 1) * (H * c) < (2 * M + 1) * 2 ^ (n + S) :=
          Nat.mul_lt_mul_of_pos_left hsmall (by omega)
    _ ≤ 2 ^ (n + S + 8) * 2 ^ (n + S) := Nat.mul_le_mul_right _ hpower
    _ = 4 ^ M := by
          rw [← pow_add, show (4 : ℕ) = 2 ^ 2 by norm_num, ← pow_mul]
          congr 1
          dsimp [M]
          omega

/-- The proposed product of level multipliers has a concrete finite cutoff
for detecting a positive count. The rectangle-count condition and the
count-size bound are explicit hypotheses, not conclusions about an evaluator. -/
theorem multiplier_product_detected {ι : Type*} (s : Finset ι) (g t r : ι → ℕ)
    (hr : ∀ i ∈ s, r i ≤ g i ^ 4) {c n : ℕ} (hc : 0 < c) (hcbound : c < 2 ^ n) :
    ∃ m, 2 ≤ m ∧ m ≤ 2 * (n + (∑ i ∈ s, 4 * g i ^ 4 * t i) + 4) ∧
      ¬m ∣ (∏ i ∈ s, levelMultiplier (g i) (t i) (r i)) * c := by
  exact scaled_count_detected_explicit (multipliers_pos s g t r) hc (multipliers_le s g t r hr) hcbound

#print axioms levelMultiplier_le
#print axioms multipliers_le
#print axioms levelMultiplier_pos
#print axioms multipliers_pos
#print axioms scaled_count_detected
#print axioms scaled_count_detected_explicit
#print axioms multiplier_product_detected

end Factoring.Improvements
