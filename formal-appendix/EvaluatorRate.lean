/-
SPDX-License-Identifier: Apache-2.0

Convert the evaluator fragment's natural step exponent to the square-root
logarithmic scale. The time to construct its input remains unaccounted for.
-/
import HalfSchedule

open Levels.Gen

namespace IntegerFactoring.HalfAppendix

theorem quotient_le_sqrt_rate (ρ : ℝ) (x l n : ℕ)
    (ha : AppRunningTime.PartA ρ x) (hn : 1 < n)
    (hnhi : n ≤ (K x (l + 1)) ^ 40) :
    ((n / 2 ^ (l / 2) : ℕ) : ℝ) ≤
      Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * (n : ℝ) / Real.sqrt (Real.log n) := by
  let r : ℕ := 2 ^ (l / 2)
  let D : ℝ := Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2)
  have hs : Real.sqrt (Real.log n) ≤ D * r := sqrt_log_le_half_level ρ x l n ha (by omega) hnhi
  have hD : 0 ≤ D := Real.sqrt_nonneg _
  have hden : 0 < Real.sqrt (Real.log n) := Real.sqrt_pos.2 (Real.log_pos (by exact_mod_cast hn))
  have hdiv : ((n / r : ℕ) : ℝ) * r ≤ n := by
    exact_mod_cast Nat.div_mul_le_self n r
  change ((n / r : ℕ) : ℝ) ≤ D * n / Real.sqrt (Real.log n)
  apply (le_div_iff₀ hden).2
  calc
    ((n / r : ℕ) : ℝ) * Real.sqrt (Real.log n) ≤ ((n / r : ℕ) : ℝ) * (D * r) :=
      mul_le_mul_of_nonneg_left hs (by positivity)
    _ = D * (((n / r : ℕ) : ℝ) * r) := by ring
    _ ≤ D * n := mul_le_mul_of_nonneg_left hdiv hD

theorem evaluator_step_exponent_sqrt (ρ : ℝ) (x l n D : ℕ)
    (ha : AppRunningTime.PartA ρ x) (hn : 1 < n)
    (hnhi : n ≤ (K x (l + 1)) ^ 40) :
    ((D * (n / 2 ^ (l / 2) + 1) : ℕ) : ℝ) ≤
      (D : ℝ) * Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) *
        (n : ℝ) / Real.sqrt (Real.log n) + D := by
  have h := quotient_le_sqrt_rate ρ x l n ha hn hnhi
  have h' := mul_le_mul_of_nonneg_left h (show 0 ≤ (D : ℝ) by positivity)
  have h'' : (D : ℝ) * ((n / 2 ^ (l / 2) : ℕ) : ℝ) ≤
      (D : ℝ) * Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * n / Real.sqrt (Real.log n) := by
    calc
      _ ≤ (D : ℝ) * (Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * n / Real.sqrt (Real.log n)) := h'
      _ = _ := by ring
  push_cast
  linarith

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.quotient_le_sqrt_rate
#print axioms IntegerFactoring.HalfAppendix.evaluator_step_exponent_sqrt
