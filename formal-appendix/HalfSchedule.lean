/-
SPDX-License-Identifier: Apache-2.0

The square-root logarithmic conversion for the Appendix B level schedule.
This is arithmetic only and does not establish a factoring machine bound.
-/
import ExplicitPad

open Levels.Gen

namespace IntegerFactoring.HalfAppendix

theorem two_pow_le_twice_half_square (l : ℕ) :
    2 ^ l ≤ 2 * (2 ^ (l / 2)) ^ 2 := by
  have h : l ≤ 2 * (l / 2) + 1 := by omega
  calc
    2 ^ l ≤ 2 ^ (2 * (l / 2) + 1) := Nat.pow_le_pow_right (by norm_num) h
    _ = 2 * (2 ^ (l / 2)) ^ 2 := by
      rw [show 2 * (l / 2) + 1 = (l / 2) * 2 + 1 by omega, pow_add, pow_mul]
      ring

theorem two_pow_le_K (x l : ℕ) (hx : 1 ≤ x) : 2 ^ l ≤ K x l := by
    induction l with
    | zero =>
        simp only [pow_zero, K]
        exact Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
    | succ d ih =>
        rw [pow_succ]
        have htwo : 2 ≤ K x d := K_two x d hx
        have hbase : 2 * 2 ^ d ≤ K x d * K x d := Nat.mul_le_mul htwo ih
        have hpow : 1 ≤ 2 ^ (3 * d + x) := Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
        have hm : K x d ≤ 2 ^ (3 * d + x) * K x d := by nlinarith
        simp only [K]
        nlinarith

theorem half_level_le_K (x l : ℕ) (hx : 1 ≤ x) : 2 ^ (l / 2) ≤ K x l :=
  (Nat.pow_le_pow_right (by norm_num) (by omega : l / 2 ≤ l)).trans (two_pow_le_K x l hx)

/-- If the next scheduled level exceeds n^(1/40), the integer parameter
2^(floor(l/2)) controls the logarithm of n. -/
theorem log_le_half_level (ρ : ℝ) (x l n : ℕ)
    (ha : AppRunningTime.PartA ρ x) (hn : 0 < n)
    (hmax : n ≤ (K x (l + 1)) ^ 40) :
    Real.logb 2 n ≤ 160 * (3 * (x : ℝ) + 6) * ((2 ^ (l / 2) : ℕ) : ℝ) ^ 2 := by
  have hlog : Real.logb 2 n ≤ Real.logb 2 ((K x (l + 1) : ℝ) ^ 40) := by
    apply Real.logb_le_logb_of_le (by norm_num) (by exact_mod_cast hn)
    exact_mod_cast hmax
  rw [Real.logb_pow, ha.logK (l + 1)] at hlog
  have hp : (2 : ℝ) ^ l ≤ 2 * ((2 : ℝ) ^ (l / 2)) ^ 2 := by
    exact_mod_cast two_pow_le_twice_half_square l
  have hx : 0 ≤ (x : ℝ) := by positivity
  have hl : 0 ≤ (l : ℝ) := by positivity
  have hmult := mul_le_mul_of_nonneg_left hp (show 0 ≤ 3 * (x : ℝ) + 6 by positivity)
  push_cast at hlog ⊢
  rw [pow_succ] at hlog
  nlinarith

/-- The same bound in natural logarithms, with the square root exposed. -/
theorem sqrt_log_le_half_level (ρ : ℝ) (x l n : ℕ)
    (ha : AppRunningTime.PartA ρ x) (hn : 0 < n)
    (hmax : n ≤ (K x (l + 1)) ^ 40) :
    Real.sqrt (Real.log n) ≤
      Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * ((2 ^ (l / 2) : ℕ) : ℝ) := by
  have hlog := log_le_half_level ρ x l n ha hn hmax
  have htwo : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hnlog : Real.log (n : ℝ) ≤
      (160 * (3 * (x : ℝ) + 6) * Real.log 2) * ((2 ^ (l / 2) : ℕ) : ℝ) ^ 2 := by
    rw [Real.logb] at hlog
    have := (div_le_iff₀ htwo).mp hlog
    nlinarith
  calc
    Real.sqrt (Real.log n) ≤ Real.sqrt
        ((160 * (3 * (x : ℝ) + 6) * Real.log 2) * ((2 ^ (l / 2) : ℕ) : ℝ) ^ 2) :=
      Real.sqrt_le_sqrt hnlog
    _ = Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * ((2 ^ (l / 2) : ℕ) : ℝ) := by
      rw [Real.sqrt_mul (by positivity), Real.sqrt_sq_eq_abs, abs_of_nonneg (by positivity)]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.two_pow_le_twice_half_square
#print axioms IntegerFactoring.HalfAppendix.two_pow_le_K
#print axioms IntegerFactoring.HalfAppendix.half_level_le_K
#print axioms IntegerFactoring.HalfAppendix.log_le_half_level
#print axioms IntegerFactoring.HalfAppendix.sqrt_log_le_half_level
