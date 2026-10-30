/-
SPDX-License-Identifier: Apache-2.0

Real-valued interpretation of the exact natural-number time bound.
This file is independent of the factoring program and its modified sources.
-/
module

public import FactoringK6.RateArithmetic
public import Mathlib.Analysis.SpecialFunctions.Pow.Real
public import Mathlib.Data.Nat.Cast.Order.Field

@[expose] public section

namespace IntegerFactoring.K6RateInterpretation

open IntegerFactoring.K6RateArithmetic

/-- The positive integer denominator does not exceed the input length. -/
theorem denominator_le_input {n : ℕ} (hn : 2 ≤ n) : rateRoot n + 1 ≤ n := by
  have hroot := rateRoot_le_log n
  have hlog := Nat.log_lt_self 2 (show n ≠ 0 by omega)
  omega

/-- Rounding the binary-logarithm root upward also bounds the root of the
natural logarithm, with no additional constant factor. -/
theorem real_log_root_le_denominator {n : ℕ} (hn : 2 ≤ n) :
    (Real.log (n : ℝ)) ^ (1 / 30 : ℝ) ≤ (rateRoot n + 1 : ℕ) := by
  have hn0 : (0 : ℝ) < n := by exact_mod_cast (by omega : 0 < n)
  have hn1 : (1 : ℝ) < n := by exact_mod_cast (by omega : 1 < n)
  have hroot : Nat.log 2 n + 1 ≤ (rateRoot n + 1) ^ 30 := by
    have h := Nat.lt_pow_nthRoot_add_one (by norm_num : 30 ≠ 0) (Nat.log 2 n)
    change Nat.log 2 n < (rateRoot n + 1) ^ 30 at h
    omega
  have hpow : (n : ℝ) ≤ (2 : ℝ) ^ (Nat.log 2 n + 1) := by
    exact_mod_cast (Nat.lt_pow_succ_log_self (by norm_num : 1 < 2) n).le
  have hlog2 : Real.log 2 ≤ 1 := by
    have h := Real.log_le_sub_one_of_pos (by norm_num : (0 : ℝ) < 2)
    norm_num at h ⊢
    exact h
  have hlog : Real.log (n : ℝ) ≤ ((rateRoot n + 1 : ℕ) : ℝ) ^ 30 := by
    calc
      Real.log (n : ℝ) ≤ Real.log ((2 : ℝ) ^ (Nat.log 2 n + 1)) :=
        Real.log_le_log hn0 hpow
      _ = ((Nat.log 2 n + 1 : ℕ) : ℝ) * Real.log 2 := Real.log_pow _ _
      _ ≤ ((Nat.log 2 n + 1 : ℕ) : ℝ) := by
        have h := mul_le_mul_of_nonneg_left hlog2
          (show (0 : ℝ) ≤ ((Nat.log 2 n + 1 : ℕ) : ℝ) by positivity)
        simpa only [mul_one] using h
      _ ≤ ((rateRoot n + 1 : ℕ) : ℝ) ^ 30 := by exact_mod_cast hroot
  calc
    (Real.log (n : ℝ)) ^ (1 / 30 : ℝ)
        ≤ (((rateRoot n + 1 : ℕ) : ℝ) ^ 30) ^ (1 / 30 : ℝ) :=
      Real.rpow_le_rpow (Real.log_pos hn1).le hlog (by norm_num)
    _ = ((rateRoot n + 1 : ℕ) : ℝ) := by
      simpa using Real.pow_rpow_inv_natCast
        (show (0 : ℝ) ≤ ((rateRoot n + 1 : ℕ) : ℝ) by positivity)
        (by norm_num : (30 : ℕ) ≠ 0)

/-- The exact exponent, including its additive one, is bounded by the usual
real-valued logarithmic-power expression. -/
theorem natural_exponent_le_real {n : ℕ} (hn : 2 ≤ n) :
    ((n / (rateRoot n + 1) + 1 : ℕ) : ℝ) ≤
      2 * (n : ℝ) / (Real.log (n : ℝ)) ^ (1 / 30 : ℝ) := by
  let D : ℝ := ((rateRoot n + 1 : ℕ) : ℝ)
  let H : ℝ := (Real.log (n : ℝ)) ^ (1 / 30 : ℝ)
  have hD : 0 < D := by dsimp [D]; positivity
  have hn1 : (1 : ℝ) < n := by exact_mod_cast (by omega : 1 < n)
  have hH : 0 < H := Real.rpow_pos_of_pos (Real.log_pos hn1) _
  have hHD : H ≤ D := real_log_root_le_denominator hn
  have hDn : D ≤ n := by
    dsimp [D]
    exact_mod_cast denominator_le_input hn
  have hq : ((n / (rateRoot n + 1) : ℕ) : ℝ) ≤ (n : ℝ) / D :=
    Nat.cast_div_le
  have hone : (1 : ℝ) ≤ (n : ℝ) / D := (le_div_iff₀ hD).mpr (by simpa using hDn)
  calc
    ((n / (rateRoot n + 1) + 1 : ℕ) : ℝ)
        = ((n / (rateRoot n + 1) : ℕ) : ℝ) + 1 := by push_cast; rfl
    _ ≤ (n : ℝ) / D + (n : ℝ) / D := add_le_add hq hone
    _ = 2 * (n : ℝ) / D := by ring
    _ ≤ 2 * (n : ℝ) / H := div_le_div_of_nonneg_left (by positivity) hH hHD

/-- An explicit constant for the real exponential bound. -/
theorem real_rate_bound (f : ℕ → ℕ) (c : ℕ)
    (hf : ∀ n, f n ≤ c * 2 ^ (c * (n / (rateRoot n + 1) + 1))) :
    ∀ n, 2 ≤ n → (f n : ℝ) ≤
      (2 * (c : ℝ) + 1) * Real.exp
        ((2 * (c : ℝ) + 1) * (n : ℝ) / (Real.log (n : ℝ)) ^ (1 / 30 : ℝ)) := by
  intro n hn
  let H : ℝ := (Real.log (n : ℝ)) ^ (1 / 30 : ℝ)
  let q : ℕ := n / (rateRoot n + 1) + 1
  have hc : (0 : ℝ) ≤ c := by positivity
  have hn1 : (1 : ℝ) < n := by exact_mod_cast (by omega : 1 < n)
  have hH : 0 < H := Real.rpow_pos_of_pos (Real.log_pos hn1) _
  have hq : (q : ℝ) ≤ 2 * (n : ℝ) / H := natural_exponent_le_real hn
  have hlog2 : Real.log 2 ≤ 1 := by
    have h := Real.log_le_sub_one_of_pos (by norm_num : (0 : ℝ) < 2)
    norm_num at h ⊢
    exact h
  have hlogq : Real.log 2 * ((c * q : ℕ) : ℝ) ≤
      (2 * (c : ℝ) + 1) * (n : ℝ) / H := by
    calc
      Real.log 2 * ((c * q : ℕ) : ℝ) ≤ ((c * q : ℕ) : ℝ) := by
        have h := mul_le_mul_of_nonneg_right hlog2
          (show (0 : ℝ) ≤ ((c * q : ℕ) : ℝ) by positivity)
        simpa only [one_mul] using h
      _ = (c : ℝ) * (q : ℝ) := by push_cast; rfl
      _ ≤ (c : ℝ) * (2 * (n : ℝ) / H) := mul_le_mul_of_nonneg_left hq hc
      _ ≤ (2 * (c : ℝ) + 1) * (n : ℝ) / H := by
        have hnonneg : (0 : ℝ) ≤ (n : ℝ) / H := by positivity
        calc
          (c : ℝ) * (2 * (n : ℝ) / H) = (2 * (c : ℝ)) * ((n : ℝ) / H) := by ring
          _ ≤ (2 * (c : ℝ) + 1) * ((n : ℝ) / H) :=
            mul_le_mul_of_nonneg_right (by linarith) hnonneg
          _ = (2 * (c : ℝ) + 1) * (n : ℝ) / H := by ring
  have hpow : (2 : ℝ) ^ (c * q) = Real.exp (Real.log 2 * ((c * q : ℕ) : ℝ)) := by
    rw [mul_comm (Real.log 2) ((c * q : ℕ) : ℝ), Real.exp_nat_mul,
      Real.exp_log (by norm_num : (0 : ℝ) < 2)]
  calc
    (f n : ℝ) ≤ (c : ℝ) * (2 : ℝ) ^ (c * q) := by exact_mod_cast hf n
    _ = (c : ℝ) * Real.exp (Real.log 2 * ((c * q : ℕ) : ℝ)) := by rw [hpow]
    _ ≤ (c : ℝ) * Real.exp ((2 * (c : ℝ) + 1) * (n : ℝ) / H) :=
      mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr hlogq) hc
    _ ≤ (2 * (c : ℝ) + 1) * Real.exp ((2 * (c : ℝ) + 1) * (n : ℝ) / H) :=
      mul_le_mul_of_nonneg_right (by linarith) (Real.exp_pos _).le

/-- The integer time bound entails the conventional real asymptotic bound,
uniformly for every input length at least two. -/
theorem exists_real_rate_bound (f : ℕ → ℕ) (c : ℕ)
    (hf : ∀ n, f n ≤ c * 2 ^ (c * (n / (rateRoot n + 1) + 1))) :
    ∃ C : ℝ, 0 < C ∧ ∀ n : ℕ, 2 ≤ n →
      (f n : ℝ) ≤ C * Real.exp (C * (n : ℝ) / (Real.log (n : ℝ)) ^ (1 / 30 : ℝ)) := by
  refine ⟨2 * (c : ℝ) + 1, by positivity, real_rate_bound f c hf⟩

/-- The new saving exponent is strictly larger than the degree-sixteen
reference rate even at its limiting contraction 11/12. -/
theorem baseline_log_rate_lt_one_thirtieth :
    Real.log (12 / 11 : ℝ) / Real.log 16 < (1 / 30 : ℝ) := by
  have hnum : (12 : ℝ) ^ 30 < 16 * (11 : ℝ) ^ 30 := by
    exact_mod_cast baseline_rate_comparison
  have hlog := Real.log_lt_log (by positivity : (0 : ℝ) < (12 : ℝ) ^ 30) hnum
  rw [Real.log_pow, Real.log_mul (by norm_num) (by positivity), Real.log_pow] at hlog
  rw [div_lt_iff₀ (Real.log_pos (by norm_num : (1 : ℝ) < 16)),
    Real.log_div (by norm_num) (by norm_num)]
  norm_num at hlog ⊢
  linarith

end IntegerFactoring.K6RateInterpretation

#print axioms IntegerFactoring.K6RateInterpretation.real_log_root_le_denominator
#print axioms IntegerFactoring.K6RateInterpretation.exists_real_rate_bound
#print axioms IntegerFactoring.K6RateInterpretation.baseline_log_rate_lt_one_thirtieth
