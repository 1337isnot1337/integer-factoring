/-
SPDX-License-Identifier: Apache-2.0

Arithmetic for the explicit pad in the half-rate composition. This file
does not construct the factoring machine or prove the composition.
-/
import Mathlib

namespace IntegerFactoring.HalfPad

def padExponent (P lam B C r : ℕ) : ℕ :=
  P ^ 2 * lam / r + 1 + C * P ^ 2 * (Nat.log 2 (P * (lam + B)) + 1)

theorem half_power_le (l : ℕ) :
    (1 / Real.sqrt 2 : ℝ) ^ l ≤ 1 / (2 ^ (l / 2) : ℝ) := by
  have hs : 0 < Real.sqrt 2 := Real.sqrt_pos.2 (by norm_num)
  have hs2 : Real.sqrt 2 ^ 2 = (2 : ℝ) := Real.sq_sqrt (by norm_num)
  have hp : (1 / Real.sqrt 2 : ℝ) ^ 2 = 1 / 2 := by
    field_simp
    nlinarith
  have hbound : (1 / Real.sqrt 2 : ℝ) ≤ 1 := by
    rw [div_le_one hs]
    nlinarith
  have hr : (1 / Real.sqrt 2 : ℝ) ^ (l % 2) ≤ 1 :=
    pow_le_one₀ (by positivity) hbound
  have hdecomp : l = 2 * (l / 2) + l % 2 := by omega
  calc
    (1 / Real.sqrt 2 : ℝ) ^ l = ((1 / Real.sqrt 2 : ℝ) ^ 2) ^ (l / 2) *
        (1 / Real.sqrt 2 : ℝ) ^ (l % 2) := by conv_lhs => rw [hdecomp, pow_add, pow_mul]
    _ ≤ (1 / 2 : ℝ) ^ (l / 2) := by rw [hp]; nlinarith [pow_nonneg (by norm_num : (0 : ℝ) ≤ 1 / 2) (l / 2)]
    _ = 1 / (2 ^ (l / 2) : ℝ) := by rw [div_pow]; simp

theorem log_le_bitBound (n : ℕ) (hn : 0 < n) :
    Real.logb 2 n ≤ (Nat.log 2 n + 1 : ℕ) := by
  apply le_of_lt
  rw [Real.logb_lt_iff_lt_rpow (by norm_num) (by exact_mod_cast hn)]
  rw [Real.rpow_natCast]
  exact_mod_cast Nat.lt_pow_succ_log_self (by norm_num : 1 < 2) n

theorem need_le_padExponent (P lam B C r need : ℕ) (hr : 0 < r)
    (hn : 0 < P * (lam + B)) (e : ℝ)
    (he : e ≤ 1 / (r : ℝ))
    (hneed : (need : ℝ) ≤ e * (P : ℝ) ^ 2 * lam +
      C * (P : ℝ) ^ 2 * Real.logb 2 (P * (lam + B))) :
    need ≤ padExponent P lam B C r := by
  have hquot : ((P ^ 2 * lam : ℕ) : ℝ) / r <
      ((P ^ 2 * lam / r + 1 : ℕ) : ℝ) := by
    apply (div_lt_iff₀ (by exact_mod_cast hr : (0 : ℝ) < r)).2
    have hm := Nat.mod_lt (P ^ 2 * lam) hr
    have hd := Nat.mod_add_div (P ^ 2 * lam) r
    have h : P ^ 2 * lam < (P ^ 2 * lam / r + 1) * r := by nlinarith
    exact_mod_cast h
  have hlog := log_le_bitBound (P * (lam + B)) hn
  have hterm : e * (P : ℝ) ^ 2 * lam ≤ ((P ^ 2 * lam : ℕ) : ℝ) / r := by
    have h := mul_le_mul_of_nonneg_right he
      (mul_nonneg (sq_nonneg (P : ℝ)) (Nat.cast_nonneg lam))
    calc
      e * (P : ℝ) ^ 2 * lam = e * ((P : ℝ) ^ 2 * lam) := by ring
      _ ≤ 1 / (r : ℝ) * ((P : ℝ) ^ 2 * lam) := h
      _ = ((P ^ 2 * lam : ℕ) : ℝ) / r := by push_cast; ring
  have hterm2 : (C : ℝ) * (P : ℝ) ^ 2 * Real.logb 2 (P * (lam + B)) ≤
      (C : ℝ) * (P : ℝ) ^ 2 * (Nat.log 2 (P * (lam + B)) + 1 : ℕ) := by
    have hlog' : Real.logb 2 ((P : ℝ) * (lam + B)) ≤
        (Nat.log 2 (P * (lam + B)) + 1 : ℕ) := by simpa only [Nat.cast_mul, Nat.cast_add] using hlog
    exact mul_le_mul_of_nonneg_left hlog' (by positivity)
  apply_mod_cast (show (need : ℝ) ≤ (padExponent P lam B C r : ℝ) by
    unfold padExponent
    push_cast at *
    linarith)

/-- The rate and pad bounds suffice to replace an existential required pad
with the explicit integer exponent. This does not establish the constructor
or the running time of a composed factoring machine. -/
theorem half_need_le_padExponent (P lam B C l need : ℕ)
    (hn : 0 < P * (lam + B)) (e : ℝ)
    (he : e ≤ (0.85 : ℝ) * (1 / Real.sqrt 2) ^ l)
    (hneed : (need : ℝ) ≤ e * (P : ℝ) ^ 2 * lam +
      C * (P : ℝ) ^ 2 * Real.logb 2 (P * (lam + B))) :
    need ≤ padExponent P lam B C (2 ^ (l / 2)) := by
  apply need_le_padExponent P lam B C _ need (by positivity) hn e _ hneed
  have hp : 0 ≤ (1 / Real.sqrt 2 : ℝ) ^ l := by positivity
  have h := half_power_le l
  push_cast
  linarith

end IntegerFactoring.HalfPad

#print axioms IntegerFactoring.HalfPad.half_power_le
#print axioms IntegerFactoring.HalfPad.log_le_bitBound
#print axioms IntegerFactoring.HalfPad.need_le_padExponent
#print axioms IntegerFactoring.HalfPad.half_need_le_padExponent
