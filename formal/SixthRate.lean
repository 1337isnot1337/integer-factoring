/-
SPDX-License-Identifier: Apache-2.0

Arithmetic rate bridge for the sixth-power candidate. The schedule is adapted
from the pinned reference's Kb / runtime depth construction, with the proposed
recurrence K_(d+1) = 25 K_d^12. Reference authors: Isaac Rudich and Louis-Martin
Rousseau; reference commit c91cc3f4947caca8977e81a5d14244f73b5ae2b1.

These theorems do not assume a cost bound for an unimplemented evaluator.
They state precisely the cutoff, volume, and padding hypotheses needed from
that evaluator. Prepared with OpenAI Codex assistance.
-/
module

public import Mathlib.Data.Nat.Log
public import Mathlib.Analysis.SpecialFunctions.Pow.Real
public import Mathlib.Tactic.NormNum
public import Mathlib.Tactic.Ring
public import Mathlib.Tactic.Positivity
public import Mathlib.Tactic.Linarith

@[expose] public section

namespace Factoring.Improvements.SixthRate

/-- The proposed group-count schedule, independently of any evaluator. -/
def groups (K0 : ℕ) : ℕ → ℕ
  | 0 => K0
  | d + 1 => 25 * groups K0 d ^ 12

/-- Connect any implementation of this recurrence to the rate schedule. -/
theorem groups_eq_of_recurrence (K0 : ℕ) (K : ℕ → ℕ) (h0 : K 0 = K0)
    (hstep : ∀ d, K (d + 1) = 25 * K d ^ 12) : ∀ d, K d = groups K0 d
  | 0 => h0
  | d + 1 => by rw [hstep, groups, groups_eq_of_recurrence K0 K h0 hstep d]

/-- A deliberately loose bound with an explicit absolute constant. -/
theorem groups_scaled_le {K0 : ℕ} (h0 : K0 ≤ 1024) (d : ℕ) :
    2 * groups K0 d ≤ 2 ^ (11 * 12 ^ d) := by
  induction d with
  | zero => norm_num [groups]; omega
  | succ d ih =>
    calc
      2 * groups K0 (d + 1) = 50 * groups K0 d ^ 12 := by simp [groups]; ring
      _ ≤ (2 * groups K0 d) ^ 12 := by
        rw [mul_pow]
        exact Nat.mul_le_mul_right _ (by norm_num : 50 ≤ 2 ^ 12)
      _ ≤ (2 ^ (11 * 12 ^ d)) ^ 12 := Nat.pow_le_pow_left ih 12
      _ = 2 ^ (11 * 12 ^ (d + 1)) := by rw [← pow_mul, pow_succ]; congr 1; ring

/-- A next-level cutoff supplies the lower bound on selected depth that is
essential for a quantitative saving. A bound on the current group count alone
does not supply this hypothesis. -/
theorem log_le_of_cutoff {K0 n d : ℕ} (h0 : K0 ≤ 1024)
    (hn : n < groups K0 (d + 1) ^ 40) :
    Nat.log 2 n + 1 ≤ 5280 * 12 ^ d := by
  have hK : groups K0 (d + 1) ≤ 2 ^ (11 * 12 ^ (d + 1)) :=
    le_trans (by omega) (groups_scaled_le h0 (d + 1))
  have hp : groups K0 (d + 1) ^ 40 ≤ 2 ^ (5280 * 12 ^ d) := by
    calc
      groups K0 (d + 1) ^ 40 ≤ (2 ^ (11 * 12 ^ (d + 1))) ^ 40 :=
        Nat.pow_le_pow_left hK 40
      _ = 2 ^ (5280 * 12 ^ d) := by rw [← pow_mul, pow_succ]; congr 1; ring
  have hlog := Nat.log_lt_of_lt_pow' (by positivity : 5280 * 12 ^ d ≠ 0) (hn.trans_le hp)
  omega

/-- The strict integer certificate giving a saving exponent 1/30. -/
theorem contraction_certificate : 12 * 23 ^ 30 < (25 : ℕ) ^ 30 := by norm_num

/-- The audited reference exponent log(12/11)/log(16) is strictly below
1/30; this comparison is expressed without rounded real logarithms. -/
theorem reference_comparison : (12 : ℕ) ^ 30 < 16 * 11 ^ 30 := by norm_num

/-- Core rate bridge. It uses only the actual depth and a cross-multiplied
volume estimate; there are no divisions, real approximations, or asymptotic
constants hidden in its conclusion. -/
theorem volume_rate {n d v A : ℕ}
    (hlog : Nat.log 2 n + 1 ≤ 5280 * 12 ^ d)
    (hv : 25 ^ d * v ≤ A * 23 ^ d * n) :
    v ^ 30 * (Nat.log 2 n + 1) ≤ 5280 * A ^ 30 * n ^ 30 := by
  have hp := Nat.pow_le_pow_left hv 30
  have hc := Nat.pow_le_pow_left contraction_certificate.le d
  have hscale : (25 ^ d) ^ 30 * (v ^ 30 * (Nat.log 2 n + 1)) ≤
      (25 ^ d) ^ 30 * (5280 * A ^ 30 * n ^ 30) := by
    calc
      (25 ^ d) ^ 30 * (v ^ 30 * (Nat.log 2 n + 1)) =
          (25 ^ d * v) ^ 30 * (Nat.log 2 n + 1) := by rw [mul_pow]; ring
      _ ≤ (A * 23 ^ d * n) ^ 30 * (5280 * 12 ^ d) := Nat.mul_le_mul hp hlog
      _ = 5280 * A ^ 30 * n ^ 30 * (12 * 23 ^ 30) ^ d := by
        simp only [mul_pow, ← pow_mul]; ring
      _ ≤ 5280 * A ^ 30 * n ^ 30 * (25 ^ 30) ^ d := Nat.mul_le_mul_left _ hc
      _ = (25 ^ d) ^ 30 * (5280 * A ^ 30 * n ^ 30) := by
        simp only [← pow_mul]; ring
  exact Nat.le_of_mul_le_mul_left hscale (by positivity)

/-- Apply the rate bridge at any depth satisfying the actual next-level
cutoff. Small inputs, including zero, are covered by the same statement. -/
theorem volume_rate_of_cutoff {K0 n d v A : ℕ} (h0 : K0 ≤ 1024)
    (hn : n < groups K0 (d + 1) ^ 40)
    (hv : 25 ^ d * v ≤ A * 23 ^ d * n) :
    v ^ 30 * (Nat.log 2 n + 1) ≤ 5280 * A ^ 30 * n ^ 30 :=
  volume_rate (log_le_of_cutoff h0 hn) hv

/-- Adding a separately bounded error term preserves the power bound. In
particular, this theorem requires an explicit budget for all padding and
bookkeeping; it does not silently discard those terms. -/
theorem padded_rate {n m v p E B C : ℕ}
    (hE : E ≤ v + p)
    (hv : v ^ 30 * (Nat.log 2 n + 1) ≤ B * m ^ 30)
    (hp : p ^ 30 * (Nat.log 2 n + 1) ≤ C * m ^ 30) :
    E ^ 30 * (Nat.log 2 n + 1) ≤ 2 ^ 29 * (B + C) * m ^ 30 := by
  calc
    E ^ 30 * (Nat.log 2 n + 1) ≤ (v + p) ^ 30 * (Nat.log 2 n + 1) :=
      Nat.mul_le_mul_right _ (Nat.pow_le_pow_left hE 30)
    _ ≤ (2 ^ 29 * (v ^ 30 + p ^ 30)) * (Nat.log 2 n + 1) :=
      Nat.mul_le_mul_right _ (add_pow_le (Nat.zero_le v) (Nat.zero_le p) 30)
    _ = 2 ^ 29 * (v ^ 30 * (Nat.log 2 n + 1) + p ^ 30 * (Nat.log 2 n + 1)) := by ring
    _ ≤ 2 ^ 29 * (B * m ^ 30 + C * m ^ 30) :=
      Nat.mul_le_mul_left _ (Nat.add_le_add hv hp)
    _ = 2 ^ 29 * (B + C) * m ^ 30 := by ring

/-- The explicit real-valued logarithmic saving follows from the integral
power bound. The numerator m can be n+1 to cover positive finite-input costs;
the logarithm plus one keeps the denominator positive even at zero. -/
theorem real_rate_of_power {n m v B : ℕ}
    (hv : v ^ 30 * (Nat.log 2 n + 1) ≤ B * m ^ 30) :
    (v : ℝ) ≤ (B : ℝ) ^ ((30 : ℝ)⁻¹) * m /
      ((Nat.log 2 n + 1 : ℕ) : ℝ) ^ ((30 : ℝ)⁻¹) := by
  let L : ℝ := (Nat.log 2 n + 1 : ℕ)
  have hL : 0 < L := by dsimp [L]; positivity
  have hroot : 0 < L ^ ((30 : ℝ)⁻¹) := Real.rpow_pos_of_pos hL _
  have hcast : (v : ℝ) ^ 30 * L ≤ (B : ℝ) * (m : ℝ) ^ 30 := by
    dsimp [L]
    exact_mod_cast hv
  have hrootpow := Real.rpow_inv_natCast_pow (n := 30) hL.le (by norm_num)
  simp only [Nat.cast_ofNat] at hrootpow
  have hpow : ((v : ℝ) * L ^ ((30 : ℝ)⁻¹)) ^ 30 ≤ (B : ℝ) * (m : ℝ) ^ 30 := by
    rw [mul_pow, hrootpow]
    exact hcast
  have hp30 (x : ℝ) : x ^ (30 : ℝ) = x ^ (30 : ℕ) := by
    simpa only [Nat.cast_ofNat] using Real.rpow_natCast x 30
  have h := (Real.le_rpow_inv_iff_of_pos
    (by positivity : 0 ≤ (v : ℝ) * L ^ ((30 : ℝ)⁻¹))
    (by positivity : 0 ≤ (B : ℝ) * (m : ℝ) ^ 30)
    (by norm_num : (0 : ℝ) < 30)).2 (by rw [hp30]; exact hpow)
  have hnroot := Real.pow_rpow_inv_natCast (n := 30) (Nat.cast_nonneg m) (by norm_num)
  simp only [Nat.cast_ofNat] at hnroot
  rw [Real.mul_rpow (by positivity) (by positivity), hnroot] at h
  exact (le_div_iff₀ hroot).2 h

theorem real_volume_rate_of_cutoff {K0 n d v A : ℕ} (h0 : K0 ≤ 1024)
    (hn : n < groups K0 (d + 1) ^ 40)
    (hv : 25 ^ d * v ≤ A * 23 ^ d * n) :
    (v : ℝ) ≤ (5280 * A ^ 30 : ℕ) ^ ((30 : ℝ)⁻¹) * n /
      ((Nat.log 2 n + 1 : ℕ) : ℝ) ^ ((30 : ℝ)⁻¹) :=
  real_rate_of_power (volume_rate_of_cutoff h0 hn hv)

/-- Carry an actual, separately established exponential time bound through
the rate conversion. The hypothesis must bound the complete machine time,
including its compiler and simulation overhead, when t is that machine time. -/
theorem time_rate_of_power {n E B C t : ℕ}
    (hE : E ^ 30 * (Nat.log 2 n + 1) ≤ B * (n + 1) ^ 30)
    (ht : t ≤ C * 2 ^ (C * (E + 1))) :
    (t : ℝ) ≤ (C : ℝ) * (2 : ℝ) ^ ((C : ℝ) *
      ((B : ℝ) ^ ((30 : ℝ)⁻¹) * (n + 1) /
        ((Nat.log 2 n + 1 : ℕ) : ℝ) ^ ((30 : ℝ)⁻¹) + 1)) := by
  have htc : (t : ℝ) ≤ (C : ℝ) * (2 : ℝ) ^ (C * (E + 1) : ℕ) := by exact_mod_cast ht
  rw [← Real.rpow_natCast] at htc
  push_cast at htc
  refine le_trans htc (mul_le_mul_of_nonneg_left ?_ (Nat.cast_nonneg C))
  apply Real.rpow_le_rpow_of_exponent_le (by norm_num : (1 : ℝ) ≤ 2)
  have hbound := real_rate_of_power hE
  push_cast at hbound ⊢
  exact mul_le_mul_of_nonneg_left (by linarith) (Nat.cast_nonneg C)

/-- The target exponent is strictly better than the audited reference rate;
this is an exact real-logarithm comparison, not a decimal comparison. -/
theorem exponent_comparison :
    Real.log (12 / 11) / Real.log 16 < (30 : ℝ)⁻¹ ∧
      (30 : ℝ)⁻¹ < Real.log (25 / 23) / Real.log 12 := by
  constructor
  · have hc : (12 : ℝ) ^ 30 < 16 * 11 ^ 30 := by exact_mod_cast reference_comparison
    have hr : (12 / 11 : ℝ) ^ 30 < 16 := by
      rw [div_pow]
      exact (div_lt_iff₀ (by positivity)).2 hc
    have hlog := Real.log_lt_log (by positivity : (0 : ℝ) < (12 / 11 : ℝ) ^ 30) hr
    rw [Real.log_pow] at hlog
    apply (div_lt_iff₀ (Real.log_pos (by norm_num : (1 : ℝ) < 16))).2
    norm_num at hlog ⊢
    linarith
  · have hc : (12 : ℝ) * 23 ^ 30 < 25 ^ 30 := by exact_mod_cast contraction_certificate
    have hr : (12 : ℝ) < (25 / 23 : ℝ) ^ 30 := by
      rw [div_pow]
      exact (lt_div_iff₀ (by positivity)).2 hc
    have hlog := Real.log_lt_log (by norm_num : (0 : ℝ) < 12) hr
    rw [Real.log_pow] at hlog
    apply (lt_div_iff₀ (Real.log_pos (by norm_num : (1 : ℝ) < 12))).2
    norm_num at hlog ⊢
    linarith

end Factoring.Improvements.SixthRate

#print axioms Factoring.Improvements.SixthRate.groups_scaled_le
#print axioms Factoring.Improvements.SixthRate.groups_eq_of_recurrence
#print axioms Factoring.Improvements.SixthRate.log_le_of_cutoff
#print axioms Factoring.Improvements.SixthRate.contraction_certificate
#print axioms Factoring.Improvements.SixthRate.reference_comparison
#print axioms Factoring.Improvements.SixthRate.volume_rate
#print axioms Factoring.Improvements.SixthRate.volume_rate_of_cutoff
#print axioms Factoring.Improvements.SixthRate.padded_rate
#print axioms Factoring.Improvements.SixthRate.real_rate_of_power
#print axioms Factoring.Improvements.SixthRate.real_volume_rate_of_cutoff
#print axioms Factoring.Improvements.SixthRate.time_rate_of_power
#print axioms Factoring.Improvements.SixthRate.exponent_comparison
