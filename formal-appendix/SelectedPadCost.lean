/-
SPDX-License-Identifier: Apache-2.0

Numerical pad bounds for the selected top parameters. These arithmetic
statements do not certify the machine that writes the pads or tables.
-/
import BoundedTop
import TopStateArithmetic
import PadCost
import HalfSchedule

open Levels.Gen Factoring.Aux IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem selected_written_exponents (P n r : ℕ) (hP : 2 ≤ P)
    (hn : P ^ 40 ≤ n) (hrP : r ≤ P) :
    r * ((P - 1) * (topWidth (sOf n) P + 1)) ≤ 6000 * n ∧
      r * topEntryBits (sOf n) P ≤ 8000 * n := by
  constructor
  · calc
      r * ((P - 1) * (topWidth (sOf n) P + 1)) ≤
          P * (P * (topWidth (sOf n) P + 1)) := by gcongr <;> omega
      _ = P ^ 2 * (topWidth (sOf n) P + 1) := by ring
      _ ≤ 6000 * n := selected_width_le P n hP hn
  · exact (Nat.mul_le_mul_right _ hrP).trans (selected_state_exponent_le P n hP hn)

theorem selected_log_overhead (P n : ℕ) (hP : 2 ≤ P)
    (hnlo : P ^ 40 ≤ n) (hnhi : n ≤ P ^ 400) :
    Nat.log 2 (P * ((topWidth (sOf n) P + 1) + topEntryBits (sOf n) P)) + 1 ≤ 500 * P := by
  let L := topWidth (sOf n) P + 1
  let B := topEntryBits (sOf n) P
  have hw : P ^ 2 * L ≤ 6000 * n := selected_width_le P n hP hnlo
  have hb : P * B ≤ 8000 * n := selected_state_exponent_le P n hP hnlo
  have hPL : P * L ≤ P ^ 2 * L :=
    Nat.mul_le_mul_right L (Nat.le_self_pow (by norm_num : 2 ≠ 0) P)
  have hQ : P * (L + B) ≤ 14000 * n := by nlinarith
  have hPexp : P ^ 400 ≤ 2 ^ (400 * P) := by
    calc
      P ^ 400 ≤ (2 ^ P) ^ 400 := Nat.pow_le_pow_left (Nat.lt_two_pow_self.le) 400
      _ = 2 ^ (400 * P) := by rw [← pow_mul]; congr 1; ring
  have hbig : 14000 * n ≤ 2 ^ (400 * P + 14) := by
    calc
      14000 * n ≤ 14000 * P ^ 400 := Nat.mul_le_mul_left _ hnhi
      _ ≤ 14000 * 2 ^ (400 * P) := Nat.mul_le_mul_left _ hPexp
      _ ≤ 2 ^ 14 * 2 ^ (400 * P) := Nat.mul_le_mul_right _ (by norm_num)
      _ = 2 ^ (400 * P + 14) := by rw [pow_add]; ring
  have hlt : P * (L + B) < 2 ^ (500 * P) :=
    hQ.trans_lt (hbig.trans_lt (Nat.pow_lt_pow_right (by norm_num) (by omega)))
  have hpos : P * (L + B) ≠ 0 := by
    have hL : 0 < L := by dsimp [L]; omega
    exact Nat.ne_of_gt (Nat.mul_pos (by omega) (by omega))
  have hlog := (Nat.log_lt_iff_lt_pow (by norm_num : 1 < 2) hpos).2 hlt
  change Nat.log 2 (P * (L + B)) + 1 ≤ 500 * P
  omega

theorem selected_pad_cost (P n C r : ℕ) (hP : 2 ≤ P)
    (hnlo : P ^ 40 ≤ n) (hnhi : n ≤ P ^ 400) (hrP : r ≤ P) :
    r * padExponent P (topWidth (sOf n) P + 1) (topEntryBits (sOf n) P) C r ≤
      (6001 + C * 500) * n := by
  have hP4 : P ^ 4 ≤ n :=
    (Nat.pow_le_pow_right (by omega : 0 < P) (by omega : 4 ≤ 40)).trans hnlo
  have h := pad_exponent_mul_le P (topWidth (sOf n) P + 1) (topEntryBits (sOf n) P)
    C r n 6000 500 hrP hP4 (selected_width_le P n hP hnlo) (selected_log_overhead P n hP hnlo hnhi)
  simpa using h

theorem scheduled_pad_cost (ρ : ℝ) (x l n C : ℕ)
    (hx : 1 ≤ x) (ha : AppRunningTime.PartA ρ x)
    (hnlo : (K x l) ^ 40 ≤ n) (hnhi : n ≤ (K x (l + 1)) ^ 40) :
    2 ^ (l / 2) * padExponent (K x l) (topWidth (sOf n) (K x l) + 1)
      (topEntryBits (sOf n) (K x l)) C (2 ^ (l / 2)) ≤ (6001 + C * 500) * n := by
  have hP := K_two x l hx
  have hhi : n ≤ (K x l) ^ 400 := by
    calc
      n ≤ (K x (l + 1)) ^ 40 := hnhi
      _ ≤ ((K x l) ^ 10) ^ 40 := Nat.pow_le_pow_left (ha.K_lt_le l).2 40
      _ = (K x l) ^ 400 := by rw [← pow_mul]
  exact selected_pad_cost (K x l) n C (2 ^ (l / 2)) hP hnlo hhi (half_level_le_K x l hx)

theorem scheduled_pad_sqrt_cost (ρ : ℝ) (x l n C : ℕ)
    (hx : 1 ≤ x) (ha : AppRunningTime.PartA ρ x) (hn : 1 < n)
    (hnlo : (K x l) ^ 40 ≤ n) (hnhi : n ≤ (K x (l + 1)) ^ 40) :
    (padExponent (K x l) (topWidth (sOf n) (K x l) + 1)
      (topEntryBits (sOf n) (K x l)) C (2 ^ (l / 2)) : ℝ) ≤
      ((6001 + C * 500 : ℕ) : ℝ) *
        Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * (n : ℝ) / Real.sqrt (Real.log n) := by
  let q := padExponent (K x l) (topWidth (sOf n) (K x l) + 1)
    (topEntryBits (sOf n) (K x l)) C (2 ^ (l / 2))
  let r := 2 ^ (l / 2)
  let D : ℝ := Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2)
  have hq : (r : ℝ) * q ≤ ((6001 + C * 500 : ℕ) : ℝ) * n := by
    exact_mod_cast scheduled_pad_cost ρ x l n C hx ha hnlo hnhi
  have hs : Real.sqrt (Real.log n) ≤ D * r := sqrt_log_le_half_level ρ x l n ha (by omega) hnhi
  have hD : 0 ≤ D := Real.sqrt_nonneg _
  have hden : 0 < Real.sqrt (Real.log n) := Real.sqrt_pos.2 (Real.log_pos (by exact_mod_cast hn))
  change (q : ℝ) ≤ ((6001 + C * 500 : ℕ) : ℝ) * D * n / Real.sqrt (Real.log n)
  apply (le_div_iff₀ hden).2
  calc
    (q : ℝ) * Real.sqrt (Real.log n) ≤ q * (D * r) := mul_le_mul_of_nonneg_left hs (by positivity)
    _ = D * ((r : ℝ) * q) := by ring
    _ ≤ D * (((6001 + C * 500 : ℕ) : ℝ) * n) := mul_le_mul_of_nonneg_left hq hD
    _ = ((6001 + C * 500 : ℕ) : ℝ) * D * n := by ring

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_log_overhead
#print axioms IntegerFactoring.HalfAppendix.selected_written_exponents
#print axioms IntegerFactoring.HalfAppendix.selected_pad_cost
#print axioms IntegerFactoring.HalfAppendix.scheduled_pad_cost
#print axioms IntegerFactoring.HalfAppendix.scheduled_pad_sqrt_cost
