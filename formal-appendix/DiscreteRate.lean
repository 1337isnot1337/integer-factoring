/-
SPDX-License-Identifier: Apache-2.0

Convert the selected-level padding rate to the challenge's discrete logarithmic scale.
-/
import MasterPad

open Groundwork PaddedDescent Levels.Gen
open IntegerFactoring.HalfPad

set_option maxRecDepth 4000

namespace IntegerFactoring.HalfAppendix

def rateDenom (n : ℕ) : ℕ := Nat.sqrt (Nat.log 2 n) + 1
def rateQuot (n : ℕ) : ℕ := n / rateDenom n + 1
def machineScale (n : ℕ) : ℕ := rateQuot n + Nat.log 2 n + 1

theorem square_le_four_pow (m : ℕ) : (m + 1) ^ 2 ≤ 4 * 2 ^ m := by
  induction m with
  | zero => norm_num
  | succ m ih =>
    by_cases hm : m ≤ 1
    · interval_cases m <;> norm_num
    · have hm' : 2 ≤ m := by omega
      have h : (m + 2) ^ 2 ≤ 2 * (m + 1) ^ 2 := by nlinarith
      rw [show (2 : ℕ) ^ (m + 1) = 2 ^ m * 2 by rw [pow_succ]]
      nlinarith

theorem log_le_rateQuot (n : ℕ) : Nat.log 2 n ≤ 4 * rateQuot n := by
  by_cases hn : n = 0
  · subst n; simp [rateQuot, rateDenom]
  let m := Nat.log 2 n
  let d := rateDenom n
  have hpow : 2 ^ m ≤ n := by
    have h := Nat.log2_self_le hn
    simpa [m, Nat.log2_eq_log_two] using h
  have hsqrt : Nat.sqrt m ≤ m := Nat.sqrt_le_self m
  have hsquare := square_le_four_pow m
  have hmd : m * d ≤ 4 * n := by
    dsimp [d, rateDenom]
    nlinarith
  have hd : 0 < d := by simp [d, rateDenom]
  have hdiv := Nat.div_add_mod n d
  have hmod := Nat.mod_lt n hd
  have hmul : d * m ≤ d * (4 * (n / d + 1)) := by nlinarith
  have h := Nat.le_of_mul_le_mul_left hmul hd
  simpa [m, d, rateQuot] using h

theorem machineScale_pos (n : ℕ) : 1 ≤ machineScale n := by simp [machineScale, rateQuot]

theorem machineScale_large (n : ℕ) : n ≤ 2 ^ machineScale n := by
  have h := Nat.lt_pow_succ_log_self Nat.one_lt_two n
  have hm : Nat.log 2 n + 1 ≤ machineScale n := by simp [machineScale, rateQuot]
  exact (Nat.le_of_lt h).trans (Nat.pow_le_pow_right (by norm_num) hm)

theorem machineScale_le (n : ℕ) : machineScale n ≤ 6 * rateQuot n := by
  have h := log_le_rateQuot n
  have hq : 1 ≤ rateQuot n := by simp [rateQuot]
  simp only [machineScale]
  omega

theorem rateDenom_le_halfLevel (ρ : ℝ) (x n : ℕ) (hx : 1 ≤ x)
    (ha : AppRunningTime.PartA ρ x) (hn : 0 < n) (hlarge : (K x 0) ^ 40 ≤ n) :
    rateDenom n ≤ (160 * (3 * x + 6) + 1) * 2 ^ (selectedLevel x n / 2) := by
  let l := selectedLevel x n
  let r := 2 ^ (l / 2)
  let A := 160 * (3 * x + 6)
  let m := Nat.log 2 n
  have hs := selectedLevel_spec x n hx hlarge
  have hlog : (m : ℝ) ≤ (A : ℝ) * (r : ℝ) ^ 2 := by
    have h := (Real.natLog_le_logb n 2).trans (log_le_half_level ρ x l n ha hn (le_of_lt hs.2))
    simpa [m, A, r, Nat.cast_add, Nat.cast_mul] using h
  have hlogN : m ≤ A * r ^ 2 := by exact_mod_cast hlog
  have hA : 1 ≤ A := by dsimp [A]; omega
  have hr : 1 ≤ r := Nat.one_le_iff_ne_zero.mpr (by dsimp [r]; exact pow_ne_zero _ (by norm_num))
  have hAA : A ≤ A * A := Nat.le_mul_of_pos_right A (by omega)
  have hbound : m ≤ (A * r) ^ 2 := by
    calc m ≤ A * r ^ 2 := hlogN
      _ ≤ (A * A) * r ^ 2 := Nat.mul_le_mul_right _ hAA
      _ = (A * r) ^ 2 := by ring
  have hsqrt : Nat.sqrt m ≤ A * r := by
    apply Nat.le_of_lt_succ
    apply Nat.sqrt_lt'.2
    nlinarith
  change Nat.sqrt m + 1 ≤ (A + 1) * r
  nlinarith

theorem masterExponent_le_rateQuot (ρ : ℝ) (x C : ℕ) (hx : 1 ≤ x)
    (ha : AppRunningTime.PartA ρ x) : ∃ A, ∀ n, masterExponent x C n ≤ A * rateQuot n := by
  let D := 160 * (3 * x + 6) + 1
  let E := 20004 + C * 500
  let N0 := (K x 0) ^ 40
  let small := (Finset.range N0).sup (masterExponent x C)
  refine ⟨max (D * E) small, ?_⟩
  intro n
  by_cases hlarge : N0 ≤ n
  · have hn : 0 < n := by
      have hK : 2 ≤ K x 0 := K_two x 0 hx
      have hp : 0 < N0 := by dsimp [N0]; positivity
      omega
    let r := 2 ^ (selectedLevel x n / 2)
    let d := rateDenom n
    have hd : 0 < d := by simp [d, rateDenom]
    have hr : 1 ≤ r := Nat.one_le_iff_ne_zero.mpr (by dsimp [r]; exact pow_ne_zero _ (by norm_num))
    have hden : d ≤ D * r := rateDenom_le_halfLevel ρ x n hx ha hn hlarge
    have hcost : r * masterExponent x C n ≤ E * n := masterExponent_cost ρ x C n hx ha hlarge
    have hmul : d * masterExponent x C n ≤ (D * E) * n := by
      calc d * masterExponent x C n ≤ (D * r) * masterExponent x C n := Nat.mul_le_mul_right _ hden
        _ = D * (r * masterExponent x C n) := by ring
        _ ≤ D * (E * n) := Nat.mul_le_mul_left _ hcost
        _ = (D * E) * n := by ring
    have hdiv := Nat.div_add_mod n d
    have hmod := Nat.mod_lt n hd
    have hn_le : n ≤ d * (n / d + 1) := by
      rw [Nat.mul_add, Nat.mul_one]
      omega
    have hcancel : d * masterExponent x C n ≤ d * ((D * E) * (n / d + 1)) := by
      calc d * masterExponent x C n ≤ (D * E) * n := hmul
        _ ≤ (D * E) * (d * (n / d + 1)) := Nat.mul_le_mul_left _ hn_le
        _ = d * ((D * E) * (n / d + 1)) := by ring
    have hq := Nat.le_of_mul_le_mul_left hcancel hd
    have hDE : D * E ≤ max (D * E) small := Nat.le_max_left _ _
    exact hq.trans (Nat.mul_le_mul_right _ hDE)
  · have hn : n ∈ Finset.range N0 := Finset.mem_range.mpr (Nat.lt_of_not_ge hlarge)
    have hsmall : masterExponent x C n ≤ small := Finset.le_sup hn
    have hsmall' : small ≤ max (D * E) small := Nat.le_max_right _ _
    have hq : 1 ≤ rateQuot n := by simp [rateQuot]
    exact hsmall.trans (hsmall'.trans (Nat.le_mul_of_pos_right _ hq))

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.masterExponent_le_rateQuot
#print axioms IntegerFactoring.HalfAppendix.log_le_rateQuot
#print axioms IntegerFactoring.HalfAppendix.machineScale_large
