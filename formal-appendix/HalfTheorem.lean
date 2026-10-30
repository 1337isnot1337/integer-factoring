/-
SPDX-License-Identifier: Apache-2.0

The reference factoring-machine square-root-log bound.
-/
import SelectedFactor
import DiscreteRate
import PaddedMachineRate

open Groundwork Groundwork.Prog PaddedDescent PaddedDescent.F Levels.Gen

namespace IntegerFactoring.HalfAppendix

private theorem size_bools_half (w : List Bool) : size w = 2 * w.length + 1 := by
  induction w with
  | nil => rfl
  | cons b w ih => simp only [size_list, List.map_cons, List.sum_cons, size_bool, List.length_cons] at ih ⊢; omega

private theorem factor_output_size (N : ℕ) :
    size (Computability.encodeNat (Factoring.factor N)) ≤ 2 * (Computability.encodeNat N).length + 5 := by
  rw [size_bools_half, encodeNat_length, encodeNat_length]
  by_cases h : N = 0
  · subst N
    change 2 * (Nat.minFac 0).size + 1 ≤ 2 * (Nat.size 0) + 5
    decide
  · have hpos : 0 < N := Nat.pos_of_ne_zero h
    have hsize : (Factoring.factor N).size ≤ N.size := Nat.size_le_size (Nat.minFac_le hpos)
    omega

private theorem output_within_exp (n A : ℕ) (hA : 4 ≤ A) :
    2 * n + 5 ≤ 2 ^ (A * machineScale n) := by
  let t := machineScale n
  let H := 2 ^ t
  have ht : 1 ≤ t := machineScale_pos n
  have hn : n ≤ H := machineScale_large n
  have hH : 2 ≤ H := by
    dsimp [H]
    simpa using Nat.pow_le_pow_right (by norm_num : 0 < (2 : ℕ)) ht
  have h1 : 2 * n + 5 ≤ 8 * H := by omega
  have h2 : 8 * H ≤ 2 ^ (4 * t) := by
    calc 8 * H = 2 ^ (t + 3) := by dsimp [H]; rw [pow_add]; norm_num; ring
      _ ≤ 2 ^ (4 * t) := Nat.pow_le_pow_right (by norm_num) (by omega)
  exact h1.trans (h2.trans (Nat.pow_le_pow_right (by norm_num) (Nat.mul_le_mul_right _ hA)))

theorem machineScale_time_to_challenge {D : ℕ} (M : Factoring.FactoringMachine)
    (hfin : M.FiniteAlphabets) (hM : ∀ n, M.time n ≤ 2 ^ (D * machineScale n)) :
    Factoring.SqrtLogLogConjecture := by
  let C := max (2 ^ (6 * D)) (6 * D)
  have hC1 : 2 ^ (6 * D) ≤ C := Nat.le_max_left _ _
  have hC2 : 6 * D ≤ C := Nat.le_max_right _ _
  refine ⟨C, M, hfin, ?_⟩
  intro n
  let d := rateDenom n
  let q := n / d
  have hd : 0 < d := by simp [d, rateDenom]
  have hscale : D * machineScale n ≤ (6 * D) * (q + 1) := by
    have h := Nat.mul_le_mul_left D (machineScale_le n)
    dsimp [q, d, rateQuot] at *
    nlinarith
  have hq : q * d ≤ n := by
    have h := Nat.div_mul_le_self n d
    simpa [q, Nat.mul_comm] using h
  have hdiv : (6 * D) * q ≤ C * n / d := by
    apply (Nat.le_div_iff_mul_le hd).2
    calc ((6 * D) * q) * d = (6 * D) * (q * d) := by ring
      _ ≤ C * n := Nat.mul_le_mul hC2 hq
  calc M.time n ≤ 2 ^ (D * machineScale n) := hM n
    _ ≤ 2 ^ ((6 * D) * (q + 1)) := Nat.pow_le_pow_right (by norm_num) hscale
    _ = 2 ^ (6 * D) * 2 ^ ((6 * D) * q) := by
      rw [Nat.mul_add, Nat.mul_one, pow_add, Nat.mul_comm]
    _ ≤ C * 2 ^ (C * n / d) := Nat.mul_le_mul hC1 (Nat.pow_le_pow_right (by norm_num) hdiv)

theorem sqrtLogLog_factoring : Factoring.SqrtLogLogConjecture := by
  obtain ⟨x, C, run, hx, hrun, hcorrect⟩ := selectedFactor_correct
  have hx' : 1 ≤ x := by omega
  have ha : AppRunningTime.PartA (3 / 4 : ℝ) x :=
    AppRunningTime.thmRunningTime.a (3 / 4 : ℝ) x (by norm_num) (by norm_num) hx
  obtain ⟨A0, hA0⟩ := masterExponent_le_rateQuot (3 / 4 : ℝ) x C hx' ha
  let A := max A0 4
  have hA0' : A0 ≤ A := Nat.le_max_left _ _
  have hA : 4 ≤ A := Nat.le_max_right _ _
  have he := selectedExponent_polyT x C hx'
  have hF := selectedFactor_polyT x C hx' run hrun
  have hE : ∀ N, selectedExponent x C (Computability.encodeNat N) ≤ A * machineScale (Computability.encodeNat N).length := by
    intro N
    let n := (Computability.encodeNat N).length
    have h := hA0 n
    change masterExponent x C n ≤ A * machineScale n
    calc masterExponent x C n ≤ A0 * rateQuot n := h
      _ ≤ A * rateQuot n := Nat.mul_le_mul_right _ hA0'
      _ ≤ A * machineScale n := Nat.mul_le_mul_left _ (by simp only [machineScale]; omega)
  have hc : ∀ N, selectedFactor x C run (Computability.encodeNat N, pad (selectedExponent x C (Computability.encodeNat N))) =
      Computability.encodeNat (Factoring.factor N) := by
    intro N
    apply hcorrect N
    simp [pad]
  have hW : ∀ N, size (selectedFactor x C run (Computability.encodeNat N, pad (selectedExponent x C (Computability.encodeNat N)))) ≤
      2 ^ (A * machineScale (Computability.encodeNat N).length) := by
    intro N
    rw [hc N]
    exact (factor_output_size N).trans (output_within_exp _ A hA)
  have hn : ∀ n, n ≤ 2 ^ (A * machineScale n) := by
    intro n
    exact (machineScale_large n).trans (Nat.pow_le_pow_right (by norm_num) (by nlinarith [machineScale_pos n]))
  obtain ⟨D, M, hfin, htime⟩ := padded_factoring_exp A machineScale he hF machineScale_pos hn hE hW hc
  exact machineScale_time_to_challenge M hfin htime

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.sqrtLogLog_factoring
