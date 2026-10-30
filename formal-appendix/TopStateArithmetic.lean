/-
SPDX-License-Identifier: Apache-2.0

An integer bound for the exponent in the selected top state count. It does
not by itself bound table entries or prove an algorithmic running time.
-/
import TopArithmetic

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem selected_state_exponent_le (K n : ℕ) (hK : 2 ≤ K) (hn : K ^ 40 ≤ n) :
    K * (Wc (sOf n) (3 * K ^ 8) * (2 * (3 * K ^ 8) / K) + 1) ≤ 8000 * n := by
  let s := sOf n
  let S := 2 ^ s
  let g := 3 * K ^ 8
  have hKS : K ^ 9 ≤ S := ninth_power_le_S K n hK hn
  have hprod : K * (g * S ^ 3) ≤ 3 * S ^ 4 := by
    calc
      K * (g * S ^ 3) = K ^ 9 * (3 * S ^ 3) := by dsimp [g]; ring
      _ ≤ S * (3 * S ^ 3) := Nat.mul_le_mul_right _ hKS
      _ = 3 * S ^ 4 := by ring
  have hprod2 : 2 * (g * S ^ 3) ≤ 3 * S ^ 4 :=
    (Nat.mul_le_mul_right _ hK).trans hprod
  have hwc := Wc_mul_le s g
  change g * Wc s g ≤ 10 * (pA s * pA s * fw s) +
    g * (16 * (pA s * fw s) + 4 * fw s) at hwc
  have hA : pA s * pA s * fw s ≤ 36 * S ^ 4 := A2fw_le s
  have hAf : pA s * fw s ≤ 18 * S ^ 3 := Afw_le s
  have hf : fw s ≤ 9 * S ^ 3 := fw_le s
  have hsmall : 16 * (pA s * fw s) + 4 * fw s ≤ 324 * S ^ 3 := by omega
  have hmul := Nat.mul_le_mul_left g hsmall
  have hscaled : 2 * (g * Wc s g) ≤ 1692 * S ^ 4 := by nlinarith
  have hdiv : K * (2 * g / K) ≤ 2 * g := Nat.mul_div_le (2 * g) K
  have hgroup : K * (Wc s g * (2 * g / K)) ≤ 2 * (g * Wc s g) := by nlinarith
  have hSn : S ^ 4 ≤ 4 * n + 4 := S_le_n n
  have hKle : K ≤ n := (Nat.le_self_pow (by norm_num : 40 ≠ 0) K).trans hn
  have hn16 : 16 ≤ n := by
    calc
      16 = 2 ^ 4 := by norm_num
      _ ≤ K ^ 4 := Nat.pow_le_pow_left hK 4
      _ ≤ K ^ 40 := Nat.pow_le_pow_right (by omega : 0 < K) (by omega)
      _ ≤ n := hn
  change K * (Wc s g * (2 * g / K) + 1) ≤ 8000 * n
  nlinarith

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_state_exponent_le
