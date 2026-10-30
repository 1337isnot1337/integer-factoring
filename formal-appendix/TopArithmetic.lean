/-
SPDX-License-Identifier: Apache-2.0

Integer width estimates for the selected factoring top table. This is only
the table-width component of the half-rate analysis, not a runtime proof.
-/
import TopGrouping
import FactoringAux.EngineCostRate

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem topWidth_bound (s K : ℕ) (hK : 1 ≤ K) :
    K ^ 2 * topWidth s K ≤ 720 * (2 ^ s) ^ 4 +
      306 * ((3 * K ^ 8) * (2 ^ s) ^ 3) := by
  simpa [topWidth, Ltop] using Ntop_S K (3 * K ^ 8) s hK

theorem ninth_power_le_S (K n : ℕ) (hK : 2 ≤ K) (hn : K ^ 40 ≤ n) :
    K ^ 9 ≤ 2 ^ sOf n := by
  have hK4 : 4 ≤ K ^ 4 := by
    calc
      4 ≤ 2 ^ 4 := by norm_num
      _ ≤ K ^ 4 := Nat.pow_le_pow_left hK 4
  have hf : 4 * (K ^ 9) ^ 4 ≤ n := by
    calc
      4 * (K ^ 9) ^ 4 ≤ K ^ 4 * (K ^ 9) ^ 4 := Nat.mul_le_mul_right _ hK4
      _ = K ^ 40 := by ring
      _ ≤ n := hn
  have hnS := n_le_S n
  have hS : (K ^ 9) ^ 4 ≤ (2 ^ sOf n) ^ 4 := by omega
  by_contra h
  have hlt : 2 ^ sOf n < K ^ 9 := by omega
  have hp := Nat.pow_lt_pow_left hlt (by norm_num : (4 : ℕ) ≠ 0)
  omega

/-- A conservative uniform bound for the selected table width. -/
theorem selected_width_le (K n : ℕ) (hK : 2 ≤ K) (hn : K ^ 40 ≤ n) :
    K ^ 2 * (topWidth (sOf n) K + 1) ≤ 6000 * n := by
  let S := 2 ^ sOf n
  have hKS : K ^ 9 ≤ S := by
    exact ninth_power_le_S K n hK hn
  have hprod : K * ((3 * K ^ 8) * S ^ 3) ≤ 3 * S ^ 4 := by
    calc
      K * ((3 * K ^ 8) * S ^ 3) = K ^ 9 * (3 * S ^ 3) := by ring
      _ ≤ S * (3 * S ^ 3) := Nat.mul_le_mul_right _ hKS
      _ = 3 * S ^ 4 := by ring
  have hprod2 : 2 * ((3 * K ^ 8) * S ^ 3) ≤ 3 * S ^ 4 :=
    (Nat.mul_le_mul_right _ hK).trans hprod
  have hw := topWidth_bound (sOf n) K (by omega)
  have hw' : 2 * (K ^ 2 * topWidth (sOf n) K) ≤ 2358 * S ^ 4 := by
    change K ^ 2 * topWidth (sOf n) K ≤ 720 * S ^ 4 + 306 * ((3 * K ^ 8) * S ^ 3) at hw
    omega
  have hSn : S ^ 4 ≤ 4 * n + 4 := by simpa [S] using S_le_n n
  have hK2 : K ^ 2 ≤ n := (Nat.pow_le_pow_right (by omega : 0 < K) (by omega : 2 ≤ 40)).trans hn
  have hn1 : 1 ≤ n := by
    have hp : 0 < K ^ 40 := by positivity
    omega
  nlinarith

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.topWidth_bound
#print axioms IntegerFactoring.HalfAppendix.selected_width_le
