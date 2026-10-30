/-
SPDX-License-Identifier: Apache-2.0

The reference mixed-radix coding transported to the integer network. Its
polynomial-time decoder and row summand still need to be certified.
-/
import SelectedRowBridge
import FactoringAux.EngineMul

open Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def mulCodingZ (s N k y : ℕ) : Coding (mulNet s ℤ N k y) where
  s := mulSz s
  pos _ := mulSz_pos _ _
  e a := (gcode s ((MG.enum (pA s)).symm a)).trans (finCongr (mulSz_fin s a).symm)

def clusterCodingZ (s g N k y : ℕ) (hg : 0 < g) : Coding (clusterNet s ℤ N k y g hg) :=
  (mulCodingZ s N k y).coarsen (fun a => clMG s g hg ((MG.enum (pA s)).symm a))

def selectedCodingZ (s N k y P : ℕ) (hP : 2 ≤ P) :
    Coding ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen (selectedTopGrouping s N k y P hP)) :=
  (clusterCodingZ s (3 * P ^ 8) N k y (by positivity)).coarsen (selectedTopGrouping s N k y P hP)

theorem selected_explicit_code_bound (s N k y P : ℕ) (hP : 2 ≤ P) (γ : Fin P) :
    (selectedCodingZ s N k y P hP).s γ ≤ 2 ^ topEntryBits s P := by
  have h := selectedTop_state_card s N k y P hP γ
  have hc := Fintype.card_congr ((selectedCodingZ s N k y P hP).e γ)
  simp only [Fintype.card_fin] at hc
  have hp : Wc s (3 * P ^ 8) * (2 * (3 * P ^ 8) / P) ≤ topEntryBits s P := by simp [topEntryBits]
  exact hc ▸ h.trans (Nat.pow_le_pow_right (by norm_num) hp)

theorem selected_explicit_row_as_padded_sum (x C n N k y : ℕ)
    (hP : 2 ≤ K x (selectedLevel x n)) (pad : List Unit)
    (hpad : 2 ^ masterExponent x C n ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n)))
    (θ : Fin (K x (selectedLevel x n)) → Bits (topWidth (sOf n) (K x (selectedLevel x n)) + 1)) :
    let P := K x (selectedLevel x n)
    let M := (clusterNet (sOf n) ℤ N k y (3 * P ^ 8) (by positivity)).coarsen
      (selectedTopGrouping (sOf n) N k y P hP)
    let codec := (selectedCodingZ (sOf n) N k y P hP).e γ
    paddedRowSumZ (2 ^ topEntryBits (sOf n) P) pad (cliqueSummand M γ θ codec) =
      M.toClique (topWidth (sOf n) P + 1) γ θ := by
  dsimp only
  let P := K x (selectedLevel x n)
  let M := (clusterNet (sOf n) ℤ N k y (3 * P ^ 8) (by positivity)).coarsen
    (selectedTopGrouping (sOf n) N k y P hP)
  have hb := selected_explicit_code_bound (sOf n) N k y P hP γ
  have hpadded := (master_pad_suffices x C n pad hpad).2.2.1
  exact paddedRowSumZ_bound_eq_clique M γ θ ((selectedCodingZ (sOf n) N k y P hP).e γ) pad hb hpadded

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_explicit_code_bound
#print axioms IntegerFactoring.HalfAppendix.selected_explicit_row_as_padded_sum
