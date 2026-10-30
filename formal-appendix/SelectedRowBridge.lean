/-
SPDX-License-Identifier: Apache-2.0

The selected top clique row is a padded sum over its bounded state space.
The chosen finite-state codec and summand are not certified executable here.
-/
import MasterPad
import RowSumBridge

open Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

theorem selected_row_as_padded_sum (x C n N k y : ℕ)
    (hP : 2 ≤ K x (selectedLevel x n)) (pad : List Unit)
    (hpad : 2 ^ masterExponent x C n ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n)))
    (θ : Fin (K x (selectedLevel x n)) → Bits (topWidth (sOf n) (K x (selectedLevel x n)) + 1)) :
    let P := K x (selectedLevel x n)
    let M := (clusterNet (sOf n) ℤ N k y (3 * P ^ 8) (by positivity)).coarsen
      (selectedTopGrouping (sOf n) N k y P hP)
    let codec := Fintype.equivFin (M.S γ)
    paddedRowSumZ (2 ^ topEntryBits (sOf n) P) pad (cliqueSummand M γ θ codec) =
      M.toClique (topWidth (sOf n) P + 1) γ θ := by
  dsimp only
  let P := K x (selectedLevel x n)
  let M := (clusterNet (sOf n) ℤ N k y (3 * P ^ 8) (by positivity)).coarsen
    (selectedTopGrouping (sOf n) N k y P hP)
  have hstates : Fintype.card (M.S γ) ≤ 2 ^ topEntryBits (sOf n) P := by
    have h := selectedTop_state_card (sOf n) N k y P hP γ
    have hp : Wc (sOf n) (3 * P ^ 8) * (2 * (3 * P ^ 8) / P) ≤ topEntryBits (sOf n) P := by
      simp [topEntryBits]
    exact h.trans (Nat.pow_le_pow_right (by norm_num) hp)
  have hpadded := (master_pad_suffices x C n pad hpad).2.2.1
  exact paddedRowSumZ_bound_eq_clique M γ θ (Fintype.equivFin (M.S γ)) pad hstates hpadded

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_row_as_padded_sum
