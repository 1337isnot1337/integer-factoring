/-
SPDX-License-Identifier: Apache-2.0

Apply the explicit level selection to the evaluator-pad rate. There is no
machine or full runtime certificate.
-/
import LevelSelect
import SelectedPadCost

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

theorem selected_level_pad_sqrt_cost (ρ : ℝ) (x n C : ℕ)
    (hx : 1 ≤ x) (ha : AppRunningTime.PartA ρ x) (hn : 1 < n)
    (h0 : (K x 0) ^ 40 ≤ n) :
    let l := selectedLevel x n
    (padExponent (K x l) (topWidth (sOf n) (K x l) + 1)
      (topEntryBits (sOf n) (K x l)) C (2 ^ (l / 2)) : ℝ) ≤
      ((6001 + C * 500 : ℕ) : ℝ) *
        Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) * (n : ℝ) / Real.sqrt (Real.log n) := by
  have hs := selectedLevel_spec x n hx h0
  exact scheduled_pad_sqrt_cost ρ x (selectedLevel x n) n C hx ha hn hs.1 hs.2.le

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_level_pad_sqrt_cost
