/-
SPDX-License-Identifier: Apache-2.0

Bound the Appendix evaluator fragment's own steps from the selected
instance-length envelope. Writing the instance and composing the factoring
machine are separate obligations.
-/
import InputSize
import FragmentCost

open Groundwork Groundwork.TM PaddedDescent.F Levels.Gen Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem evaluator_steps_at_schedule (ρ : ℝ) (x C : ℕ) (run : U × U × PInst → ℤ)
    (hrun : PolyT run) :
    ∃ (Pgm : Fr) (D : ℕ), ∀ (l n : ℕ) (I : PInst),
      1 ≤ x → AppRunningTime.PartA ρ x →
      (K x l) ^ 40 ≤ n → n ≤ (K x (l + 1)) ^ 40 →
      len I + 2 ≤ PaddedDescent.F.instLen (K x l)
        (topWidth (sOf n) (K x l) + 1) (topEntryBits (sOf n) (K x l)) →
      let r := 2 ^ (l / 2)
      let p := padExponent (K x l) (topWidth (sOf n) (K x l) + 1)
        (topEntryBits (sOf n) (K x l)) C r
      let input := (List.replicate (2 ^ p) (), List.replicate l (), I)
      ∀ (M S T N : List Sym), ∃ t ≤ 2 ^ (D * (n / r + 1)),
        Pgm.Runs none (st (wd input ++ M) S T N) t (st (wd (run input) ++ M) S T N) := by
  obtain ⟨Pgm, D, hfrag⟩ := polyT_fragment_exp run hrun (13006 + C * 500)
  refine ⟨Pgm, D, ?_⟩
  intro l n I hx ha hnlo hnhi hI
  dsimp only
  intro M S T N
  apply hfrag _ (n / 2 ^ (l / 2) + 1) (Nat.le_add_left 1 _) _ M S T N
  exact scheduled_evaluator_input_len ρ x l n C I hx ha hnlo hnhi hI

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.evaluator_steps_at_schedule
