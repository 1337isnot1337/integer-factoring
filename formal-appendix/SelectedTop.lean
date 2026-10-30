/-
SPDX-License-Identifier: Apache-2.0

Specialize the top evaluator to the reference's deterministic first-good
grouping. Row bounds, tabulation pad, and all machine costs remain open.
-/
import TopGrouping
import TopEvaluator

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem evaluator_counts_selected_top :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l s N k y B : ℕ) (hK : 2 ≤ K x l) (pad : U),
        let P := K x l
        let L := topWidth s P + 1
        let χ := selectedTopGrouping s N k y P hK
        let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
        1 ≤ B → 2 ^ ((P - 1) * L) ≤ pad.length →
        (∀ (v : Fin P) (θ : Fin P → List Bool), WFQ L v θ →
          (cliqueQuery T (v.1, List.ofFn θ)).natAbs ≤ 2 ^ B) →
        1 ≤ N → N < 2 ^ (pL s * pL s) →
        run (List.replicate (2 ^ padExponent P L B C (2 ^ (l / 2))) (),
          List.replicate l (), writtenInstance pad P L (cliqueQuery T)) =
          ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  obtain ⟨x, C, run, hx, hpoly, heval⟩ := evaluator_counts_top
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l s N k y B hK pad
  dsimp only
  intro hB hpad hbound hN1 hN
  apply heval l (topWidth s (K x l) + 1) B s N k y (3 * K x l ^ 8)
    (by positivity) (selectedTopGrouping s N k y (K x l) hK) pad
  · omega
  · exact hB
  · exact hpad
  · exact selectedTopGrouping_fit s N k y (K x l) hK
  · exact hbound
  · exact hN1
  · exact hN

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_selected_top
