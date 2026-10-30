/-
SPDX-License-Identifier: Apache-2.0

The selected top evaluator with direct clique-row lookup. The reference
clique's entry bound and the complete constructor cost are still hypotheses
or open obligations; this theorem does not prove the factoring rate.
-/
import SelectedTop
import RowSimplify

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem evaluator_counts_selected_fast :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l s N k y B : ℕ) (hK : 2 ≤ K x l) (pad : U),
        let P := K x l
        let L := topWidth s P + 1
        let χ := selectedTopGrouping s N k y P hK
        let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
        1 ≤ B → 2 ^ ((P - 1) * L) ≤ pad.length →
        (∀ a z, (T a z).natAbs ≤ 2 ^ B) →
        1 ≤ N → N < 2 ^ (pL s * pL s) →
        run (List.replicate (2 ^ padExponent P L B C (2 ^ (l / 2))) (),
          List.replicate l (), writtenInstance pad P L (cliqueQueryFast T)) =
          ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  obtain ⟨x, C, run, hx, hpoly, heval⟩ := evaluator_counts_selected_top
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l s N k y B hK pad
  dsimp only
  intro hB hpad hbound hN1 hN
  rw [cliqueQueryFast_eq]
  apply heval l s N k y B hK pad hB hpad
  · intro a θ _
    exact cliqueQuery_bound _ hbound a θ
  · exact hN1
  · exact hN

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_selected_fast
