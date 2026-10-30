/-
SPDX-License-Identifier: Apache-2.0

Instantiate the selected evaluator with the proved integer entry bound. The
tabulation-pad length and the cost of constructing the input are not proved
here, nor is an end-to-end factoring machine or asymptotic time bound.
-/
import FastTop
import SelectedEntry

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

def topEntryBits (s K : ℕ) : ℕ :=
  Wc s (3 * K ^ 8) * (2 * (3 * K ^ 8) / K) + 1

theorem evaluator_counts_bounded_top :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l s N k y : ℕ) (hK : 2 ≤ K x l) (pad : U),
        let P := K x l
        let L := topWidth s P + 1
        let B := topEntryBits s P
        let χ := selectedTopGrouping s N k y P hK
        let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
        2 ^ ((P - 1) * L) ≤ pad.length →
        1 ≤ N → N < 2 ^ (pL s * pL s) →
        run (List.replicate (2 ^ padExponent P L B C (2 ^ (l / 2))) (),
          List.replicate l (), writtenInstance pad P L (cliqueQueryFast T)) =
          ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  obtain ⟨x, C, run, hx, hpoly, heval⟩ := evaluator_counts_selected_fast
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l s N k y hK pad
  dsimp only
  intro hpad hN1 hN
  apply heval l s N k y (topEntryBits s (K x l)) hK pad
  · simp [topEntryBits]
  · exact hpad
  · intro a z
    exact (selectedTop_entry_bound s N k y (K x l) (topWidth s (K x l) + 1) hK a z).trans
      (Nat.pow_le_pow_right (by norm_num) (by simp [topEntryBits]))
  · exact hN1
  · exact hN

/-- The top table and evaluator both use explicit, sufficient pads. This is
a value-level theorem; no cost certificate for writing the table is asserted. -/
theorem evaluator_counts_explicit_top :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l s N k y : ℕ) (hK : 2 ≤ K x l),
        let P := K x l
        let L := topWidth s P + 1
        let B := topEntryBits s P
        let χ := selectedTopGrouping s N k y P hK
        let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
        let tabPad := List.replicate (2 ^ ((P - 1) * L)) ()
        1 ≤ N → N < 2 ^ (pL s * pL s) →
        run (List.replicate (2 ^ padExponent P L B C (2 ^ (l / 2))) (),
          List.replicate l (), writtenInstance tabPad P L (cliqueQueryFast T)) =
          ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  obtain ⟨x, C, run, hx, hpoly, heval⟩ := evaluator_counts_bounded_top
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l s N k y hK
  dsimp only
  intro hN1 hN
  apply heval l s N k y hK _
  · simp
  · exact hN1
  · exact hN

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_bounded_top
#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_explicit_top
