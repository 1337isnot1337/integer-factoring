/-
SPDX-License-Identifier: Apache-2.0

The explicit top evaluator returns the reference search's count after
removing divisor one. This remains a value theorem, not a costed factoring
machine.
-/
import BoundedTop
import HalfPrefix

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad IntegerFactoring.HalfPrefix

namespace IntegerFactoring.HalfAppendix

theorem evaluator_counts_reference_prefix :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l n N k y : ℕ) (hK : 2 ≤ K x l),
        let s := sOf n
        let P := K x l
        let L := topWidth s P + 1
        let B := topEntryBits s P
        let χ := selectedTopGrouping s N k y P hK
        let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
        let tabPad := List.replicate (2 ^ ((P - 1) * L)) ()
        2 ≤ N → N < 2 ^ n →
        run (List.replicate (2 ^ padExponent P L B C (2 ^ (l / 2))) (),
          List.replicate l (), writtenInstance tabPad P L (cliqueQueryFast T)) -
            (if 1 / 2 ^ k = y then 1 else 0) = (cntPre N k y : ℤ) := by
  obtain ⟨x, C, run, hx, hpoly, heval⟩ := evaluator_counts_explicit_top
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l n N k y hK
  dsimp only
  intro hN hNbits
  have hNlt : N < 2 ^ (pL (sOf n) * pL (sOf n)) :=
    hNbits.trans_le (Nat.pow_le_pow_right (by norm_num) (sOf_spec n))
  rw [heval l (sOf n) N k y hK (by omega) hNlt]
  have h := all_divisors_count N k y hN
  by_cases h1 : 1 / 2 ^ k = y
  · simp only [h1, ite_true] at h ⊢
    have h' : ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) = (cntPre N k y : ℤ) + 1 := by
      exact_mod_cast h
    omega
  · simp only [h1, ite_false, Nat.add_zero, sub_zero] at h ⊢
    exact_mod_cast h

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_reference_prefix
