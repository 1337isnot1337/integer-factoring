/-
SPDX-License-Identifier: Apache-2.0

Composition at the value interface. The table-value equality is an explicit
hypothesis: this file does not yet implement or verify the table constructor.
-/
import ExplicitPad
import Tabulate
import FactoringAux.MulCount

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg
open IntegerFactoring.HalfPad IntegerFactoring.HalfAppendix Factoring.Aux

namespace IntegerFactoring.HalfAppendix

/-- The actual Appendix evaluator computes the divisor count whenever the
written instance has the reference multiplication network's value. This
does not prove that such an instance can be constructed within the bound. -/
theorem evaluator_counts_of_value_bridge :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), PolyT run ∧
      ∀ (l lam B : ℕ) (I : PInst) (s N k y : ℕ),
        1 ≤ lam → 1 ≤ B → WFe (K x l) lam B I →
        1 ≤ N → N < 2 ^ (pL s * pL s) →
        zval (K x l) I = (mulNet s ℤ N k y).value →
        run (List.replicate (2 ^ padExponent (K x l) lam B C (2 ^ (l / 2))) (),
          List.replicate l (), I) =
            ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  obtain ⟨x, C, run, hpoly, heval⟩ := exists_explicit_padded_evaluator
  refine ⟨x, C, run, hpoly, ?_⟩
  intro l lam B I s N k y hlam hB hI hN1 hN hvalue
  rw [heval l lam B I hlam hB hI, hvalue]
  exact mulNet_value s N k y hN1 hN

/-- Tabulating bounded row functions supplies the required well-formed
written input. Equality of the row-function contraction with the original
multiplication network is still a hypothesis, as is a cost bound for F. -/
theorem evaluator_counts_tabulated_of_value_bridge :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), PolyT run ∧
      ∀ (l lam B : ℕ) (pad : U) (F : IQ → ℤ) (s N k y : ℕ),
        1 ≤ lam → 1 ≤ B → 2 ^ ((K x l - 1) * lam) ≤ pad.length →
        (∀ (v : Fin (K x l)) (θ : Fin (K x l) → List Bool), WFQ lam v θ →
          (F (v.1, List.ofFn θ)).natAbs ≤ 2 ^ B) →
        1 ≤ N → N < 2 ^ (pL s * pL s) →
        Inst.ZW lam [] (fun (v : Fin (K x l)) θ => F (v.1, List.ofFn θ)) =
          (mulNet s ℤ N k y).value →
        run (List.replicate (2 ^ padExponent (K x l) lam B C (2 ^ (l / 2))) (),
          List.replicate l (), writtenInstance pad (K x l) lam F) =
            ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  obtain ⟨x, C, run, hpoly, hcount⟩ := evaluator_counts_of_value_bridge
  refine ⟨x, C, run, hpoly, ?_⟩
  intro l lam B pad F s N k y hlam hB hpad hbound hN1 hN hvalue
  apply hcount l lam B (writtenInstance pad (K x l) lam F) s N k y
    hlam hB (writtenInstance_wfe pad _ _ _ F hpad hlam hbound) hN1 hN
  rw [writtenInstance_value pad _ _ F hpad]
  exact hvalue

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_of_value_bridge
#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_tabulated_of_value_bridge
