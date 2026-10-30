/-
SPDX-License-Identifier: Apache-2.0

The Appendix evaluator applied to the tabulated, grouped factoring network.
This is a value theorem with explicit size assumptions, not a machine or time
bound for the constructor and full factoring composition.
-/
import ExplicitPad
import TopValue

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem evaluator_counts_top :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l L B : ℕ) (s N k y g : ℕ) (hg : 0 < g)
        (χ : Fin g → Fin (K x l)) (pad : U),
        1 ≤ L → 1 ≤ B → 2 ^ ((K x l - 1) * L) ≤ pad.length →
        (∀ a b, a < b → ((clusterNet s ℤ N k y g hg).coarsen χ).c a b ≤ L) →
        (∀ (v : Fin (K x l)) (θ : Fin (K x l) → List Bool), WFQ L v θ →
          (cliqueQuery (((clusterNet s ℤ N k y g hg).coarsen χ).toClique L)
            (v.1, List.ofFn θ)).natAbs ≤ 2 ^ B) →
        1 ≤ N → N < 2 ^ (pL s * pL s) →
        run (List.replicate (2 ^ padExponent (K x l) L B C (2 ^ (l / 2))) (),
          List.replicate l (),
          writtenInstance pad (K x l) L
            (cliqueQuery (((clusterNet s ℤ N k y g hg).coarsen χ).toClique L))) =
          ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  obtain ⟨x, C, run, hx, hpoly, heval⟩ := exists_explicit_padded_evaluator_ge
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l L B s N k y g hg χ pad hL hB hpad hfit hbound hN1 hN
  rw [heval l L B _ hL hB (writtenInstance_wfe pad _ _ _ _ hpad hL hbound)]
  exact tabulated_top_value s N k y g (K x l) L hg χ pad hpad hfit hN1 hN

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.evaluator_counts_top
