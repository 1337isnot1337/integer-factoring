/-
SPDX-License-Identifier: Apache-2.0

Exact value of the tabulated, grouped multiplication network. The grouping
and pad-size hypotheses remain explicit, as do construction and cost bounds.
-/
import CliqueRows
import FactoringAux.Cluster

open Groundwork PaddedDescent.F Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem tabulated_top_value (s N k y g K L : ℕ) (hg : 0 < g)
    (χ : Fin g → Fin K) (pad : U)
    (hpad : 2 ^ ((K - 1) * L) ≤ pad.length)
    (hfit : ∀ a b, a < b → ((clusterNet s ℤ N k y g hg).coarsen χ).c a b ≤ L)
    (hN1 : 1 ≤ N) (hN : N < 2 ^ (pL s * pL s)) :
    zval K (writtenInstance pad K L
      (cliqueQuery (((clusterNet s ℤ N k y g hg).coarsen χ).toClique L))) =
        ((N.divisors.filter fun d => d / 2 ^ k = y).card : ℤ) := by
  rw [tabulatedClique_value _ pad hpad]
  rw [Net.value_group _ χ L hfit]
  exact clusterNet_value s ℤ N k y g hg hN1 hN

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.tabulated_top_value
