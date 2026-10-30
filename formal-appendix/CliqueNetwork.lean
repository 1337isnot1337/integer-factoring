/-
SPDX-License-Identifier: Apache-2.0

View the factoring reference's integer clique as an Appendix network.
-/
import FactoringAux.Clique
import Groundwork.Net.Net

open Finset

namespace IntegerFactoring.HalfAppendix

def cliqueNetwork {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) :
    Groundwork.Network (Fin K) (Fin K → Factoring.Aux.Bits L) (Factoring.Aux.Bits L) where
  S _ := univ
  W a z := if z a = Factoring.Aux.zeroBits L then T a z else 0
  Θ _ _ := univ
  f _ b z t := if z b = t then 1 else 0
  c _ _ := L

theorem cliqueNetwork_value {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) :
    (cliqueNetwork T).value = T.value := by
  rw [Factoring.Aux.Clique.value_eq_toNet]
  unfold Groundwork.Network.value Factoring.Aux.Net.value
  simp only [cliqueNetwork, Fintype.piFinset_univ]
  rfl

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.cliqueNetwork_value
