/-
SPDX-License-Identifier: Apache-2.0

Convert the reference clique network to Appendix B row functions on words.
This proves a value identity, not a bound on computing the row functions.
-/
import CliqueNetwork
import StepMachines.Math.SideCodes
import Tabulate

open Finset Groundwork Groundwork.Grid PaddedDescent.F

namespace IntegerFactoring.HalfAppendix

def cliqueRows {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) : Inst (Fin K) (List Bool) :=
  fun a θ => (cliqueNetwork T).toInst (fun _ _ => univ) (fun _ _ t => t)
    a (fun b => Groundwork.Grid.bitsOf L (θ b))

theorem cliqueRows_value {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) :
    Inst.ZW L [] (cliqueRows T) = T.value := by
  rw [Groundwork.Grid.ZW_eq_Z L (cliqueRows T) (Factoring.Aux.zeroBits L)]
  have hround (t : Factoring.Aux.Bits L) :
      Groundwork.Grid.bitsOf L ((bitsC L).enc t) = t := by
    exact (bitsC L).dec_enc (mem_univ t)
  have hinst : instW L (cliqueRows T) =
      (cliqueNetwork T).toInst (fun _ _ => univ) (fun _ _ t => t) := by
    funext a θ
    unfold instW cliqueRows Groundwork.Network.toInst
    apply sum_congr rfl
    intro z hz
    congr 1
    apply prod_congr rfl
    intro b hb
    have hne : b ≠ a := (mem_erase.mp hb).1
    simp [hne, hround]
  rw [hinst]
  rw [Groundwork.Network.Z_toInst]
  · exact cliqueNetwork_value T
  · intro a b hab
    exact Set.bijOn_id _

/-- Query a row by its numeric node and the list of its incident words. -/
def cliqueQuery {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) (q : Groundwork.Prog.IQ) : ℤ :=
  if h : q.1 < K then cliqueRows T ⟨q.1, h⟩ (fun u => q.2.getD u.1 []) else 0

theorem cliqueQuery_ofFn {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L)
    (v : Fin K) (θ : Fin K → List Bool) :
    cliqueQuery T (v.1, List.ofFn θ) = cliqueRows T v θ := by
  unfold cliqueQuery
  simp only [v.2, dite_true]
  congr 1
  funext u
  simp [List.getD_eq_getElem?_getD, u.2]

/-- The source tabulation writes a clique with exactly the reference value. -/
theorem tabulatedClique_value {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) (pad : U)
    (hpad : 2 ^ ((K - 1) * L) ≤ pad.length) :
    PaddedDescent.F.zval K (writtenInstance pad K L (cliqueQuery T)) = T.value := by
  rw [writtenInstance_value pad K L (cliqueQuery T) hpad]
  have h : (fun (v : Fin K) θ => cliqueQuery T (v.1, List.ofFn θ)) = cliqueRows T := by
    funext v θ
    exact cliqueQuery_ofFn T v θ
  rw [h, cliqueRows_value]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.cliqueRows_value
#print axioms IntegerFactoring.HalfAppendix.tabulatedClique_value
