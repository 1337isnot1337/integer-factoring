/-
SPDX-License-Identifier: Apache-2.0

The reference coding bounds the number of settings in each selected group.
The bound concerns the finite types only; it is not yet a table-entry bound.
-/
import TopGrouping
import FactoringAux.EngineCostTop

open Finset Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem selectedTop_state_card (s N k y K : ℕ) (hK : 2 ≤ K) (a : Fin K) :
    Fintype.card (((clusterNet s ℤ N k y (3 * K ^ 8) (by positivity)).coarsen
      (selectedTopGrouping s N k y K hK)).S a) ≤
        2 ^ (Wc s (3 * K ^ 8) * (2 * (3 * K ^ 8) / K)) := by
  let g := 3 * K ^ 8
  have hg : 0 < g := by dsimp [g]; positivity
  let χ := selectedTopGrouping s N k y K hK
  let κ := (mulCoding s 2 N k y).coarsen (fun a => clMG s g hg ((MG.enum (pA s)).symm a))
  have hsize : (univ.filter fun b => χ b = a).card ≤ 2 * g / K := by
    rw [Nat.le_div_iff_mul_le (by omega : 0 < K), mul_comm]
    exact (topGroupingCode_good s N k y K hK).2.1 a
  have hb := coarsen_s_le_unif κ χ a (Wc s g) (2 * g / K)
    (fun b => cluster_s_le s 2 N k y g hg b) hsize
  have hcard : Fintype.card (((clusterNet s ℤ N k y g hg).coarsen χ).S a) = (κ.coarsen χ).s a := by
    have heq : ((clusterNet s ℤ N k y g hg).coarsen χ).S a =
        ((clusterNet s (ZMod 2) N k y g hg).coarsen χ).S a := by rfl
    let e := (Equiv.cast heq).trans ((κ.coarsen χ).e a)
    have h := @Fintype.card_congr _ _ (((clusterNet s ℤ N k y g hg).coarsen χ).fS a)
      (inferInstance : Fintype (Fin ((κ.coarsen χ).s a))) e
    simpa only [Fintype.card_fin] using h
  simpa [g, χ] using hcard.trans_le hb

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedTop_state_card
