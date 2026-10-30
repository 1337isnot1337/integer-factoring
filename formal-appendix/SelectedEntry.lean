/-
SPDX-License-Identifier: Apache-2.0

Entry bounds for the selected grouped multiplication network.
-/
import BinaryCoarsen
import EntryBound
import TopStates

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem is01_bounds {z : ℤ} (h : Is01 z) : 0 ≤ z ∧ z ≤ 1 := by
  rcases h with rfl | rfl <;> omega

theorem padF_01 {K L : ℕ} (N : Net ℤ K) (hN : Net01 N) (a b : Fin K) :
    Fac01 (N.padF L a b) := by
  unfold Net.padF
  split_ifs with h
  · exact fac_pad_01 (N.F a b) L h (hN.factors a b)
  · constructor <;> intro x t <;> exact Or.inl rfl

theorem selectedTop_01 (s N k y K : ℕ) (hK : 2 ≤ K) :
    Net01 ((clusterNet s ℤ N k y (3 * K ^ 8) (by positivity)).coarsen
      (selectedTopGrouping s N k y K hK)) := by
  have h := coarsen_01 (mulNet s ℤ N k y) (multiplicationNet_01 s N k y)
    (fun a => clMG s (3 * K ^ 8) (by positivity) ((MG.enum (pA s)).symm a))
  have h' := coarsen_01 _ h (selectedTopGrouping s N k y K hK)
  exact h'

theorem selectedTop_entry_bound (s N k y K L : ℕ) (hK : 2 ≤ K)
    (a : Fin K) (θ : Fin K → Bits L) :
    ((((clusterNet s ℤ N k y (3 * K ^ 8) (by positivity)).coarsen
      (selectedTopGrouping s N k y K hK)).toClique L) a θ).natAbs ≤
        2 ^ (Wc s (3 * K ^ 8) * (2 * (3 * K ^ 8) / K)) := by
  let M := (clusterNet s ℤ N k y (3 * K ^ 8) (by positivity)).coarsen
    (selectedTopGrouping s N k y K hK)
  have hM : Net01 M := selectedTop_01 s N k y K hK
  have hw : ∀ a z, 0 ≤ M.W a z ∧ M.W a z ≤ 1 := fun a z => is01_bounds (hM.weights a z)
  have hl : ∀ a b z t, 0 ≤ (M.padF L a b).l z t ∧ (M.padF L a b).l z t ≤ 1 :=
    fun a b z t => is01_bounds ((padF_01 M hM a b).1 z t)
  have hr : ∀ a b z t, 0 ≤ (M.padF L a b).r z t ∧ (M.padF L a b).r z t ≤ 1 :=
    fun a b z t => is01_bounds ((padF_01 M hM a b).2 z t)
  exact (clique_entry_le_card M hw hl hr a θ).trans (selectedTop_state_card s N k y K hK a)

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedTop_entry_bound
