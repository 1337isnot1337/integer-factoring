/-
SPDX-License-Identifier: Apache-2.0

A uniform unary bound for numeric coarsened edge widths.
-/
import ClusterEdge

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem crossWidth_le (g W u v : ℕ) (color : ℕ → ℕ) (edge : ℕ → ℕ → ℕ)
    (h : ∀ i j, edge i j ≤ W) : crossWidth g color edge u v ≤ g * g * W := by
  let row := fun i => (List.range g).map fun j =>
    if i < j ∧ ((color i = u ∧ color j = v) ∨ (color i = v ∧ color j = u)) then edge i j else 0
  let l := (List.range g).flatMap row
  have he : ∀ a ∈ l, a ≤ W := by
    intro a ha
    obtain ⟨i, hi, hm⟩ := List.mem_flatMap.mp ha
    obtain ⟨j, hj, rfl⟩ := List.mem_map.mp hm
    split <;> simp [h]
  have hlen : l.length = g * g := by
    simp [l, row, List.length_flatMap, List.length_map, List.length_range]
  have hsum := List.sum_le_card_nsmul l W he
  rw [hlen, Nat.nsmul_eq_mul] at hsum
  exact hsum

theorem clusterEdge_le (A width g u v : ℕ) : clusterEdge A width g u v ≤ gateCount A * gateCount A * width := by
  apply crossWidth_le
  intro i j
  unfold gateEdge
  split <;> simp

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.crossWidth_le
#print axioms IntegerFactoring.HalfAppendix.clusterEdge_le
