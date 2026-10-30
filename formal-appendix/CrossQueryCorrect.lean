/-
SPDX-License-Identifier: Apache-2.0

Value bridge between a numeric crossing scan and the reference coarsening.
The numeric edge-width oracle itself remains to be implemented.
-/
import CrossQueryPoly
import FactoringAux.EngineGroup

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

private theorem sum_filter_map {α : Type} (xs : List α) (p : α → Bool) (f : α → ℕ) :
    ((xs.filter p).map f).sum = (xs.map fun x => if p x then f x else 0).sum := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
      simp only [List.filter_cons]
      cases h : p x <;> simp [h, ih]

theorem crossWidth_eq_coarsen {R : Type} [CommRing R] {g K : ℕ}
    (N : Net R g) (χ : Fin g → Fin K) (color : ℕ → ℕ) (edge : ℕ → ℕ → ℕ)
    (hc : ∀ i : Fin g, color i.val = (χ i).val)
    (he : ∀ i j : Fin g, edge i.val j.val = N.c i j) (u v : Fin K) :
    crossWidth g color edge u.val v.val = (N.coarsen χ).c u v := by
  unfold crossWidth Net.coarsen Net.crossList pairsList
  rw [← List.map_coe_finRange_eq_range]
  simp only [List.flatMap_map, List.map_map, Function.comp_def]
  rw [sum_filter_map, sum_filter_map]
  simp only [SProd.sprod]
  unfold List.product
  rw [List.map_flatMap]
  simp only [List.map_map, Function.comp_def]
  simp [hc, he, Net.Cross, Fin.ext_iff, ite_and]

theorem crossingOK_eq {R : Type} [CommRing R] {g K : ℕ}
    (N : Net R g) (χ : Fin g → Fin K) (color : ℕ → ℕ) (edge : ℕ → ℕ → ℕ) (L : ℕ)
    (hc : ∀ i : Fin g, color i.val = (χ i).val)
    (he : ∀ i j : Fin g, edge i.val j.val = N.c i j) :
    crossingOK K g L color edge = true ↔
      ∀ u v : Fin K, u < v → (N.coarsen χ).c u v ≤ L := by
  simp only [crossingOK, List.all_eq_true, List.mem_range, decide_eq_true_eq]
  constructor
  · intro h u v huv
    have h' := h u.val u.isLt v.val v.isLt huv
    rwa [crossWidth_eq_coarsen N χ color edge hc he u v] at h'
  · intro h u hu v hv huv
    have h' := h ⟨u, hu⟩ ⟨v, hv⟩ huv
    rw [crossWidth_eq_coarsen N χ color edge hc he ⟨u, hu⟩ ⟨v, hv⟩]
    exact h'

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.crossWidth_eq_coarsen
#print axioms IntegerFactoring.HalfAppendix.crossingOK_eq
