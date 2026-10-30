/-
SPDX-License-Identifier: Apache-2.0

Zero-one preservation for grouping factorizations and weights.
-/
import BinaryNet

open Finset Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem prod_01 {α : Type} (s : Finset α) (f : α → ℤ) (h : ∀ a ∈ s, Is01 (f a)) :
    Is01 (∏ a ∈ s, f a) :=
  Finset.prod_induction f Is01 (fun _ _ => is01_mul) (Or.inr rfl) h

theorem orient_01 {g K : ℕ} (N : Net ℤ g) (hN : Net01 N)
    (χ : Fin g → Fin K) (a b : Fin K) (p : Fin g × Fin g) :
    Fac01 (N.orient χ a b p) := by
  constructor
  · intro z t
    unfold Net.orient
    split_ifs <;> first
      | exact Or.inl rfl
      | exact (hN.factors p.1 p.2).1 _ _
      | exact (hN.factors p.1 p.2).2 _ _
  · intro z t
    unfold Net.orient
    split_ifs <;> first
      | exact Or.inl rfl
      | exact (hN.factors p.1 p.2).1 _ _
      | exact (hN.factors p.1 p.2).2 _ _

theorem coarsen_factors_01 {g K : ℕ} (N : Net ℤ g) (hN : Net01 N)
    (χ : Fin g → Fin K) (a b : Fin K) : Fac01 ((N.coarsen χ).F a b) := by
  change Fac01 (Fac.prodList _ _ _)
  exact fac_prodList_01 _ _ (fun p => orient_01 N hN χ a b p) _

theorem coarsen_weights_01 {g K : ℕ} (N : Net ℤ g) (hN : Net01 N)
    (χ : Fin g → Fin K) (a : Fin K) (z : (N.coarsen χ).S a) :
    Is01 ((N.coarsen χ).W a z) := by
  change Is01 ((∏ i : {i // χ i = a}, N.W i (z i)) *
    ∏ p ∈ pairs g, (if h : χ p.1 = a ∧ χ p.2 = a then N.Ψ p.1 p.2 (z ⟨p.1, h.1⟩) (z ⟨p.2, h.2⟩) else 1))
  apply is01_mul
  · exact prod_01 univ _ (fun i _ => hN.weights i (z i))
  · apply prod_01
    intro p hp
    split_ifs with h
    · exact hN.interactions p.1 p.2 (mem_pairs.mp hp) _ _
    · exact Or.inr rfl

theorem orient_eval_01 {g K : ℕ} (N : Net ℤ g) (hN : Net01 N)
    (χ : Fin g → Fin K) (a b : Fin K) (hab : a < b) (p : Fin g × Fin g)
    (hp : p ∈ Net.crossList χ a b)
    (z : (N.coarsen χ).S a) (z' : (N.coarsen χ).S b) :
    Is01 ((N.orient χ a b p).eval z z') := by
  have hp' : p.1 < p.2 ∧ Net.Cross χ a b p := by
    simpa [Net.crossList, mem_pairsList] using hp
  rcases hp'.2 with ⟨h1, h2⟩ | ⟨h1, h2⟩
  · have h := hN.interactions p.1 p.2 hp'.1 (z ⟨p.1, h1⟩) (z' ⟨p.2, h2⟩)
    simpa [Net.orient, Fac.eval, Net.Ψ, h1, h2] using h
  · have h := hN.interactions p.1 p.2 hp'.1 (z' ⟨p.1, h1⟩) (z ⟨p.2, h2⟩)
    have heq : (N.orient χ a b p).eval z z' = N.Ψ p.1 p.2 (z' ⟨p.1, h1⟩) (z ⟨p.2, h2⟩) := by
      simp only [Net.orient, h2, h1, hab.ne, hab.ne', dif_pos, Fac.eval, Net.Ψ]
      exact Finset.sum_congr rfl fun _ _ => mul_comm _ _
    rw [heq]
    exact h

theorem coarsen_interactions_01 {g K : ℕ} (N : Net ℤ g) (hN : Net01 N)
    (χ : Fin g → Fin K) (a b : Fin K) (hab : a < b)
    (z : (N.coarsen χ).S a) (z' : (N.coarsen χ).S b) :
    Is01 ((N.coarsen χ).Ψ a b z z') := by
  change Is01 (((N.coarsen χ).F a b).eval z z')
  rw [Net.coarsen_F_eval]
  apply List.prod_induction Is01 (fun _ _ => is01_mul) (Or.inr rfl)
  intro v hv
  obtain ⟨p, hp, rfl⟩ := List.mem_map.mp hv
  exact orient_eval_01 N hN χ a b hab p hp z z'

theorem coarsen_01 {g K : ℕ} (N : Net ℤ g) (hN : Net01 N)
    (χ : Fin g → Fin K) : Net01 (N.coarsen χ) :=
  ⟨coarsen_weights_01 N hN χ, coarsen_factors_01 N hN χ,
    coarsen_interactions_01 N hN χ⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.coarsen_factors_01
#print axioms IntegerFactoring.HalfAppendix.coarsen_weights_01
#print axioms IntegerFactoring.HalfAppendix.coarsen_01
