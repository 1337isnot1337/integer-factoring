/-
SPDX-License-Identifier: Apache-2.0

Executable cardinality of a selected top-group state space.
-/
import SelectedRadix

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedStateCount (x n : ℕ) (pad : List Unit) (γ : ℕ) : ℕ :=
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  cappedProduct (2 ^ topEntryBits (sOf n) P) (selectedRadix x n pad γ) g

set_option maxHeartbeats 2000000 in
theorem selectedStateCount_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n γ : α → ℕ} (hn : Unary n) (hγ : PolyT γ) :
    PolyT fun p : α × List Unit => selectedStateCount x (n p.1) p.2 (γ p.1) := by
  let β := α × List Unit
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let B : β → ℕ := fun p => topEntryBits (sOf (n p.1)) (P p)
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hB : Unary B := Unary.comp (selected_topEntryBits_unary x hx hn) PolyT.id'.fst
  have hr : PolyT fun p : β × ℕ => selectedRadix x (n p.1.1) p.1.2 (γ p.1.1) p.2 :=
    selectedRadix_polyT x hx hn hγ
  exact cappedProduct_polyT (B := fun p : β => 2 ^ B p) (n := g)
    (f := fun p i => selectedRadix x (n p.1) p.2 (γ p.1) i) (PolyT.two_pow hB) hg hr

theorem selectedStateCount_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    selectedStateCount x n pad γ.val = (κ.coarsen χ).s γ := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  have hrad : selectedRadix x n pad γ.val = Coding.sg κ χ γ := by
    funext i
    exact selectedRadix_eq x n N k y hP pad hpad γ i
  have hpos : ∀ i, 0 < Coding.sg κ χ γ i := Coding.sg_pos κ χ γ
  have hfull : pp (Coding.sg κ χ γ) g = (κ.coarsen χ).s γ := by
    change pp (Coding.sg κ χ γ) g = if h : γ.val < P then pp (Coding.sg κ χ ⟨γ.val, h⟩) g else 1
    rw [dif_pos γ.isLt]
  have hbound : (κ.coarsen χ).s γ ≤ 2 ^ topEntryBits s P :=
    selected_explicit_code_bound s N k y P hP γ
  unfold selectedStateCount
  dsimp only
  rw [hrad, cappedProduct_eq _ _ _ (by intro i _; exact hpos i) (by rw [range_prod_eq_pp, hfull]; exact hbound), range_prod_eq_pp, hfull]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedStateCount_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedStateCount_eq
