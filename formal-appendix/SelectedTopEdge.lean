/-
SPDX-License-Identifier: Apache-2.0

Executable edge width of the selected top grouping.
-/
import CrossWidthBound
import SelectedGroupingPoly

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedTopEdge (x n : ℕ) (pad : List Unit) (γ δ : ℕ) : ℕ :=
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  crossWidth g (selectedGroupDigit x n pad) (clusterEdge (pA s) (fw s) g) γ δ

set_option maxHeartbeats 2000000 in
theorem selectedTopEdge_unary (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n γ δ : α → ℕ} (hn : Unary n) (hγ : PolyT γ) (hδ : PolyT δ) :
    Unary fun p : α × List Unit => selectedTopEdge x (n p.1) p.2 (γ p.1) (δ p.1) := by
  let β := α × List Unit
  let A : β → ℕ := fun p => pA (sOf (n p.1))
  let F : β → ℕ := fun p => fw (sOf (n p.1))
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  have hA : Unary A := Unary.comp (selected_pA_unary hn) PolyT.id'.fst
  have hF : Unary F := Unary.comp (selected_fw_unary hn) PolyT.id'.fst
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hnodes : Unary fun p : β => gateCount (A p) :=
    Unary.add hA (Unary.add hA (Unary.add hA (Unary.add hA (Unary.mul hA hA))))
  have hW : Unary fun p : β => g p * g p * (gateCount (A p) * gateCount (A p) * F p) :=
    Unary.mul (Unary.mul hg hg) (Unary.mul (Unary.mul hnodes hnodes) hF)
  have hcolor : PolyT fun p : β × ℕ => selectedGroupDigit x (n p.1.1) p.1.2 p.2 :=
    selectedGroupDigit_polyT x hx hn
  have hedge : PolyT fun p : (β × ℕ) × ℕ => clusterEdge (A p.1.1) (F p.1.1) (g p.1.1) p.1.2 p.2 :=
    clusterEdge_polyT (Unary.comp hA PolyT.id'.fst.fst) (PolyT.comp hF.poly PolyT.id'.fst.fst)
      (PolyT.comp hg.poly PolyT.id'.fst.fst) PolyT.id'.fst.snd PolyT.id'.snd
  have hp : PolyT fun p : β => selectedTopEdge x (n p.1) p.2 (γ p.1) (δ p.1) :=
    crossWidth_polyT (g := g) (u := fun p : β => γ p.1) (v := fun p : β => δ p.1)
      (color := fun p i => selectedGroupDigit x (n p.1) p.2 i)
      (edge := fun p i j => clusterEdge (A p) (F p) (g p) i j)
      hg (PolyT.comp hγ PolyT.id'.fst) (PolyT.comp hδ PolyT.id'.fst) hcolor hedge
  apply Unary.of_le hW hp
  intro p
  exact crossWidth_le _ _ _ _ _ _ (fun i j => clusterEdge_le _ _ _ _ _)

theorem selectedTopEdge_symm (x n : ℕ) (pad : List Unit) (γ δ : ℕ) :
    selectedTopEdge x n pad γ δ = selectedTopEdge x n pad δ γ := by
  simp only [selectedTopEdge, crossWidth, or_comm]

theorem selectedTopEdge_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ δ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    selectedTopEdge x n pad γ.val δ.val = ((clusterNet (sOf n) ℤ N k y g (by positivity)).coarsen χ).c γ δ := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  apply crossWidth_eq_coarsen (clusterNet s ℤ N k y g (by positivity)) χ
    (selectedGroupDigit x n pad) (clusterEdge (pA s) (fw s) g)
  · intro a
    exact selectedGroupDigit_eq x n N k y hP pad hpad a
  · intro a b
    exact clusterEdge_eq s g N k y (by positivity) a b

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedTopEdge_unary
#print axioms IntegerFactoring.HalfAppendix.selectedTopEdge_eq

#print axioms IntegerFactoring.HalfAppendix.selectedTopEdge_symm
