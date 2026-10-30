/-
SPDX-License-Identifier: Apache-2.0

Executable top-grouping predicate for the reference multiplication network.
-/
import ClusterEdge
import GroupingPredicate

open Groundwork Levels.Gen Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def selectedGood (x n code : ℕ) : Bool :=
  let s := sOf n
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  goodCode P g (topWidth s P) code (clusterEdge (pA s) (fw s) g)

theorem selectedGood_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : PolyT fun p : α × ℕ => selectedGood x (n p.1) p.2 := by
  let β := α × ℕ
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let A : β → ℕ := fun p => pA (sOf (n p.1))
  let w : β → ℕ := fun p => fw (sOf (n p.1))
  let L : β → ℕ := fun p => topWidth (sOf (n p.1)) (P p)
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hA : Unary A := Unary.comp (selected_pA_unary hn) PolyT.id'.fst
  have hw : PolyT w := (Unary.comp (selected_fw_unary hn) PolyT.id'.fst).poly
  have hL : PolyT L := by
    have ht : PolyT fun p : β => topWidth (sOf (n p.1)) (P p) + 1 :=
      (Unary.comp (selected_topWidth_unary x hx hn) PolyT.id'.fst).poly
    exact (PolyT.sub ht (PolyT.const 1)).congr (fun p => by simp [L, P])
  have he : PolyT fun q : (β × ℕ) × ℕ => clusterEdge (A q.1.1) (w q.1.1) (g q.1.1) q.1.2 q.2 :=
    clusterEdge_polyT (Unary.comp hA PolyT.id'.fst.fst) (PolyT.comp hw PolyT.id'.fst.fst)
      (PolyT.comp hg.poly PolyT.id'.fst.fst) PolyT.id'.fst.snd PolyT.id'.snd
  exact goodCode_polyT (K := P) (g := g) (L := L) (code := fun p : β => p.2)
    (edge := fun p i j => clusterEdge (A p) (w p) (g p) i j) hP hg hL PolyT.id'.snd he

theorem selectedGood_eq_okG (x n code N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n)) :
    selectedGood x n code = true ↔
      okG (clusterNet (sOf n) ℤ N k y (3 * K x (selectedLevel x n) ^ 8) (by positivity))
        (K x (selectedLevel x n)) (topWidth (sOf n) (K x (selectedLevel x n))) (by omega) code := by
  unfold selectedGood
  apply goodCode_eq_okG
  intro i j
  exact clusterEdge_eq _ _ N k y (by positivity) i j

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedGood_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedGood_eq_okG
