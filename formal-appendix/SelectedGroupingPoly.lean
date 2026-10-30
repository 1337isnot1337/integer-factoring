/-
SPDX-License-Identifier: Apache-2.0

Executable top-grouping search using a supplied unary pad.
-/
import SelectedGood
import GroupingSearchPoly

open Groundwork Levels.Gen Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def selectedCode (x n : ℕ) (pad : List Unit) : ℕ :=
  let P := K x (selectedLevel x n)
  let count := min (P ^ (3 * P ^ 8)) pad.length
  firstGood (fun code => selectedGood x n code = true) count count

theorem selectedCode_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : PolyT fun p : α × List Unit => selectedCode x (n p.1) p.2 := by
  let P : α → ℕ := fun a => K x (selectedLevel x (n a))
  let g : α → ℕ := fun a => 3 * P a ^ 8
  have hP : Unary P := selected_K_unary x hx hn
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hcount : PolyT fun a => P a ^ g a := PolyT.pow_unary hP.poly hg
  have hgood : PolyT fun p : α × ℕ => selectedGood x (n p.1) p.2 := selectedGood_polyT x hx hn
  have hmap : PolyT fun p : (α × List Unit) × ℕ => (p.1.1, p.2) :=
    PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd
  have htest : PolyT fun p : (α × List Unit) × ℕ => selectedGood x (n p.1.1) p.2 :=
    PolyT.comp (f := fun q : α × ℕ => selectedGood x (n q.1) q.2)
      (g := fun p : (α × List Unit) × ℕ => (p.1.1, p.2)) hgood hmap
  have hsearch : PolyT fun p : α × List Unit => firstGood (fun code => selectedGood x (n p.1) code = true)
      (min (P p.1 ^ g p.1) p.2.length) (min (P p.1 ^ g p.1) p.2.length) :=
    padded_firstGood_of_count_polyT (S := fun a => P a ^ g a)
      (ok := fun p code => selectedGood x (n p.1) code) hcount htest
  apply hsearch.congr
  intro p
  rfl

theorem selectedCode_eq_topGroupingCode (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length) :
    selectedCode x n pad = topGroupingCode (sOf n) N k y (K x (selectedLevel x n)) hP := by
  unfold selectedCode topGroupingCode
  dsimp only
  rw [Nat.min_eq_left hpad]
  have hp : (fun code => selectedGood x n code = true) =
      okG (clusterNet (sOf n) ℤ N k y (3 * K x (selectedLevel x n) ^ 8) (by positivity))
        (K x (selectedLevel x n)) (topWidth (sOf n) (K x (selectedLevel x n))) (by omega) := by
    funext code
    exact propext (selectedGood_eq_okG x n code N k y hP)
  simp only [hp]

def selectedGroupDigit (x n : ℕ) (pad : List Unit) (i : ℕ) : ℕ :=
  let P := K x (selectedLevel x n)
  cappedGroupingDigit P (3 * P ^ 8) (selectedCode x n pad) i

theorem selectedGroupDigit_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : PolyT fun p : (α × List Unit) × ℕ => selectedGroupDigit x (n p.1.1) p.1.2 p.2 := by
  let β := α × List Unit
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let code : β → ℕ := fun p => selectedCode x (n p.1) p.2
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hc : PolyT code := selectedCode_polyT x hx hn
  exact cappedGroupingDigit_polyT (K := P) (g := g) (code := code) hP.poly hg hc

theorem selectedGroupDigit_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (i : Fin (3 * (K x (selectedLevel x n)) ^ 8)) :
    selectedGroupDigit x n pad i.val = (selectedTopGrouping (sOf n) N k y (K x (selectedLevel x n)) hP i).val := by
  unfold selectedGroupDigit selectedTopGrouping
  dsimp only
  rw [selectedCode_eq_topGroupingCode x n N k y hP pad hpad]
  exact cappedGroupingDigit_chiOf _ _ _ (by omega) i

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedCode_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedCode_eq_topGroupingCode
#print axioms IntegerFactoring.HalfAppendix.selectedGroupDigit_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedGroupDigit_eq
