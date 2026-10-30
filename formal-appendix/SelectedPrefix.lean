/-
SPDX-License-Identifier: Apache-2.0

Concrete selected-level prefix counter, using a single supplied master pad.
-/
import SelectedTable
import PureFactor

open Groundwork Groundwork.Prog PaddedDescent.F Factoring.Aux Levels.Gen
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

def selectedPrefix (x C : ℕ) (run : U × U × PaddedDescent.F.PInst → ℤ)
    (n : ℕ) (pad : List Unit) (N k y : ℕ) : ℕ :=
  let l := selectedLevel x n
  let P := K x l
  let L := topWidth (sOf n) P + 1
  let B := topEntryBits (sOf n) P
  let ep := pad.take (2 ^ padExponent P L B C (2 ^ (l / 2)))
  let depth := List.replicate l ()
  (run (ep, depth, selectedTable x n pad N k y) -
    (if (if k = 0 then 1 else 0) = y then 1 else 0)).toNat

set_option maxHeartbeats 2000000 in
theorem selectedPrefix_polyT (x C : ℕ) (hx : 1 ≤ x) (run : U × U × PaddedDescent.F.PInst → ℤ) (hrun : PolyT run)
    {α : Type} [Enc α] {n N k y : α → ℕ} (hn : Unary n) (hN : PolyT N) (hk : PolyT k) (hy : PolyT y) :
    PolyT fun p : α × List Unit => selectedPrefix x C run (n p.1) p.2 (N p.1) (k p.1) (y p.1) := by
  let β := α × List Unit
  let e : α → ℕ := fun a => padExponent (K x (selectedLevel x (n a)))
    (topWidth (sOf (n a)) (K x (selectedLevel x (n a))) + 1)
    (topEntryBits (sOf (n a)) (K x (selectedLevel x (n a)))) C (2 ^ (selectedLevel x (n a) / 2))
  have he : Unary e := selected_padExponent_unary x C hx hn
  have hep : PolyT fun p : β => p.2.take (2 ^ e p.1) := bounded_exp_pad_polyT he
  have hd : PolyT fun p : β => List.replicate (selectedLevel x (n p.1)) () :=
    Unary.comp (selectedLevel_unary x hn) PolyT.id'.fst
  have ht : PolyT fun p : β => selectedTable x (n p.1) p.2 (N p.1) (k p.1) (y p.1) :=
    selectedTable_polyT x hx (n := n) (N := N) (k := k) (y := y) hn hN hk hy
  have hr : PolyT fun p : β => run (p.2.take (2 ^ e p.1), List.replicate (selectedLevel x (n p.1)) (), selectedTable x (n p.1) p.2 (N p.1) (k p.1) (y p.1)) :=
    PolyT.comp hrun (PolyT.prodMk hep (PolyT.prodMk hd ht))
  have hk' : PolyT fun p : β => k p.1 := PolyT.comp hk PolyT.id'.fst
  have hy' : PolyT fun p : β => y p.1 := PolyT.comp hy PolyT.id'.fst
  have hbit : PolyT fun p : β => if k p.1 = 0 then (1 : ℕ) else 0 :=
    PolyT.ite' (PolyT.decide_eq hk' (PolyT.const 0)) (PolyT.const 1) (PolyT.const 0)
  have hcorr : PolyT fun p : β => if (if k p.1 = 0 then 1 else 0) = y p.1 then (1 : ℤ) else 0 :=
    PolyT.ite' (PolyT.decide_eq hbit hy') (PolyT.const (1 : ℤ)) (PolyT.const 0)
  unfold selectedPrefix
  with_reducible exact PolyT.toNat (PolyT.subZ hr hcorr)

private theorem one_div_pow (k : ℕ) : 1 / 2 ^ k = if k = 0 then 1 else 0 := by
  cases k with
  | zero => simp
  | succ k =>
    have h : 1 < 2 ^ (k + 1) := by
      rw [pow_succ]
      have hpos : 1 ≤ 2 ^ k := Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
      omega
    simp [Nat.div_eq_of_lt h]

set_option maxHeartbeats 2000000 in
theorem selectedPrefix_eq (x C n N k y : ℕ) (run : U × U × PaddedDescent.F.PInst → ℤ)
    (hP : 2 ≤ K x (selectedLevel x n)) (pad : List Unit)
    (hpad : 2 ^ masterExponent x C n ≤ pad.length) :
    selectedPrefix x C run n pad N k y = prefixFromEvaluator x C run (selectedLevel x n) n N hP k y := by
  let l := selectedLevel x n
  let P := K x l
  let L := topWidth (sOf n) P + 1
  let B := topEntryBits (sOf n) P
  let e := padExponent P L B C (2 ^ (l / 2))
  have he : pad.take (2 ^ e) = List.replicate (2 ^ e) () :=
    bounded_exp_pad_eq e pad (master_pad_suffices x C n pad hpad).2.2.2
  have ht := selectedTable_eq x C n N k y hP pad hpad
  dsimp only at ht
  unfold selectedPrefix prefixFromEvaluator
  dsimp only
  rw [he, ht, one_div_pow]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedPrefix_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedPrefix_eq
