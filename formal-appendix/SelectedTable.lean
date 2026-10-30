/-
SPDX-License-Identifier: Apache-2.0

Executable selected clique table, using the supplied master pad.
-/
import SelectedRow
import WrittenPoly

open Groundwork Groundwork.Prog Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedTable (x n : ℕ) (pad : List Unit) (N k y : ℕ) : PaddedDescent.F.PInst :=
  let P := K x (selectedLevel x n)
  let L := topWidth (sOf n) P + 1
  let tabPad := pad.take (2 ^ ((P - 1) * L))
  writtenInstance tabPad P L (fun q => selectedRow x n pad N k y q.1 q.2)

set_option maxHeartbeats 2000000 in
theorem selectedTable_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n N k y : α → ℕ} (hn : Unary n) (hN : PolyT N) (hk : PolyT k) (hy : PolyT y) :
    PolyT fun p : α × List Unit => selectedTable x (n p.1) p.2 (N p.1) (k p.1) (y p.1) := by
  let β := α × List Unit
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let L : β → ℕ := fun p => topWidth (sOf (n p.1)) (P p) + 1
  let e : α → ℕ := fun a => (K x (selectedLevel x (n a)) - 1) * (topWidth (sOf (n a)) (K x (selectedLevel x (n a))) + 1)
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hL : Unary L := Unary.comp (selected_topWidth_unary x hx hn) PolyT.id'.fst
  have he : Unary e := selected_tableExponent_unary x hx hn
  have ht : PolyT fun p : β => p.2.take (2 ^ e p.1) := bounded_exp_pad_polyT he
  have hrowbase : PolyT fun p : (α × IQ) × List Unit => selectedRow x (n p.1.1) p.2
      (N p.1.1) (k p.1.1) (y p.1.1) p.1.2.1 p.1.2.2 :=
    selectedRow_polyT x hx (n := fun p : α × IQ => n p.1) (N := fun p : α × IQ => N p.1)
      (k := fun p : α × IQ => k p.1) (y := fun p : α × IQ => y p.1)
      (γ := fun p : α × IQ => p.2.1) (row := fun p : α × IQ => p.2.2)
      (Unary.comp hn PolyT.id'.fst) (PolyT.comp hN PolyT.id'.fst) (PolyT.comp hk PolyT.id'.fst)
      (PolyT.comp hy PolyT.id'.fst) PolyT.id'.snd.fst PolyT.id'.snd.snd
  have hrow : PolyT fun p : β × IQ => selectedRow x (n p.1.1) p.1.2 (N p.1.1) (k p.1.1) (y p.1.1) p.2.1 p.2.2 :=
    PolyT.comp (f := fun p : (α × IQ) × List Unit => selectedRow x (n p.1.1) p.2
      (N p.1.1) (k p.1.1) (y p.1.1) p.1.2.1 p.1.2.2)
      (g := fun p : β × IQ => ((p.1.1, p.2), p.1.2)) hrowbase
      (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd) PolyT.id'.fst.snd)
  unfold selectedTable
  with_reducible
    exact writtenInstance_polyT (pad := fun p : β => p.2.take (2 ^ e p.1)) (K := P) (L := L)
      (F := fun p q => selectedRow x (n p.1) p.2 (N p.1) (k p.1) (y p.1) q.1 q.2) ht hP hL hrow

set_option maxHeartbeats 2000000 in
theorem selectedTable_eq (x C n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : 2 ^ masterExponent x C n ≤ pad.length) :
    let P := K x (selectedLevel x n)
    let L := topWidth (sOf n) P + 1
    let M := (clusterNet (sOf n) ℤ N k y (3 * P ^ 8) (by positivity)).coarsen
      (selectedTopGrouping (sOf n) N k y P hP)
    selectedTable x n pad N k y = writtenInstance (List.replicate (2 ^ ((P - 1) * L)) ()) P L (cliqueQueryFast (M.toClique L)) := by
  dsimp only
  let P := K x (selectedLevel x n)
  let L := topWidth (sOf n) P + 1
  let M := (clusterNet (sOf n) ℤ N k y (3 * P ^ 8) (by positivity)).coarsen
    (selectedTopGrouping (sOf n) N k y P hP)
  have ht : pad.take (2 ^ ((P - 1) * L)) = List.replicate (2 ^ ((P - 1) * L)) () :=
    bounded_exp_pad_eq _ pad (master_pad_suffices x C n pad hpad).2.1
  have hrow : (fun q : IQ => selectedRow x n pad N k y q.1 q.2) = cliqueQueryFast (M.toClique L) := by
    funext q
    have h := selectedRow_eq x C n N k y hP pad hpad q.1 q.2
    dsimp only at h
    with_reducible exact h
  unfold selectedTable
  dsimp only
  rw [ht, hrow]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedTable_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedTable_eq
