/-
SPDX-License-Identifier: Apache-2.0

Executable selected clique-row query, with padded state enumeration.
-/
import SelectedSummand
import RowSimplify

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedRow (x n : ℕ) (pad : List Unit) (N k y γ : ℕ) (row : List (List Bool)) : ℤ :=
  let P := K x (selectedLevel x n)
  if γ < P then paddedRowSumZ (2 ^ topEntryBits (sOf n) P) pad
    (fun code => selectedSummand x n pad N k y γ code row)
  else 0

set_option maxHeartbeats 2000000 in
theorem selectedRow_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n N k y γ : α → ℕ} {row : α → List (List Bool)}
    (hn : Unary n) (hN : PolyT N) (hk : PolyT k) (hy : PolyT y)
    (hγ : PolyT γ) (hrow : PolyT row) :
    PolyT fun p : α × List Unit => selectedRow x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (row p.1) := by
  let β := α × List Unit
  let P : α → ℕ := fun p => K x (selectedLevel x (n p))
  let B : α → ℕ := fun p => topEntryBits (sOf (n p)) (P p)
  have hP : Unary P := selected_K_unary x hx hn
  have hB : Unary B := selected_topEntryBits_unary x hx hn
  have hS : PolyT fun p : α => 2 ^ B p := PolyT.two_pow hB
  have hbase : PolyT fun p : (α × ℕ) × List Unit => selectedSummand x (n p.1.1) p.2
      (N p.1.1) (k p.1.1) (y p.1.1) (γ p.1.1) p.1.2 (row p.1.1) :=
    selectedSummand_polyT x hx (n := fun p : α × ℕ => n p.1) (N := fun p : α × ℕ => N p.1)
      (k := fun p : α × ℕ => k p.1) (y := fun p : α × ℕ => y p.1) (γ := fun p : α × ℕ => γ p.1)
      (code := fun p : α × ℕ => p.2) (row := fun p : α × ℕ => row p.1)
      (Unary.comp hn PolyT.id'.fst) (PolyT.comp hN PolyT.id'.fst) (PolyT.comp hk PolyT.id'.fst)
      (PolyT.comp hy PolyT.id'.fst) (PolyT.comp hγ PolyT.id'.fst) PolyT.id'.snd (PolyT.comp hrow PolyT.id'.fst)
  have hsummand : PolyT fun p : β × ℕ => selectedSummand x (n p.1.1) p.1.2
      (N p.1.1) (k p.1.1) (y p.1.1) (γ p.1.1) p.2 (row p.1.1) :=
    PolyT.comp (f := fun p : (α × ℕ) × List Unit => selectedSummand x (n p.1.1) p.2
      (N p.1.1) (k p.1.1) (y p.1.1) (γ p.1.1) p.1.2 (row p.1.1))
      (g := fun p : β × ℕ => ((p.1.1, p.2), p.1.2)) hbase
      (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd) PolyT.id'.fst.snd)
  have hsum : PolyT fun p : β => paddedRowSumZ (2 ^ B p.1) p.2
      (fun code => selectedSummand x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) code (row p.1)) :=
    paddedRowSumZ_polyT (S := fun p : α => 2 ^ B p)
      (f := fun p code => selectedSummand x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) code (row p.1)) hS hsummand
  have hguard : PolyT fun p : β => decide (γ p.1 < P p.1) :=
    PolyT.decide_lt (PolyT.comp hγ PolyT.id'.fst) (PolyT.comp hP.poly PolyT.id'.fst)
  unfold selectedRow
  with_reducible
    exact PolyT.ite' (c := fun p : β => γ p.1 < P p.1)
      (f := fun p : β => paddedRowSumZ (2 ^ B p.1) p.2 (fun code => selectedSummand x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) code (row p.1)))
      (g := fun _ : β => (0 : ℤ)) hguard hsum (PolyT.const 0)

set_option maxHeartbeats 2000000 in
theorem selectedRow_eq (x C n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : 2 ^ masterExponent x C n ≤ pad.length) (γ : ℕ) (row : List (List Bool)) :
    let P := K x (selectedLevel x n)
    let M := (clusterNet (sOf n) ℤ N k y (3 * P ^ 8) (by positivity)).coarsen
      (selectedTopGrouping (sOf n) N k y P hP)
    selectedRow x n pad N k y γ row = cliqueQueryFast (M.toClique (topWidth (sOf n) P + 1)) (γ, row) := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let L := topWidth s P + 1
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let M := (clusterNet s ℤ N k y g (by positivity)).coarsen χ
  by_cases hg : γ < P
  · let a : Fin P := ⟨γ, hg⟩
    let θ : Fin P → Factoring.Aux.Bits L := fun δ => if δ = a then zeroBits L else Groundwork.Grid.bitsOf L (row.getD δ.val [])
    have hp : P ^ g ≤ pad.length := (master_pad_suffices x C n pad hpad).1
    have hpoint : (fun code => selectedSummand x n pad N k y γ code row) =
        cliqueSummand M a θ ((κ.coarsen χ).e a) := by
      funext code
      have h := selectedSummand_eq x n N k y hP pad hp a row code
      dsimp only at h
      with_reducible exact h
    have hrow := selected_explicit_row_as_padded_sum x C n N k y hP pad hpad a θ
    dsimp only at hrow
    unfold selectedCodingZ at hrow
    unfold selectedRow cliqueQueryFast
    dsimp only
    change (if γ < P then paddedRowSumZ (2 ^ topEntryBits s P) pad (fun code => selectedSummand x n pad N k y γ code row) else 0) =
      (if h : γ < P then M.toClique L ⟨γ, h⟩ (fun δ => if δ.val = γ then zeroBits L else Groundwork.Grid.bitsOf L (row.getD δ.val [])) else 0)
    rw [if_pos hg, dif_pos hg]
    rw [hpoint]
    have hθ : (fun δ : Fin P => if δ.val = γ then zeroBits L else Groundwork.Grid.bitsOf L (row.getD δ.val [])) = θ := by
      funext δ
      simp [θ, a, Fin.ext_iff]
    rw [hθ]
    with_reducible exact hrow
  · unfold selectedRow cliqueQueryFast
    dsimp only
    rw [if_neg hg, dif_neg hg]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedRow_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedRow_eq
