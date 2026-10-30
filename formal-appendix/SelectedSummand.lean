/-
SPDX-License-Identifier: Apache-2.0

Executable selected clique-row summand.
-/
import SelectedStateCount
import SelectedTopWeight
import SelectedPadSide
import WordCode
import SelectedEntry
import CliqueTermCheck

open Groundwork Factoring.Aux Levels.Gen Finset

namespace IntegerFactoring.HalfAppendix

def selectedSummand (x n : ℕ) (pad : List Unit) (N k y γ code : ℕ) (row : List (List Bool)) : ℤ :=
  let P := K x (selectedLevel x n)
  if code < selectedStateCount x n pad γ then
    bif selectedTopWeight x n pad N k y γ code &&
      (List.range P).all (fun δ => decide (δ = γ) || selectedPadSide x n pad γ code δ (PolyT.bval (row.getD δ [])))
    then 1 else 0
  else 0

set_option maxHeartbeats 2000000 in
theorem selectedSummand_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n N k y γ code : α → ℕ} {row : α → List (List Bool)}
    (hn : Unary n) (hN : PolyT N) (hk : PolyT k) (hy : PolyT y)
    (hγ : PolyT γ) (hc : PolyT code) (hrow : PolyT row) :
    PolyT fun p : α × List Unit => selectedSummand x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) (row p.1) := by
  let β := α × List Unit
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hcount : PolyT fun p : β => selectedStateCount x (n p.1) p.2 (γ p.1) := selectedStateCount_polyT x hx (n := n) (γ := γ) hn hγ
  have hvalid : PolyT fun p : β => decide (code p.1 < selectedStateCount x (n p.1) p.2 (γ p.1)) :=
    PolyT.decide_lt (PolyT.comp hc PolyT.id'.fst) hcount
  have hweight : PolyT fun p : β => selectedTopWeight x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) :=
    selectedTopWeight_polyT x hx (n := n) (N := N) (k := k) (y := y) (γ := γ) (code := code) hn hN hk hy hγ hc
  have hword : PolyT fun p : α × ℕ => (row p.1).getD p.2 [] := PolyT.getD (PolyT.comp hrow PolyT.id'.fst) PolyT.id'.snd (PolyT.const [])
  have hq : PolyT fun p : α × ℕ => PolyT.bval ((row p.1).getD p.2 []) := PolyT.bval' hword
  have hbase : PolyT fun p : (α × ℕ) × List Unit => selectedPadSide x (n p.1.1) p.2
      (γ p.1.1) (code p.1.1) p.1.2 (PolyT.bval ((row p.1.1).getD p.1.2 [])) :=
    selectedPadSide_polyT x hx (n := fun p : α × ℕ => n p.1) (γ := fun p : α × ℕ => γ p.1)
      (code := fun p : α × ℕ => code p.1) (δ := fun p : α × ℕ => p.2)
      (q := fun p : α × ℕ => PolyT.bval ((row p.1).getD p.2 []))
      (Unary.comp hn PolyT.id'.fst) (PolyT.comp hγ PolyT.id'.fst) (PolyT.comp hc PolyT.id'.fst) PolyT.id'.snd hq
  have hside : PolyT fun p : β × ℕ => selectedPadSide x (n p.1.1) p.1.2
      (γ p.1.1) (code p.1.1) p.2 (PolyT.bval ((row p.1.1).getD p.2 [])) :=
    PolyT.comp (f := fun p : (α × ℕ) × List Unit => selectedPadSide x (n p.1.1) p.2
      (γ p.1.1) (code p.1.1) p.1.2 (PolyT.bval ((row p.1.1).getD p.1.2 [])))
      (g := fun p : β × ℕ => ((p.1.1, p.2), p.1.2)) hbase
      (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd) PolyT.id'.fst.snd)
  have htest : PolyT fun p : β × ℕ => decide (p.2 = γ p.1.1) ||
      selectedPadSide x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) p.2 (PolyT.bval ((row p.1.1).getD p.2 [])) :=
    PolyT.or' (PolyT.decide_eq PolyT.id'.snd (PolyT.comp hγ PolyT.id'.fst.fst)) hside
  have hall : PolyT fun p : β => (List.range (P p)).all (fun δ => decide (δ = γ p.1) ||
      selectedPadSide x (n p.1) p.2 (γ p.1) (code p.1) δ (PolyT.bval ((row p.1).getD δ []))) :=
    PolyT.all (f := fun p : β => fun δ => decide (δ = γ p.1) ||
      selectedPadSide x (n p.1) p.2 (γ p.1) (code p.1) δ (PolyT.bval ((row p.1).getD δ []))) htest (PolyT.range hP)
  have hinner : PolyT fun p : β => bif selectedTopWeight x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) &&
      (List.range (P p)).all (fun δ => decide (δ = γ p.1) || selectedPadSide x (n p.1) p.2 (γ p.1) (code p.1) δ (PolyT.bval ((row p.1).getD δ []))) then (1 : ℤ) else 0 :=
    PolyT.cond' (c := fun p : β => selectedTopWeight x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) &&
      (List.range (P p)).all (fun δ => decide (δ = γ p.1) || selectedPadSide x (n p.1) p.2 (γ p.1) (code p.1) δ (PolyT.bval ((row p.1).getD δ []))))
      (PolyT.and' hweight hall) (PolyT.const (1 : ℤ)) (PolyT.const 0)
  exact PolyT.ite' (c := fun p : β => code p.1 < selectedStateCount x (n p.1) p.2 (γ p.1))
    (f := fun p : β => bif selectedTopWeight x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) &&
      (List.range (P p)).all (fun δ => decide (δ = γ p.1) || selectedPadSide x (n p.1) p.2 (γ p.1) (code p.1) δ (PolyT.bval ((row p.1).getD δ []))) then (1 : ℤ) else 0)
    (g := fun _ : β => (0 : ℤ)) hvalid hinner (PolyT.const 0)


set_option maxHeartbeats 2000000 in
theorem selectedSummand_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) (row : List (List Bool)) (code : ℕ) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let L := topWidth (sOf n) P + 1
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    let M := (clusterNet (sOf n) ℤ N k y g (by positivity)).coarsen χ
    let θ : Fin P → Factoring.Aux.Bits L := fun δ => if δ = γ then zeroBits L else Groundwork.Grid.bitsOf L (row.getD δ.val [])
    selectedSummand x n pad N k y γ.val code row = cliqueSummand M γ θ ((κ.coarsen χ).e γ) code := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let L := topWidth s P + 1
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let M := (clusterNet s ℤ N k y g (by positivity)).coarsen χ
  let θ : Fin P → Factoring.Aux.Bits L := fun δ => if δ = γ then zeroBits L else Groundwork.Grid.bitsOf L (row.getD δ.val [])
  have hcount := selectedStateCount_eq x n N k y hP pad hpad γ
  change selectedStateCount x n pad γ.val = (κ.coarsen χ).s γ at hcount
  by_cases hv : code < (κ.coarsen χ).s γ
  · let c : Fin ((κ.coarsen χ).s γ) := ⟨code, hv⟩
    let z := ((κ.coarsen χ).e γ).symm c
    have hN : Net01 M := selectedTop_01 s N k y P hP
    have hweight := selectedTopWeight_eq x n N k y hP pad hpad γ c
    change selectedTopWeight x n pad N k y γ.val code = true ↔ M.W γ z = 1 at hweight
    have hleft : ∀ δ : Fin P, γ < δ →
        (selectedPadSide x n pad γ.val code δ.val (PolyT.bval (row.getD δ.val [])) = true ↔
          (M.padF L γ δ).l z (θ δ) = 1) := by
      intro δ hlt
      have hθ : θ δ = bdec L (PolyT.bval (row.getD δ.val [])) := by
        dsimp [θ]
        rw [if_neg hlt.ne', bitsOf_wordCode]
      rw [hθ]
      have hh := selectedPadSide_eq_left x n N k y hP pad hpad γ δ hlt.ne c (PolyT.bval (row.getD δ.val []))
      dsimp only at hh
      dsimp only [c] at hh
      with_reducible exact hh
    have hright : ∀ δ : Fin P, δ < γ →
        (selectedPadSide x n pad γ.val code δ.val (PolyT.bval (row.getD δ.val [])) = true ↔
          (M.padF L δ γ).r z (θ δ) = 1) := by
      intro δ hlt
      have hθ : θ δ = bdec L (PolyT.bval (row.getD δ.val [])) := by
        dsimp [θ]
        rw [if_neg hlt.ne, bitsOf_wordCode]
      rw [hθ]
      have hh := selectedPadSide_eq_right x n N k y hP pad hpad γ δ hlt.ne' c (PolyT.bval (row.getD δ.val []))
      dsimp only at hh
      dsimp only [c] at hh
      with_reducible exact hh
    let check := (List.range P).all (fun δ => decide (δ = γ.val) || selectedPadSide x n pad γ.val code δ (PolyT.bval (row.getD δ [])))
    have hall : check = true ↔
        (∀ δ : Fin P, γ < δ → (M.padF L γ δ).l z (θ δ) = 1) ∧
        (∀ δ : Fin P, δ < γ → (M.padF L δ γ).r z (θ δ) = 1) := by
      simp only [check, List.all_eq_true, List.mem_range, Bool.or_eq_true, decide_eq_true_eq]
      constructor
      · intro h
        constructor
        · intro δ hlt
          have hh := h δ.val δ.isLt
          rcases hh with he | he
          · exact False.elim (hlt.ne (Fin.ext he.symm))
          · exact (hleft δ hlt).mp he
        · intro δ hlt
          have hh := h δ.val δ.isLt
          rcases hh with he | he
          · exact False.elim (hlt.ne (Fin.ext he))
          · exact (hright δ hlt).mp he
      · rintro ⟨hl, hr⟩ i hi
        let δ : Fin P := ⟨i, hi⟩
        by_cases he : δ = γ
        · exact Or.inl (congrArg Fin.val he)
        · exact Or.inr (by
            rcases lt_or_gt_of_ne he with hlt | hgt
            · exact (hright δ hlt).mpr (hr δ hlt)
            · exact (hleft δ hgt).mpr (hl δ hgt))
    have hterm : cliqueTerm M γ θ z = 1 ↔ selectedTopWeight x n pad N k y γ.val code = true ∧ check = true := by
      rw [cliqueTerm_one_iff M hN γ θ z, ← hweight, ← hall]
    have hb := cliqueTerm_01 M hN γ θ z
    unfold selectedSummand cliqueSummand
    dsimp only
    rw [hcount, dif_pos hv, if_pos hv]
    change (bif selectedTopWeight x n pad N k y γ.val code && check then (1 : ℤ) else 0) = cliqueTerm M γ θ z
    cases hW : selectedTopWeight x n pad N k y γ.val code <;> cases hC : check <;>
      rcases hb with hb | hb <;> simp_all
  · unfold selectedSummand cliqueSummand
    dsimp only
    rw [hcount, dif_neg hv, if_neg hv]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedSummand_polyT

#print axioms IntegerFactoring.HalfAppendix.selectedSummand_eq
