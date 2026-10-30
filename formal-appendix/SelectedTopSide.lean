/-
SPDX-License-Identifier: Apache-2.0

Numeric factor side of a selected top-group state.
-/
import SelectedClusterSide
import CrossWidthBound
import BinaryCoarsen

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedTopLeft (x n : ℕ) (pad : List Unit) (γ code δ q : ℕ) : Bool :=
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let color := selectedGroupDigit x n pad
  factorCheck g q (clusterEdge (pA s) (fw s) g)
    (fun i j => decide (i < j ∧ ((color i = γ ∧ color j = δ) ∨ (color i = δ ∧ color j = γ))))
    (fun i j t => if color i = γ then selectedClusterLeft x n pad γ code i j t
      else if color j = γ then selectedClusterLeft x n pad γ code j i t else false)

set_option maxHeartbeats 2000000 in
theorem selectedTopLeft_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n γ code δ q : α → ℕ} (hn : Unary n) (hγ : PolyT γ) (hc : PolyT code)
    (hδ : PolyT δ) (hq : PolyT q) :
    PolyT fun p : α × List Unit => selectedTopLeft x (n p.1) p.2 (γ p.1) (code p.1) (δ p.1) (q p.1) := by
  let β := α × List Unit
  let A : β → ℕ := fun p => pA (sOf (n p.1))
  let F : β → ℕ := fun p => fw (sOf (n p.1))
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let W : β → ℕ := fun p => gateCount (A p) * gateCount (A p) * F p
  have hA : Unary A := Unary.comp (selected_pA_unary hn) PolyT.id'.fst
  have hF : Unary F := Unary.comp (selected_fw_unary hn) PolyT.id'.fst
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hnodes : Unary fun p : β => gateCount (A p) :=
    Unary.add hA (Unary.add hA (Unary.add hA (Unary.add hA (Unary.mul hA hA))))
  have hW : Unary W := Unary.mul (Unary.mul hnodes hnodes) hF
  have hcolor : PolyT fun p : β × ℕ => selectedGroupDigit x (n p.1.1) p.1.2 p.2 :=
    selectedGroupDigit_polyT x hx hn
  have hcolorI : PolyT fun p : (β × ℕ) × ℕ => selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.1.2 := PolyT.comp hcolor PolyT.id'.fst
  have hcolorJ : PolyT fun p : (β × ℕ) × ℕ => selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.2 :=
    PolyT.comp (f := fun p : β × ℕ => selectedGroupDigit x (n p.1.1) p.1.2 p.2)
      (g := fun p : (β × ℕ) × ℕ => (p.1.1, p.2)) hcolor (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hγ2 : PolyT fun p : (β × ℕ) × ℕ => γ p.1.1.1 := PolyT.comp hγ PolyT.id'.fst.fst.fst
  have hδ2 : PolyT fun p : (β × ℕ) × ℕ => δ p.1.1.1 := PolyT.comp hδ PolyT.id'.fst.fst.fst
  have hused : PolyT fun p : (β × ℕ) × ℕ => decide (p.1.2 < p.2 ∧
      ((selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.1.2 = γ p.1.1.1 ∧ selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.2 = δ p.1.1.1) ∨
      (selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.1.2 = δ p.1.1.1 ∧ selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.2 = γ p.1.1.1))) :=
    PolyT.decide_and (PolyT.decide_lt PolyT.id'.fst.snd PolyT.id'.snd)
      (PolyT.decide_or (PolyT.decide_and (PolyT.decide_eq hcolorI hγ2) (PolyT.decide_eq hcolorJ hδ2))
        (PolyT.decide_and (PolyT.decide_eq hcolorI hδ2) (PolyT.decide_eq hcolorJ hγ2)))
  have hwidth : PolyT fun p : (β × ℕ) × ℕ => clusterEdge (A p.1.1) (F p.1.1) (g p.1.1) p.1.2 p.2 :=
    clusterEdge_polyT (Unary.comp hA PolyT.id'.fst.fst) (PolyT.comp hF.poly PolyT.id'.fst.fst)
      (PolyT.comp hg.poly PolyT.id'.fst.fst) PolyT.id'.fst.snd PolyT.id'.snd
  let θ := ((α × ℕ) × ℕ) × ℕ
  have hbase : PolyT fun p : θ × List Unit => selectedClusterLeft x (n p.1.1.1.1) p.2
      (γ p.1.1.1.1) (code p.1.1.1.1) p.1.1.1.2 p.1.1.2 p.1.2 := by
    exact selectedClusterLeft_polyT x hx
      (n := fun p : θ => n p.1.1.1) (γ := fun p : θ => γ p.1.1.1) (code := fun p : θ => code p.1.1.1)
      (u := fun p : θ => p.1.1.2) (v := fun p : θ => p.1.2) (q := fun p : θ => p.2)
      (Unary.comp hn PolyT.id'.fst.fst.fst) (PolyT.comp hγ PolyT.id'.fst.fst.fst)
      (PolyT.comp hc PolyT.id'.fst.fst.fst) PolyT.id'.fst.fst.snd PolyT.id'.fst.snd PolyT.id'.snd
  have hfirst : PolyT fun p : ((β × ℕ) × ℕ) × ℕ => selectedClusterLeft x (n p.1.1.1.1) p.1.1.1.2
      (γ p.1.1.1.1) (code p.1.1.1.1) p.1.1.2 p.1.2 p.2 :=
    PolyT.comp (f := fun p : θ × List Unit => selectedClusterLeft x (n p.1.1.1.1) p.2
      (γ p.1.1.1.1) (code p.1.1.1.1) p.1.1.1.2 p.1.1.2 p.1.2)
      (g := fun p : ((β × ℕ) × ℕ) × ℕ => ((((p.1.1.1.1, p.1.1.2), p.1.2), p.2), p.1.1.1.2)) hbase
      (PolyT.prodMk (PolyT.prodMk (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst.fst.fst PolyT.id'.fst.fst.snd) PolyT.id'.fst.snd) PolyT.id'.snd) PolyT.id'.fst.fst.fst.snd)
  have hsecond : PolyT fun p : ((β × ℕ) × ℕ) × ℕ => selectedClusterLeft x (n p.1.1.1.1) p.1.1.1.2
      (γ p.1.1.1.1) (code p.1.1.1.1) p.1.2 p.1.1.2 p.2 :=
    PolyT.comp (f := fun p : ((β × ℕ) × ℕ) × ℕ => selectedClusterLeft x (n p.1.1.1.1) p.1.1.1.2
      (γ p.1.1.1.1) (code p.1.1.1.1) p.1.1.2 p.1.2 p.2)
      (g := fun p : ((β × ℕ) × ℕ) × ℕ => (((p.1.1.1, p.1.2), p.1.1.2), p.2)) hfirst
      (PolyT.prodMk (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst.fst PolyT.id'.fst.snd) PolyT.id'.fst.fst.snd) PolyT.id'.snd)
  have hside : PolyT fun p : ((β × ℕ) × ℕ) × ℕ =>
      if selectedGroupDigit x (n p.1.1.1.1) p.1.1.1.2 p.1.1.2 = γ p.1.1.1.1 then
        selectedClusterLeft x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) p.1.1.2 p.1.2 p.2
      else if selectedGroupDigit x (n p.1.1.1.1) p.1.1.1.2 p.1.2 = γ p.1.1.1.1 then
        selectedClusterLeft x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) p.1.2 p.1.1.2 p.2 else false :=
    PolyT.ite' (PolyT.comp (PolyT.decide_eq hcolorI hγ2) PolyT.id'.fst) hfirst
      (PolyT.ite' (PolyT.comp (PolyT.decide_eq hcolorJ hγ2) PolyT.id'.fst) hsecond (PolyT.const false))
  exact factorCheck_polyT (g := g) (W := W) (q := fun p => q p.1)
    (width := fun p i j => clusterEdge (A p) (F p) (g p) i j)
    (used := fun p i j => decide (i < j ∧ ((selectedGroupDigit x (n p.1) p.2 i = γ p.1 ∧ selectedGroupDigit x (n p.1) p.2 j = δ p.1) ∨
      (selectedGroupDigit x (n p.1) p.2 i = δ p.1 ∧ selectedGroupDigit x (n p.1) p.2 j = γ p.1))))
    (side := fun p i j t => if selectedGroupDigit x (n p.1) p.2 i = γ p.1 then selectedClusterLeft x (n p.1) p.2 (γ p.1) (code p.1) i j t
      else if selectedGroupDigit x (n p.1) p.2 j = γ p.1 then selectedClusterLeft x (n p.1) p.2 (γ p.1) (code p.1) j i t else false)
    hg hW (PolyT.comp hq PolyT.id'.fst) hwidth hused hside (fun p i j => clusterEdge_le _ _ _ _ _)


private theorem bool_value (b : Bool) (r : ℤ) (hb : Is01 r) (he : b = true ↔ r = 1) :
    (if b then (1 : ℤ) else 0) = r := by
  cases b <;> rcases hb with hb | hb <;> simp_all

set_option maxHeartbeats 2000000 in
theorem selectedTopLeft_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ δ : Fin (K x (selectedLevel x n))) (hne : γ ≠ δ) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (q : ℕ),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedTopLeft x n pad γ.val c.val δ.val q = true ↔
        (((clusterNet (sOf n) ℤ N k y g (by positivity)).coarsen χ).F γ δ).l z (bdec _ q) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let M := clusterNet s ℤ N k y g (by positivity)
  intro c q
  let z := ((κ.coarsen χ).e γ).symm c
  have h01 : Net01 M := coarsen_01 (mulNet s ℤ N k y) (multiplicationNet_01 s N k y) χ0
  have hw : ∀ a b : Fin g, clusterEdge (pA s) (fw s) g a.val b.val = M.c a b := by
    intro a b
    exact clusterEdge_eq s g N k y (by positivity) a b
  have hcolor : ∀ a : Fin g, selectedGroupDigit x n pad a.val = (χ a).val := by
    intro a
    exact selectedGroupDigit_eq x n N k y hP pad hpad a
  have hu : ∀ a b : Fin g,
      decide (a.val < b.val ∧ ((selectedGroupDigit x n pad a.val = γ.val ∧ selectedGroupDigit x n pad b.val = δ.val) ∨
        (selectedGroupDigit x n pad a.val = δ.val ∧ selectedGroupDigit x n pad b.val = γ.val))) = decide (Net.UsedP χ γ δ (a,b)) := by
    intro a b
    rw [hcolor, hcolor]
    simp only [Fin.val_inj, Net.UsedP, Net.Cross]
    rfl
  have hs : ∀ (a b : Fin g), Net.UsedP χ γ δ (a,b) → ∀ t : ℕ,
      (if (if selectedGroupDigit x n pad a.val = γ.val then selectedClusterLeft x n pad γ.val c.val a.val b.val t
        else if selectedGroupDigit x n pad b.val = γ.val then selectedClusterLeft x n pad γ.val c.val b.val a.val t else false) then (1 : ℤ) else 0) =
        (M.orient χ γ δ (a,b)).l z (bdec _ t) := by
    intro a b hab t
    rw [hcolor, hcolor]
    simp only [Fin.val_inj]
    rcases hab.2 with ⟨ha, hb⟩ | ⟨ha, hb⟩
    · rw [if_pos ha]
      have he := selectedClusterLeft_eq x n N k y hP pad hpad γ c ⟨a, ha⟩ b t
      change selectedClusterLeft x n pad γ.val c.val a.val b.val t = true ↔ (M.F a b).l (z ⟨a, ha⟩) (bdec _ t) = 1 at he
      have h := bool_value _ _ ((h01.factors a b).1 (z ⟨a, ha⟩) (bdec _ t)) he
      simpa [Net.orient, ha] using h
    · have ha' : χ a ≠ γ := ha ▸ Ne.symm hne
      have hba : b ≠ a := by
        intro he
        subst b
        exact hne (hb.symm.trans ha)
      rw [if_neg ha', if_pos hb]
      have he := selectedClusterRight_eq x n N k y hP pad hpad γ c ⟨b, hb⟩ a t hba
      change selectedClusterLeft x n pad γ.val c.val b.val a.val t = true ↔ (M.F a b).r (z ⟨b, hb⟩) (bdec _ t) = 1 at he
      have h := bool_value _ _ ((h01.factors a b).2 (z ⟨b, hb⟩) (bdec _ t)) he
      simpa [Net.orient, ha', hb] using h
  exact factorCheck_eq_left_used g P (by positivity) χ γ δ M z q
    (clusterEdge (pA s) (fw s) g)
    (fun i j => decide (i < j ∧ ((selectedGroupDigit x n pad i = γ.val ∧ selectedGroupDigit x n pad j = δ.val) ∨
      (selectedGroupDigit x n pad i = δ.val ∧ selectedGroupDigit x n pad j = γ.val))))
    (fun i j t => if selectedGroupDigit x n pad i = γ.val then selectedClusterLeft x n pad γ.val c.val i j t
      else if selectedGroupDigit x n pad j = γ.val then selectedClusterLeft x n pad γ.val c.val j i t else false) hw hu hs


set_option maxHeartbeats 2000000 in
theorem selectedTopRight_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ δ : Fin (K x (selectedLevel x n))) (hne : γ ≠ δ) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (q : ℕ),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedTopLeft x n pad γ.val c.val δ.val q = true ↔
        (((clusterNet (sOf n) ℤ N k y g (by positivity)).coarsen χ).F δ γ).r z (bdec _ q) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let M := clusterNet s ℤ N k y g (by positivity)
  intro c q
  let z := ((κ.coarsen χ).e γ).symm c
  have h01 : Net01 M := coarsen_01 (mulNet s ℤ N k y) (multiplicationNet_01 s N k y) χ0
  have hw : ∀ a b : Fin g, clusterEdge (pA s) (fw s) g a.val b.val = M.c a b := by
    intro a b
    exact clusterEdge_eq s g N k y (by positivity) a b
  have hcolor : ∀ a : Fin g, selectedGroupDigit x n pad a.val = (χ a).val := by
    intro a
    exact selectedGroupDigit_eq x n N k y hP pad hpad a
  have hu : ∀ a b : Fin g,
      decide (a.val < b.val ∧ ((selectedGroupDigit x n pad a.val = γ.val ∧ selectedGroupDigit x n pad b.val = δ.val) ∨
        (selectedGroupDigit x n pad a.val = δ.val ∧ selectedGroupDigit x n pad b.val = γ.val))) = decide (Net.UsedP χ δ γ (a,b)) := by
    intro a b
    rw [hcolor, hcolor]
    simp only [Fin.val_inj, Net.UsedP, Net.Cross]
    congr 1
    exact propext (and_congr Iff.rfl or_comm)
  have hs : ∀ (a b : Fin g), Net.UsedP χ δ γ (a,b) → ∀ t : ℕ,
      (if (if selectedGroupDigit x n pad a.val = γ.val then selectedClusterLeft x n pad γ.val c.val a.val b.val t
        else if selectedGroupDigit x n pad b.val = γ.val then selectedClusterLeft x n pad γ.val c.val b.val a.val t else false) then (1 : ℤ) else 0) =
        (M.orient χ δ γ (a,b)).r z (bdec _ t) := by
    intro a b hab t
    rw [hcolor, hcolor]
    simp only [Fin.val_inj]
    rcases hab.2 with ⟨ha, hb⟩ | ⟨ha, hb⟩
    · have ha' : χ a ≠ γ := ha ▸ Ne.symm hne
      have hba : b ≠ a := by
        intro he
        subst b
        exact hne (hb.symm.trans ha)
      rw [if_neg ha', if_pos hb]
      have he := selectedClusterRight_eq x n N k y hP pad hpad γ c ⟨b, hb⟩ a t hba
      change selectedClusterLeft x n pad γ.val c.val b.val a.val t = true ↔ (M.F a b).r (z ⟨b, hb⟩) (bdec _ t) = 1 at he
      have h := bool_value _ _ ((h01.factors a b).2 (z ⟨b, hb⟩) (bdec _ t)) he
      simpa [Net.orient, hb] using h
    · have hb' : χ b ≠ γ := hb ▸ Ne.symm hne
      rw [if_pos ha]
      have he := selectedClusterLeft_eq x n N k y hP pad hpad γ c ⟨a, ha⟩ b t
      change selectedClusterLeft x n pad γ.val c.val a.val b.val t = true ↔ (M.F a b).l (z ⟨a, ha⟩) (bdec _ t) = 1 at he
      have h := bool_value _ _ ((h01.factors a b).1 (z ⟨a, ha⟩) (bdec _ t)) he
      simpa [Net.orient, ha, hb'] using h
  exact factorCheck_eq_right_used g P (by positivity) χ δ γ M z q
    (clusterEdge (pA s) (fw s) g)
    (fun i j => decide (i < j ∧ ((selectedGroupDigit x n pad i = γ.val ∧ selectedGroupDigit x n pad j = δ.val) ∨
      (selectedGroupDigit x n pad i = δ.val ∧ selectedGroupDigit x n pad j = γ.val))))
    (fun i j t => if selectedGroupDigit x n pad i = γ.val then selectedClusterLeft x n pad γ.val c.val i j t
      else if selectedGroupDigit x n pad j = γ.val then selectedClusterLeft x n pad γ.val c.val j i t else false) hw hu hs

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedTopLeft_polyT

#print axioms IntegerFactoring.HalfAppendix.selectedTopLeft_eq

#print axioms IntegerFactoring.HalfAppendix.selectedTopRight_eq
