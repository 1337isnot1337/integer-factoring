/-
SPDX-License-Identifier: Apache-2.0

Numeric left factor of a selected cluster against another cluster.
-/
import GateSide
import SelectedGateCode
import FactorCheck
import SelectedFieldBounds

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedClusterLeft (x n : ℕ) (pad : List Unit) (γ code u v q : ℕ) : Bool :=
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let A := pA s
  let B := topEntryBits s P
  let color := gateCluster A g
  let state := fun i => selectedGateCode x n pad γ code u i
  factorCheck (gateCount A) q (gateEdge A (fw s))
    (fun i j => decide (i < j ∧ ((color i = u ∧ color j = v) ∨ (color i = v ∧ color j = u))))
    (fun i j t => if color i = u then gateSide s B i j (state i) t
      else if color j = u then gateSide s B j i (state j) t else false)

set_option maxHeartbeats 2000000 in
theorem selectedClusterLeft_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n γ code u v q : α → ℕ} (hn : Unary n) (hγ : PolyT γ) (hc : PolyT code)
    (hu : PolyT u) (hv : PolyT v) (hq : PolyT q) :
    PolyT fun p : α × List Unit => selectedClusterLeft x (n p.1) p.2 (γ p.1) (code p.1) (u p.1) (v p.1) (q p.1) := by
  let β := α × List Unit
  let A : β → ℕ := fun p => pA (sOf (n p.1))
  let F : β → ℕ := fun p => fw (sOf (n p.1))
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let B : β → ℕ := fun p => topEntryBits (sOf (n p.1)) (P p)
  have hn' : Unary fun p : β => n p.1 := Unary.comp hn PolyT.id'.fst
  have hA : Unary A := Unary.comp (selected_pA_unary hn) PolyT.id'.fst
  have hF : Unary F := Unary.comp (selected_fw_unary hn) PolyT.id'.fst
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hB : Unary B := Unary.comp (selected_topEntryBits_unary x hx hn) PolyT.id'.fst
  have hnodes : Unary fun p : β => gateCount (A p) :=
    Unary.add hA (Unary.add hA (Unary.add hA (Unary.add hA (Unary.mul hA hA))))
  have hu' : PolyT fun p : β => u p.1 := PolyT.comp hu PolyT.id'.fst
  have hbase : PolyT fun p : (β × ℕ) × ℕ =>
      selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2 :=
    selectedGateCode_polyT x hx hn hγ hc
  have hmap : PolyT fun p : β × ℕ => ((p.1, u p.1.1), p.2) :=
    PolyT.prodMk (PolyT.prodMk PolyT.id'.fst (PolyT.comp (f := fun p : β => u p.1) hu' PolyT.id'.fst)) PolyT.id'.snd
  have hstate : PolyT fun p : β × ℕ => selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (u p.1.1) p.2 :=
    PolyT.comp (f := fun p : (β × ℕ) × ℕ => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2)
      (g := fun p : β × ℕ => ((p.1, u p.1.1), p.2)) hbase hmap
  have hcolor : PolyT fun p : β × ℕ => gateCluster (A p.1) (g p.1) p.2 :=
    gateCluster_polyT (PolyT.comp hA.poly PolyT.id'.fst) (PolyT.comp hg.poly PolyT.id'.fst) PolyT.id'.snd
  have hcolorI : PolyT fun p : (β × ℕ) × ℕ => gateCluster (A p.1.1) (g p.1.1) p.1.2 := PolyT.comp hcolor PolyT.id'.fst
  have hcolorJ : PolyT fun p : (β × ℕ) × ℕ => gateCluster (A p.1.1) (g p.1.1) p.2 :=
    PolyT.comp (f := fun p : β × ℕ => gateCluster (A p.1) (g p.1) p.2)
      (g := fun p : (β × ℕ) × ℕ => (p.1.1, p.2)) hcolor (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hu2 : PolyT fun p : (β × ℕ) × ℕ => u p.1.1.1 := PolyT.comp hu PolyT.id'.fst.fst.fst
  have hv2 : PolyT fun p : (β × ℕ) × ℕ => v p.1.1.1 := PolyT.comp hv PolyT.id'.fst.fst.fst
  have hused : PolyT fun p : (β × ℕ) × ℕ => decide (p.1.2 < p.2 ∧
      ((gateCluster (A p.1.1) (g p.1.1) p.1.2 = u p.1.1.1 ∧ gateCluster (A p.1.1) (g p.1.1) p.2 = v p.1.1.1) ∨
      (gateCluster (A p.1.1) (g p.1.1) p.1.2 = v p.1.1.1 ∧ gateCluster (A p.1.1) (g p.1.1) p.2 = u p.1.1.1))) :=
    PolyT.decide_and (PolyT.decide_lt PolyT.id'.fst.snd PolyT.id'.snd)
      (PolyT.decide_or (PolyT.decide_and (PolyT.decide_eq hcolorI hu2) (PolyT.decide_eq hcolorJ hv2))
        (PolyT.decide_and (PolyT.decide_eq hcolorI hv2) (PolyT.decide_eq hcolorJ hu2)))
  have hwidth : PolyT fun p : (β × ℕ) × ℕ => gateEdge (A p.1.1) (F p.1.1) p.1.2 p.2 :=
    gateEdge_polyT (PolyT.comp hA.poly PolyT.id'.fst.fst) (PolyT.comp hF.poly PolyT.id'.fst.fst) PolyT.id'.fst.snd PolyT.id'.snd
  have hstateI : PolyT fun p : ((β × ℕ) × ℕ) × ℕ => selectedGateCode x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) (u p.1.1.1.1) p.1.1.2 :=
    PolyT.comp hstate PolyT.id'.fst.fst
  have hstateJ : PolyT fun p : ((β × ℕ) × ℕ) × ℕ => selectedGateCode x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) (u p.1.1.1.1) p.1.2 :=
    PolyT.comp (f := fun p : β × ℕ => selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (u p.1.1) p.2)
      (g := fun p : ((β × ℕ) × ℕ) × ℕ => (p.1.1.1, p.1.2)) hstate (PolyT.prodMk PolyT.id'.fst.fst.fst PolyT.id'.fst.snd)
  have hsI := gateSide_polyT (n := fun p : ((β × ℕ) × ℕ) × ℕ => n p.1.1.1.1)
    (B := fun p => B p.1.1.1) (i := fun p => p.1.1.2) (j := fun p => p.1.2)
    (code := fun p => selectedGateCode x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) (u p.1.1.1.1) p.1.1.2)
    (q := fun p => p.2) (Unary.comp hn' PolyT.id'.fst.fst.fst) (Unary.comp hB PolyT.id'.fst.fst.fst)
    PolyT.id'.fst.fst.snd PolyT.id'.fst.snd hstateI PolyT.id'.snd
  have hsJ := gateSide_polyT (n := fun p : ((β × ℕ) × ℕ) × ℕ => n p.1.1.1.1)
    (B := fun p => B p.1.1.1) (i := fun p => p.1.2) (j := fun p => p.1.1.2)
    (code := fun p => selectedGateCode x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) (u p.1.1.1.1) p.1.2)
    (q := fun p => p.2) (Unary.comp hn' PolyT.id'.fst.fst.fst) (Unary.comp hB PolyT.id'.fst.fst.fst)
    PolyT.id'.fst.snd PolyT.id'.fst.fst.snd hstateJ PolyT.id'.snd
  have hside : PolyT fun p : ((β × ℕ) × ℕ) × ℕ =>
      if gateCluster (A p.1.1.1) (g p.1.1.1) p.1.1.2 = u p.1.1.1.1 then
        gateSide (sOf (n p.1.1.1.1)) (B p.1.1.1) p.1.1.2 p.1.2 (selectedGateCode x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) (u p.1.1.1.1) p.1.1.2) p.2
      else if gateCluster (A p.1.1.1) (g p.1.1.1) p.1.2 = u p.1.1.1.1 then
        gateSide (sOf (n p.1.1.1.1)) (B p.1.1.1) p.1.2 p.1.1.2 (selectedGateCode x (n p.1.1.1.1) p.1.1.1.2 (γ p.1.1.1.1) (code p.1.1.1.1) (u p.1.1.1.1) p.1.2) p.2 else false :=
    PolyT.ite' (PolyT.comp (PolyT.decide_eq hcolorI hu2) PolyT.id'.fst) hsI
      (PolyT.ite' (PolyT.comp (PolyT.decide_eq hcolorJ hu2) PolyT.id'.fst) hsJ (PolyT.const false))
  exact factorCheck_polyT (g := fun p : β => gateCount (A p)) (W := F) (q := fun p => q p.1)
    (width := fun p i j => gateEdge (A p) (F p) i j)
    (used := fun p i j => decide (i < j ∧ ((gateCluster (A p) (g p) i = u p.1 ∧ gateCluster (A p) (g p) j = v p.1) ∨ (gateCluster (A p) (g p) i = v p.1 ∧ gateCluster (A p) (g p) j = u p.1))))
    (side := fun p i j t => if gateCluster (A p) (g p) i = u p.1 then gateSide (sOf (n p.1)) (B p) i j (selectedGateCode x (n p.1) p.2 (γ p.1) (code p.1) (u p.1) i) t
      else if gateCluster (A p) (g p) j = u p.1 then gateSide (sOf (n p.1)) (B p) j i (selectedGateCode x (n p.1) p.2 (γ p.1) (code p.1) (u p.1) j) t else false)
    hnodes hF (PolyT.comp hq PolyT.id'.fst) hwidth hused hside (fun p i j => by unfold gateEdge; split <;> simp)


set_option maxHeartbeats 2000000 in
theorem selectedClusterLeft_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (u : {a : Fin g // χ a = γ}) (v : Fin g) (q : ℕ),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedClusterLeft x n pad γ.val c.val u.val.val v.val q = true ↔
        ((clusterNet (sOf n) ℤ N k y g (by positivity)).F u.val v).l (z u) (bdec _ q) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let κ0 := mulCodingZ s N k y
  intro c u v q
  let z := ((κ.coarsen χ).e γ).symm c
  let B := topEntryBits s P
  have hsize := selected_base_sizes_bound s N k y P hP
  have hcolor : ∀ a : Fin (mulNodes s), gateCluster (pA s) g a.val = (χ0 a).val := by
    intro a
    have h := gateCluster_eq s g (by positivity) ((MG.enum (pA s)).symm a)
    have hidx : (MG.enum (pA s) ((MG.enum (pA s)).symm a)).val = a.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) a)
    rw [hidx] at h
    exact h
  have hcode : ∀ (a : Fin (mulNodes s)) (ha : χ0 a = u.val),
      selectedGateCode x n pad γ.val c.val u.val.val a.val = ((κ0.e a) (z u ⟨a, ha⟩)).val := by
    intro a ha
    exact selectedGateCode_eq x n N k y hP pad hpad γ c u ⟨a, ha⟩
  have hw : ∀ a b : Fin (mulNodes s), gateEdge (pA s) (fw s) a.val b.val = (mulNet s ℤ N k y).c a b := by
    intro a b
    obtain ⟨i, rfl⟩ := (MG.enum (pA s)).surjective a
    obtain ⟨j, rfl⟩ := (MG.enum (pA s)).surjective b
    simpa [mulNet, GNet.toNet, mulG] using gateEdge_eq s i j
  have hu : ∀ a b : Fin (mulNodes s),
      decide (a.val < b.val ∧ ((gateCluster (pA s) g a.val = u.val.val ∧ gateCluster (pA s) g b.val = v.val) ∨
        (gateCluster (pA s) g a.val = v.val ∧ gateCluster (pA s) g b.val = u.val.val))) = decide (Net.UsedP χ0 u.val v (a,b)) := by
    intro a b
    rw [hcolor, hcolor]
    simp only [Fin.val_inj, Net.UsedP, Net.Cross]
    rfl
  have hs : ∀ (a b : Fin (mulNodes s)) (t : ℕ),
      (if (if gateCluster (pA s) g a.val = u.val.val then gateSide s B a.val b.val (selectedGateCode x n pad γ.val c.val u.val.val a.val) t
        else if gateCluster (pA s) g b.val = u.val.val then gateSide s B b.val a.val (selectedGateCode x n pad γ.val c.val u.val.val b.val) t else false) then (1 : ℤ) else 0) =
        ((mulNet s ℤ N k y).orient χ0 u.val v (a,b)).l (z u) (bdec _ t) := by
    intro a b t
    rw [hcolor, hcolor]
    simp only [Fin.val_inj]
    by_cases ha : χ0 a = u.val
    · rw [if_pos ha, hcode a ha]
      have h := gateSide_eq_left s B N k y a b (z u ⟨a, ha⟩) t hsize.1 hsize.2
      simpa [Net.orient, ha] using h
    · rw [if_neg ha]
      by_cases hb : χ0 b = u.val
      · rw [if_pos hb, hcode b hb]
        have h := gateSide_eq_right s B N k y a b (z u ⟨b, hb⟩) t hsize.1 hsize.2
        simpa [Net.orient, ha, hb] using h
      · simp [Net.orient, ha, hb]
  have hnodes : 0 < mulNodes s := by
    have hA := pA_pos s
    simp only [mulNodes]
    omega
  exact factorCheck_eq_left (mulNodes s) g hnodes χ0 u.val v (mulNet s ℤ N k y) (z u) q
    (gateEdge (pA s) (fw s))
    (fun i j => decide (i < j ∧ ((gateCluster (pA s) g i = u.val.val ∧ gateCluster (pA s) g j = v.val) ∨
      (gateCluster (pA s) g i = v.val ∧ gateCluster (pA s) g j = u.val.val))))
    (fun i j t => if gateCluster (pA s) g i = u.val.val then gateSide s B i j (selectedGateCode x n pad γ.val c.val u.val.val i) t
      else if gateCluster (pA s) g j = u.val.val then gateSide s B j i (selectedGateCode x n pad γ.val c.val u.val.val j) t else false) hw hu hs


set_option maxHeartbeats 2000000 in
theorem selectedClusterRight_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (u : {a : Fin g // χ a = γ}) (v : Fin g) (q : ℕ), u.val ≠ v →
      let z := ((κ.coarsen χ).e γ).symm c
      selectedClusterLeft x n pad γ.val c.val u.val.val v.val q = true ↔
        ((clusterNet (sOf n) ℤ N k y g (by positivity)).F v u.val).r (z u) (bdec _ q) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let κ0 := mulCodingZ s N k y
  intro c u v q hne
  let z := ((κ.coarsen χ).e γ).symm c
  let B := topEntryBits s P
  have hsize := selected_base_sizes_bound s N k y P hP
  have hcolor : ∀ a : Fin (mulNodes s), gateCluster (pA s) g a.val = (χ0 a).val := by
    intro a
    have h := gateCluster_eq s g (by positivity) ((MG.enum (pA s)).symm a)
    have hidx : (MG.enum (pA s) ((MG.enum (pA s)).symm a)).val = a.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) a)
    rw [hidx] at h
    exact h
  have hcode : ∀ (a : Fin (mulNodes s)) (ha : χ0 a = u.val),
      selectedGateCode x n pad γ.val c.val u.val.val a.val = ((κ0.e a) (z u ⟨a, ha⟩)).val := by
    intro a ha
    exact selectedGateCode_eq x n N k y hP pad hpad γ c u ⟨a, ha⟩
  have hw : ∀ a b : Fin (mulNodes s), gateEdge (pA s) (fw s) a.val b.val = (mulNet s ℤ N k y).c a b := by
    intro a b
    obtain ⟨i, rfl⟩ := (MG.enum (pA s)).surjective a
    obtain ⟨j, rfl⟩ := (MG.enum (pA s)).surjective b
    simpa [mulNet, GNet.toNet, mulG] using gateEdge_eq s i j
  have hu : ∀ a b : Fin (mulNodes s),
      decide (a.val < b.val ∧ ((gateCluster (pA s) g a.val = u.val.val ∧ gateCluster (pA s) g b.val = v.val) ∨
        (gateCluster (pA s) g a.val = v.val ∧ gateCluster (pA s) g b.val = u.val.val))) = decide (Net.UsedP χ0 v u.val (a,b)) := by
    intro a b
    rw [hcolor, hcolor]
    simp only [Fin.val_inj, Net.UsedP, Net.Cross]
    congr 1
    exact propext (and_congr Iff.rfl or_comm)
  have hs : ∀ (a b : Fin (mulNodes s)), Net.UsedP χ0 v u.val (a,b) → ∀ t : ℕ,
      (if (if gateCluster (pA s) g a.val = u.val.val then gateSide s B a.val b.val (selectedGateCode x n pad γ.val c.val u.val.val a.val) t
        else if gateCluster (pA s) g b.val = u.val.val then gateSide s B b.val a.val (selectedGateCode x n pad γ.val c.val u.val.val b.val) t else false) then (1 : ℤ) else 0) =
        ((mulNet s ℤ N k y).orient χ0 v u.val (a,b)).r (z u) (bdec _ t) := by
    intro a b hab t
    rw [hcolor, hcolor]
    simp only [Fin.val_inj]
    rcases hab.2 with ⟨ha, hb⟩ | ⟨ha, hb⟩
    · have ha' : χ0 a ≠ u.val := ha ▸ Ne.symm hne
      rw [if_neg ha', if_pos hb, hcode b hb]
      have h := gateSide_eq_right s B N k y a b (z u ⟨b, hb⟩) t hsize.1 hsize.2
      simpa [Net.orient, hb] using h
    · have hb' : χ0 b ≠ u.val := hb ▸ Ne.symm hne
      rw [if_pos ha, hcode a ha]
      have h := gateSide_eq_left s B N k y a b (z u ⟨a, ha⟩) t hsize.1 hsize.2
      simpa [Net.orient, ha, hb'] using h
  have hnodes : 0 < mulNodes s := by
    have hA := pA_pos s
    simp only [mulNodes]
    omega
  exact factorCheck_eq_right_used (mulNodes s) g hnodes χ0 v u.val (mulNet s ℤ N k y) (z u) q
    (gateEdge (pA s) (fw s))
    (fun i j => decide (i < j ∧ ((gateCluster (pA s) g i = u.val.val ∧ gateCluster (pA s) g j = v.val) ∨
      (gateCluster (pA s) g i = v.val ∧ gateCluster (pA s) g j = u.val.val))))
    (fun i j t => if gateCluster (pA s) g i = u.val.val then gateSide s B i j (selectedGateCode x n pad γ.val c.val u.val.val i) t
      else if gateCluster (pA s) g j = u.val.val then gateSide s B j i (selectedGateCode x n pad γ.val c.val u.val.val j) t else false) hw hu hs

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedClusterLeft_polyT

#print axioms IntegerFactoring.HalfAppendix.selectedClusterLeft_eq

#print axioms IntegerFactoring.HalfAppendix.selectedClusterRight_eq
