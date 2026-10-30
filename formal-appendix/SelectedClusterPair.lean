/-
SPDX-License-Identifier: Apache-2.0

Numeric cross-pair check for two selected cluster states.
-/
import CrossCheck
import GatePair
import SelectedGateCode
import SelectedFieldBounds

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedClusterPair (x n : ℕ) (pad : List Unit) (γ code u v : ℕ) : Bool :=
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let A := pA s
  let B := topEntryBits s P
  let color := gateCluster A g
  let state := fun i => selectedGateCode x n pad γ code (color i) i
  crossCheck (gateCount A)
    (fun i j => decide (i < j ∧ ((color i = u ∧ color j = v) ∨ (color i = v ∧ color j = u))))
    (fun i j => gatePair s B i j (state i) (state j))

set_option maxHeartbeats 2000000 in
theorem selectedClusterPair_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n γ code u v : α → ℕ} (hn : Unary n) (hγ : PolyT γ) (hc : PolyT code)
    (hu : PolyT u) (hv : PolyT v) :
    PolyT fun p : α × List Unit => selectedClusterPair x (n p.1) p.2 (γ p.1) (code p.1) (u p.1) (v p.1) := by
  let β := α × List Unit
  let A : β → ℕ := fun p => pA (sOf (n p.1))
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let B : β → ℕ := fun p => topEntryBits (sOf (n p.1)) (P p)
  have hn' : Unary fun p : β => n p.1 := Unary.comp hn PolyT.id'.fst
  have hA : Unary A := Unary.comp (selected_pA_unary hn) PolyT.id'.fst
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hB : Unary B := Unary.comp (selected_topEntryBits_unary x hx hn) PolyT.id'.fst
  have hnodes : Unary fun p : β => gateCount (A p) :=
    Unary.add hA (Unary.add hA (Unary.add hA (Unary.add hA (Unary.mul hA hA))))
  have hcolor : PolyT fun p : β × ℕ => gateCluster (A p.1) (g p.1) p.2 :=
    gateCluster_polyT (PolyT.comp hA.poly PolyT.id'.fst) (PolyT.comp hg.poly PolyT.id'.fst) PolyT.id'.snd
  have hbase : PolyT fun p : (β × ℕ) × ℕ =>
      selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2 :=
    selectedGateCode_polyT x hx hn hγ hc
  have hmap : PolyT fun p : β × ℕ => ((p.1, gateCluster (A p.1) (g p.1) p.2), p.2) :=
    PolyT.prodMk (PolyT.prodMk PolyT.id'.fst hcolor) PolyT.id'.snd
  have hstate : PolyT fun p : β × ℕ => selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (gateCluster (A p.1) (g p.1) p.2) p.2 :=
    PolyT.comp (f := fun p : (β × ℕ) × ℕ => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2)
      (g := fun p : β × ℕ => ((p.1, gateCluster (A p.1) (g p.1) p.2), p.2)) hbase hmap
  have hcolorI : PolyT fun p : (β × ℕ) × ℕ => gateCluster (A p.1.1) (g p.1.1) p.1.2 := PolyT.comp hcolor PolyT.id'.fst
  have hcolorJ : PolyT fun p : (β × ℕ) × ℕ => gateCluster (A p.1.1) (g p.1.1) p.2 :=
    PolyT.comp (f := fun p : β × ℕ => gateCluster (A p.1) (g p.1) p.2)
      (g := fun p : (β × ℕ) × ℕ => (p.1.1, p.2)) hcolor (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hu' : PolyT fun p : (β × ℕ) × ℕ => u p.1.1.1 := PolyT.comp hu PolyT.id'.fst.fst.fst
  have hv' : PolyT fun p : (β × ℕ) × ℕ => v p.1.1.1 := PolyT.comp hv PolyT.id'.fst.fst.fst
  have hused : PolyT fun p : (β × ℕ) × ℕ => decide (p.1.2 < p.2 ∧
      ((gateCluster (A p.1.1) (g p.1.1) p.1.2 = u p.1.1.1 ∧ gateCluster (A p.1.1) (g p.1.1) p.2 = v p.1.1.1) ∨
      (gateCluster (A p.1.1) (g p.1.1) p.1.2 = v p.1.1.1 ∧ gateCluster (A p.1.1) (g p.1.1) p.2 = u p.1.1.1))) :=
    PolyT.decide_and (PolyT.decide_lt PolyT.id'.fst.snd PolyT.id'.snd)
      (PolyT.decide_or (PolyT.decide_and (PolyT.decide_eq hcolorI hu') (PolyT.decide_eq hcolorJ hv'))
        (PolyT.decide_and (PolyT.decide_eq hcolorI hv') (PolyT.decide_eq hcolorJ hu')))
  have hci : PolyT fun p : (β × ℕ) × ℕ => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (gateCluster (A p.1.1) (g p.1.1) p.1.2) p.1.2 :=
    PolyT.comp hstate PolyT.id'.fst
  have hcj : PolyT fun p : (β × ℕ) × ℕ => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (gateCluster (A p.1.1) (g p.1.1) p.2) p.2 :=
    PolyT.comp (f := fun p : β × ℕ => selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (gateCluster (A p.1) (g p.1) p.2) p.2)
      (g := fun p : (β × ℕ) × ℕ => (p.1.1, p.2)) hstate (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hpair : PolyT fun p : (β × ℕ) × ℕ => gatePair (sOf (n p.1.1.1)) (B p.1.1) p.1.2 p.2
      (selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (gateCluster (A p.1.1) (g p.1.1) p.1.2) p.1.2)
      (selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (gateCluster (A p.1.1) (g p.1.1) p.2) p.2) :=
    gatePair_polyT (n := fun p : (β × ℕ) × ℕ => n p.1.1.1) (B := fun p : (β × ℕ) × ℕ => B p.1.1)
      (i := fun p : (β × ℕ) × ℕ => p.1.2) (j := fun p : (β × ℕ) × ℕ => p.2)
      (ci := fun p => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (gateCluster (A p.1.1) (g p.1.1) p.1.2) p.1.2)
      (cj := fun p => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (gateCluster (A p.1.1) (g p.1.1) p.2) p.2)
      (Unary.comp hn' PolyT.id'.fst.fst) (Unary.comp hB PolyT.id'.fst.fst) PolyT.id'.fst.snd PolyT.id'.snd hci hcj
  exact crossCheck_polyT (g := fun p : β => gateCount (A p))
    (used := fun p i j => decide (i < j ∧ ((gateCluster (A p) (g p) i = u p.1 ∧ gateCluster (A p) (g p) j = v p.1) ∨ (gateCluster (A p) (g p) i = v p.1 ∧ gateCluster (A p) (g p) j = u p.1))))
    (pair := fun p i j => gatePair (sOf (n p.1)) (B p) i j (selectedGateCode x (n p.1) p.2 (γ p.1) (code p.1) (gateCluster (A p) (g p) i) i)
      (selectedGateCode x (n p.1) p.2 (γ p.1) (code p.1) (gateCluster (A p) (g p) j) j)) hnodes hused hpair


set_option maxHeartbeats 2000000 in
theorem selectedClusterPair_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (u v : {a : Fin g // χ a = γ}), u.val ≠ v.val →
      let z := ((κ.coarsen χ).e γ).symm c
      selectedClusterPair x n pad γ.val c.val u.val.val v.val.val = true ↔
        (clusterNet (sOf n) ℤ N k y g (by positivity)).Ψ u.val v.val (z u) (z v) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let κ0 := mulCodingZ s N k y
  intro c u v huv
  let z := ((κ.coarsen χ).e γ).symm c
  let B := topEntryBits s P
  have hsize := selected_base_sizes_bound s N k y P hP
  have hcolor : ∀ a : Fin (mulNodes s), gateCluster (pA s) g a.val = (χ0 a).val := by
    intro a
    have h := gateCluster_eq s g (by positivity) ((MG.enum (pA s)).symm a)
    have hidx : (MG.enum (pA s) ((MG.enum (pA s)).symm a)).val = a.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) a)
    rw [hidx] at h
    exact h
  have hcode : ∀ (cl : {a : Fin g // χ a = γ}) (a : Fin (mulNodes s)) (ha : χ0 a = cl.val),
      selectedGateCode x n pad γ.val c.val cl.val.val a.val = ((κ0.e a) (z cl ⟨a, ha⟩)).val := by
    intro cl a ha
    exact selectedGateCode_eq x n N k y hP pad hpad γ c cl ⟨a, ha⟩
  have hu : ∀ a b : Fin (mulNodes s),
      decide (a.val < b.val ∧ ((gateCluster (pA s) g a.val = u.val.val ∧ gateCluster (pA s) g b.val = v.val.val) ∨
        (gateCluster (pA s) g a.val = v.val.val ∧ gateCluster (pA s) g b.val = u.val.val))) = decide (Net.UsedP χ0 u.val v.val (a,b)) := by
    intro a b
    rw [hcolor, hcolor]
    simp only [Fin.val_inj, Net.UsedP, Net.Cross]
    rfl
  have hp : ∀ a b : Fin (mulNodes s), Net.UsedP χ0 u.val v.val (a,b) →
      (if gatePair s B a.val b.val (selectedGateCode x n pad γ.val c.val (gateCluster (pA s) g a.val) a.val)
        (selectedGateCode x n pad γ.val c.val (gateCluster (pA s) g b.val) b.val) then (1 : ℤ) else 0) =
        ((mulNet s ℤ N k y).orient χ0 u.val v.val (a,b)).eval (z u) (z v) := by
    intro a b hab
    rw [hcolor, hcolor]
    rcases hab.2 with ⟨ha, hb⟩ | ⟨ha, hb⟩
    · rw [ha, hb, hcode u a ha, hcode v b hb]
      have h := gatePair_eq s B N k y a b (z u ⟨a, ha⟩) (z v ⟨b, hb⟩) hsize.1 hsize.2
      simpa only [Net.Ψ, Fac.eval, Net.orient, dif_pos ha, dif_pos hb] using h
    · have ha' : χ0 a ≠ u.val := ha ▸ Ne.symm huv
      have hb' : χ0 b ≠ v.val := hb ▸ huv
      rw [ha, hb, hcode v a ha, hcode u b hb]
      have h := gatePair_eq s B N k y a b (z v ⟨a, ha⟩) (z u ⟨b, hb⟩) hsize.1 hsize.2
      simp only [Net.Ψ, Fac.eval, Net.orient, dif_neg ha', dif_pos hb, dif_neg hb', dif_pos ha] at h ⊢
      exact h.trans (Finset.sum_congr rfl fun _ _ => mul_comm _ _)
  have hnodes : 0 < mulNodes s := by
    have hA := pA_pos s
    simp only [mulNodes]
    omega
  exact crossCheck_eq (mulNodes s) g hnodes χ0 u.val v.val (mulNet s ℤ N k y) (z u) (z v)
    (fun i j => decide (i < j ∧ ((gateCluster (pA s) g i = u.val.val ∧ gateCluster (pA s) g j = v.val.val) ∨
      (gateCluster (pA s) g i = v.val.val ∧ gateCluster (pA s) g j = u.val.val))))
    (fun i j => gatePair s B i j (selectedGateCode x n pad γ.val c.val (gateCluster (pA s) g i) i)
      (selectedGateCode x n pad γ.val c.val (gateCluster (pA s) g j) j)) hu hp

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedClusterPair_polyT

#print axioms IntegerFactoring.HalfAppendix.selectedClusterPair_eq
