/-
SPDX-License-Identifier: Apache-2.0

Numeric weight of a selected cluster state. The coarsened-factor and complete
row construction remain separate obligations.
-/
import WeightCheck
import GateWeight
import GatePair
import SelectedFieldBounds
import SelectedGateCode

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedClusterWeight (x n : ℕ) (pad : List Unit) (N k y γ code cluster : ℕ) : Bool :=
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let A := pA s
  let B := topEntryBits s P
  let color := gateCluster A g
  let state := fun i => selectedGateCode x n pad γ code cluster i
  weightCheck (gateCount A) (fun i => decide (color i = cluster))
    (fun i => gateWeight s B N k y i (state i))
    (fun i j => decide (i < j ∧ color i = cluster ∧ color j = cluster))
    (fun i j => gatePair s B i j (state i) (state j))

set_option maxHeartbeats 2000000 in
theorem selectedClusterWeight_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n N k y γ code cluster : α → ℕ} (hn : Unary n)
    (hN : PolyT N) (hk : PolyT k) (hy : PolyT y) (hγ : PolyT γ) (hc : PolyT code) (hcluster : PolyT cluster) :
    PolyT fun p : α × List Unit => selectedClusterWeight x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) (cluster p.1) := by
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
  have hcluster' : PolyT fun p : β => cluster p.1 := PolyT.comp hcluster PolyT.id'.fst
  have hbase : PolyT fun p : (β × ℕ) × ℕ =>
      selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2 :=
    selectedGateCode_polyT x hx hn hγ hc
  have hmap : PolyT fun p : β × ℕ => ((p.1, cluster p.1.1), p.2) :=
    PolyT.prodMk (PolyT.prodMk PolyT.id'.fst
      (PolyT.comp (f := fun p : β => cluster p.1) hcluster' PolyT.id'.fst)) PolyT.id'.snd
  have hstate : PolyT fun p : β × ℕ => selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (cluster p.1.1) p.2 :=
    PolyT.comp (f := fun p : (β × ℕ) × ℕ => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2)
      (g := fun p : β × ℕ => ((p.1, cluster p.1.1), p.2)) hbase hmap
  have hcolor : PolyT fun p : β × ℕ => gateCluster (A p.1) (g p.1) p.2 :=
    gateCluster_polyT (PolyT.comp hA.poly PolyT.id'.fst) (PolyT.comp hg.poly PolyT.id'.fst) PolyT.id'.snd
  have hmember : PolyT fun p : β × ℕ => decide (gateCluster (A p.1) (g p.1) p.2 = cluster p.1.1) :=
    PolyT.decide_eq hcolor (PolyT.comp (f := fun p : β => cluster p.1) hcluster' PolyT.id'.fst)
  have hnode : PolyT fun p : β × ℕ => gateWeight (sOf (n p.1.1)) (B p.1) (N p.1.1) (k p.1.1) (y p.1.1) p.2
      (selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (cluster p.1.1) p.2) :=
    gateWeight_polyT (n := fun p : β × ℕ => n p.1.1) (B := fun p : β × ℕ => B p.1)
      (N := fun p : β × ℕ => N p.1.1) (k := fun p : β × ℕ => k p.1.1) (y := fun p : β × ℕ => y p.1.1)
      (i := fun p : β × ℕ => p.2)
      (code := fun p : β × ℕ => selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (cluster p.1.1) p.2)
      (Unary.comp hn' PolyT.id'.fst) (Unary.comp hB PolyT.id'.fst)
      (PolyT.comp hN PolyT.id'.fst.fst) (PolyT.comp hk PolyT.id'.fst.fst) (PolyT.comp hy PolyT.id'.fst.fst) PolyT.id'.snd hstate
  have hcolorI : PolyT fun p : (β × ℕ) × ℕ => gateCluster (A p.1.1) (g p.1.1) p.1.2 := PolyT.comp hcolor PolyT.id'.fst
  have hcolorJ : PolyT fun p : (β × ℕ) × ℕ => gateCluster (A p.1.1) (g p.1.1) p.2 :=
    PolyT.comp (f := fun p : β × ℕ => gateCluster (A p.1) (g p.1) p.2)
      (g := fun p : (β × ℕ) × ℕ => (p.1.1, p.2)) hcolor (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hcl2 : PolyT fun p : (β × ℕ) × ℕ => cluster p.1.1.1 := PolyT.comp (f := fun p : β => cluster p.1) hcluster' PolyT.id'.fst.fst
  have hused : PolyT fun p : (β × ℕ) × ℕ => decide (p.1.2 < p.2 ∧ gateCluster (A p.1.1) (g p.1.1) p.1.2 = cluster p.1.1.1 ∧
      gateCluster (A p.1.1) (g p.1.1) p.2 = cluster p.1.1.1) :=
    PolyT.decide_and (PolyT.decide_lt PolyT.id'.fst.snd PolyT.id'.snd)
      (PolyT.decide_and (PolyT.decide_eq hcolorI hcl2) (PolyT.decide_eq hcolorJ hcl2))
  have hci : PolyT fun p : (β × ℕ) × ℕ => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (cluster p.1.1.1) p.1.2 :=
    PolyT.comp hstate PolyT.id'.fst
  have hcj : PolyT fun p : (β × ℕ) × ℕ => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (cluster p.1.1.1) p.2 :=
    PolyT.comp (f := fun p : β × ℕ => selectedGateCode x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) (cluster p.1.1) p.2)
      (g := fun p : (β × ℕ) × ℕ => (p.1.1, p.2)) hstate (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hpair : PolyT fun p : (β × ℕ) × ℕ => gatePair (sOf (n p.1.1.1)) (B p.1.1) p.1.2 p.2
      (selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (cluster p.1.1.1) p.1.2)
      (selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (cluster p.1.1.1) p.2) :=
    gatePair_polyT (n := fun p : (β × ℕ) × ℕ => n p.1.1.1) (B := fun p : (β × ℕ) × ℕ => B p.1.1)
      (i := fun p : (β × ℕ) × ℕ => p.1.2) (j := fun p : (β × ℕ) × ℕ => p.2)
      (ci := fun p => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (cluster p.1.1.1) p.1.2)
      (cj := fun p => selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) (cluster p.1.1.1) p.2)
      (Unary.comp hn' PolyT.id'.fst.fst) (Unary.comp hB PolyT.id'.fst.fst) PolyT.id'.fst.snd PolyT.id'.snd hci hcj
  exact weightCheck_polyT (g := fun p : β => gateCount (A p))
    (member := fun p i => decide (gateCluster (A p) (g p) i = cluster p.1))
    (node := fun p i => gateWeight (sOf (n p.1)) (B p) (N p.1) (k p.1) (y p.1) i (selectedGateCode x (n p.1) p.2 (γ p.1) (code p.1) (cluster p.1) i))
    (used := fun p i j => decide (i < j ∧ gateCluster (A p) (g p) i = cluster p.1 ∧ gateCluster (A p) (g p) j = cluster p.1))
    (pair := fun p i j => gatePair (sOf (n p.1)) (B p) i j (selectedGateCode x (n p.1) p.2 (γ p.1) (code p.1) (cluster p.1) i)
      (selectedGateCode x (n p.1) p.2 (γ p.1) (code p.1) (cluster p.1) j)) hnodes hmember hnode hused hpair


set_option maxHeartbeats 2000000 in
theorem selectedClusterWeight_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (cluster : {a : Fin g // χ a = γ}),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedClusterWeight x n pad N k y γ.val c.val cluster.val.val = true ↔
        (clusterNet (sOf n) ℤ N k y g (by positivity)).W cluster.val (z cluster) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let κ0 := mulCodingZ s N k y
  intro c cluster
  let z := ((κ.coarsen χ).e γ).symm c
  let B := topEntryBits s P
  have hsize := selected_base_sizes_bound s N k y P hP
  have hcolor : ∀ a : Fin (mulNodes s), gateCluster (pA s) g a.val = (χ0 a).val := by
    intro a
    have h := gateCluster_eq s g (by positivity) ((MG.enum (pA s)).symm a)
    have hidx : (MG.enum (pA s) ((MG.enum (pA s)).symm a)).val = a.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) a)
    rw [hidx] at h
    exact h
  have hcode : ∀ (a : Fin (mulNodes s)) (ha : χ0 a = cluster.val),
      selectedGateCode x n pad γ.val c.val cluster.val.val a.val = ((κ0.e a) (z cluster ⟨a, ha⟩)).val := by
    intro a ha
    exact selectedGateCode_eq x n N k y hP pad hpad γ c cluster ⟨a, ha⟩
  have hm : ∀ a : Fin (mulNodes s), decide (gateCluster (pA s) g a.val = cluster.val.val) = decide (χ0 a = cluster.val) := by
    intro a
    rw [hcolor]
    simp only [Fin.val_inj]
  have hn : ∀ (a : Fin (mulNodes s)) (ha : χ0 a = cluster.val),
      (if gateWeight s B N k y a.val (selectedGateCode x n pad γ.val c.val cluster.val.val a.val) then (1 : ℤ) else 0) =
        (mulNet s ℤ N k y).W a (z cluster ⟨a, ha⟩) := by
    intro a ha
    rw [hcode a ha]
    let gate := (MG.enum (pA s)).symm a
    let state := z cluster ⟨a, ha⟩
    have hw := gateWeight_eq s B N k y gate state (by intro r hr; exact hsize.1)
    have hw' : gateWeight s B N k y a.val ((κ0.e a) state).val =
        decide ((mulNet s ℤ N k y).W a state = 1) := by
      have hidx : (MG.enum (pA s) gate).val = a.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) a)
      have hci : ((κ0.e a) state).val = (gcode s gate state).val := rfl
      rw [hidx, ← hci] at hw
      exact hw
    rw [hw']
    have hbinary : (mulNet s ℤ N k y).W a state = 0 ∨ (mulNet s ℤ N k y).W a state = 1 := by
      exact mW_01 (R := ℤ) (s := s) N k y gate state
    change (if decide ((mulNet s ℤ N k y).W a state = 1) then (1 : ℤ) else 0) = (mulNet s ℤ N k y).W a state
    rcases hbinary with h | h <;> simp [h]
  have hu : ∀ a b : Fin (mulNodes s),
      decide (a.val < b.val ∧ gateCluster (pA s) g a.val = cluster.val.val ∧ gateCluster (pA s) g b.val = cluster.val.val) =
        decide (a < b ∧ χ0 a = cluster.val ∧ χ0 b = cluster.val) := by
    intro a b
    rw [hcolor, hcolor]
    simp only [Fin.val_inj, Fin.lt_def]
  have hp : ∀ (a b : Fin (mulNodes s)) (ha : χ0 a = cluster.val) (hb : χ0 b = cluster.val),
      (if gatePair s B a.val b.val (selectedGateCode x n pad γ.val c.val cluster.val.val a.val)
        (selectedGateCode x n pad γ.val c.val cluster.val.val b.val) then (1 : ℤ) else 0) =
        (mulNet s ℤ N k y).Ψ a b (z cluster ⟨a, ha⟩) (z cluster ⟨b, hb⟩) := by
    intro a b ha hb
    rw [hcode a ha, hcode b hb]
    exact gatePair_eq s B N k y a b (z cluster ⟨a, ha⟩) (z cluster ⟨b, hb⟩) hsize.1 hsize.2
  have hnodes : 0 < mulNodes s := by
    have hA := pA_pos s
    simp only [mulNodes]
    omega
  exact weightCheck_eq (mulNodes s) g hnodes χ0 cluster.val (mulNet s ℤ N k y) (z cluster)
    (fun i => decide (gateCluster (pA s) g i = cluster.val.val))
    (fun i => gateWeight s B N k y i (selectedGateCode x n pad γ.val c.val cluster.val.val i))
    (fun i j => decide (i < j ∧ gateCluster (pA s) g i = cluster.val.val ∧ gateCluster (pA s) g j = cluster.val.val))
    (fun i j => gatePair s B i j (selectedGateCode x n pad γ.val c.val cluster.val.val i) (selectedGateCode x n pad γ.val c.val cluster.val.val j)) hm hn hu hp

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedClusterWeight_polyT

#print axioms IntegerFactoring.HalfAppendix.selectedClusterWeight_eq
