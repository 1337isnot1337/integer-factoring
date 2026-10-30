/-
SPDX-License-Identifier: Apache-2.0

Executable weight of a selected top-group state.
-/
import SelectedClusterWeight
import SelectedClusterPair
import BinaryCoarsen

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedTopWeight (x n : ℕ) (pad : List Unit) (N k y γ code : ℕ) : Bool :=
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  let color := selectedGroupDigit x n pad
  weightCheck g (fun i => decide (color i = γ))
    (fun i => selectedClusterWeight x n pad N k y γ code i)
    (fun i j => decide (i < j ∧ color i = γ ∧ color j = γ))
    (fun i j => selectedClusterPair x n pad γ code i j)

set_option maxHeartbeats 2000000 in
theorem selectedTopWeight_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n N k y γ code : α → ℕ} (hn : Unary n) (hN : PolyT N) (hk : PolyT k)
    (hy : PolyT y) (hγ : PolyT γ) (hc : PolyT code) :
    PolyT fun p : α × List Unit => selectedTopWeight x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) := by
  let β := α × List Unit
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hcolor : PolyT fun p : β × ℕ => selectedGroupDigit x (n p.1.1) p.1.2 p.2 :=
    selectedGroupDigit_polyT x hx hn
  have hγ' : PolyT fun p : (β × ℕ) × ℕ => γ p.1.1.1 := PolyT.comp hγ PolyT.id'.fst.fst.fst
  have hmember : PolyT fun p : β × ℕ => decide (selectedGroupDigit x (n p.1.1) p.1.2 p.2 = γ p.1.1) :=
    PolyT.decide_eq hcolor (PolyT.comp hγ PolyT.id'.fst.fst)
  have hbase : PolyT fun p : (α × ℕ) × List Unit => selectedClusterWeight x (n p.1.1) p.2
      (N p.1.1) (k p.1.1) (y p.1.1) (γ p.1.1) (code p.1.1) p.1.2 :=
    selectedClusterWeight_polyT x hx (n := fun p : α × ℕ => n p.1) (N := fun p : α × ℕ => N p.1)
      (k := fun p : α × ℕ => k p.1) (y := fun p : α × ℕ => y p.1)
      (γ := fun p : α × ℕ => γ p.1) (code := fun p : α × ℕ => code p.1) (cluster := fun p : α × ℕ => p.2)
      (Unary.comp hn PolyT.id'.fst) (PolyT.comp hN PolyT.id'.fst)
      (PolyT.comp hk PolyT.id'.fst) (PolyT.comp hy PolyT.id'.fst) (PolyT.comp hγ PolyT.id'.fst)
      (PolyT.comp hc PolyT.id'.fst) PolyT.id'.snd
  have hnode : PolyT fun p : β × ℕ => selectedClusterWeight x (n p.1.1) p.1.2
      (N p.1.1) (k p.1.1) (y p.1.1) (γ p.1.1) (code p.1.1) p.2 :=
    PolyT.comp (f := fun p : (α × ℕ) × List Unit => selectedClusterWeight x (n p.1.1) p.2
      (N p.1.1) (k p.1.1) (y p.1.1) (γ p.1.1) (code p.1.1) p.1.2)
      (g := fun p : β × ℕ => ((p.1.1, p.2), p.1.2)) hbase
      (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd) PolyT.id'.fst.snd)
  have hcolorI : PolyT fun p : (β × ℕ) × ℕ => selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.1.2 :=
    PolyT.comp hcolor PolyT.id'.fst
  have hcolorJ : PolyT fun p : (β × ℕ) × ℕ => selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.2 :=
    PolyT.comp (f := fun p : β × ℕ => selectedGroupDigit x (n p.1.1) p.1.2 p.2)
      (g := fun p : (β × ℕ) × ℕ => (p.1.1, p.2)) hcolor (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hused : PolyT fun p : (β × ℕ) × ℕ => decide (p.1.2 < p.2 ∧
      selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.1.2 = γ p.1.1.1 ∧
      selectedGroupDigit x (n p.1.1.1) p.1.1.2 p.2 = γ p.1.1.1) :=
    PolyT.decide_and (PolyT.decide_lt PolyT.id'.fst.snd PolyT.id'.snd)
      (PolyT.decide_and (PolyT.decide_eq hcolorI hγ') (PolyT.decide_eq hcolorJ hγ'))
  have hpairbase : PolyT fun p : ((α × ℕ) × ℕ) × List Unit => selectedClusterPair x (n p.1.1.1) p.2
      (γ p.1.1.1) (code p.1.1.1) p.1.1.2 p.1.2 :=
    selectedClusterPair_polyT x hx (n := fun p : (α × ℕ) × ℕ => n p.1.1)
      (γ := fun p : (α × ℕ) × ℕ => γ p.1.1) (code := fun p : (α × ℕ) × ℕ => code p.1.1)
      (u := fun p : (α × ℕ) × ℕ => p.1.2) (v := fun p : (α × ℕ) × ℕ => p.2)
      (Unary.comp hn PolyT.id'.fst.fst)
      (PolyT.comp hγ PolyT.id'.fst.fst) (PolyT.comp hc PolyT.id'.fst.fst) PolyT.id'.fst.snd PolyT.id'.snd
  have hpair : PolyT fun p : (β × ℕ) × ℕ => selectedClusterPair x (n p.1.1.1) p.1.1.2
      (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2 :=
    PolyT.comp (f := fun p : ((α × ℕ) × ℕ) × List Unit => selectedClusterPair x (n p.1.1.1) p.2
      (γ p.1.1.1) (code p.1.1.1) p.1.1.2 p.1.2)
      (g := fun p : (β × ℕ) × ℕ => (((p.1.1.1, p.1.2), p.2), p.1.1.2)) hpairbase
      (PolyT.prodMk (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst.fst PolyT.id'.fst.snd) PolyT.id'.snd) PolyT.id'.fst.fst.snd)
  exact weightCheck_polyT (g := g)
    (member := fun p i => decide (selectedGroupDigit x (n p.1) p.2 i = γ p.1))
    (node := fun p i => selectedClusterWeight x (n p.1) p.2 (N p.1) (k p.1) (y p.1) (γ p.1) (code p.1) i)
    (used := fun p i j => decide (i < j ∧ selectedGroupDigit x (n p.1) p.2 i = γ p.1 ∧ selectedGroupDigit x (n p.1) p.2 j = γ p.1))
    (pair := fun p i j => selectedClusterPair x (n p.1) p.2 (γ p.1) (code p.1) i j) hg hmember hnode hused hpair


set_option maxHeartbeats 2000000 in
theorem selectedTopWeight_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ c : Fin ((κ.coarsen χ).s γ),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedTopWeight x n pad N k y γ.val c.val = true ↔
        ((clusterNet (sOf n) ℤ N k y g (by positivity)).coarsen χ).W γ z = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let M := clusterNet s ℤ N k y g (by positivity)
  intro c
  let z := ((κ.coarsen χ).e γ).symm c
  have h01 : Net01 M := coarsen_01 (mulNet s ℤ N k y) (multiplicationNet_01 s N k y) χ0
  have hm : ∀ a : Fin g, decide (selectedGroupDigit x n pad a.val = γ.val) = decide (χ a = γ) := by
    intro a
    rw [selectedGroupDigit_eq x n N k y hP pad hpad a]
    simp only [Fin.val_inj]
    rfl
  have hn : ∀ (a : Fin g) (ha : χ a = γ),
      (if selectedClusterWeight x n pad N k y γ.val c.val a.val then (1 : ℤ) else 0) = M.W a (z ⟨a, ha⟩) := by
    intro a ha
    have he := selectedClusterWeight_eq x n N k y hP pad hpad γ c ⟨a, ha⟩
    change selectedClusterWeight x n pad N k y γ.val c.val a.val = true ↔ M.W a (z ⟨a, ha⟩) = 1 at he
    have hb := h01.weights a (z ⟨a, ha⟩)
    cases hw : selectedClusterWeight x n pad N k y γ.val c.val a.val <;> rcases hb with hb | hb <;> simp_all
  have hu : ∀ a b : Fin g, decide (a.val < b.val ∧ selectedGroupDigit x n pad a.val = γ.val ∧ selectedGroupDigit x n pad b.val = γ.val) =
      decide (a < b ∧ χ a = γ ∧ χ b = γ) := by
    intro a b
    rw [selectedGroupDigit_eq x n N k y hP pad hpad a, selectedGroupDigit_eq x n N k y hP pad hpad b]
    simp only [Fin.val_inj]
    rfl
  have hp : ∀ (a b : Fin g) (hab : a < b) (ha : χ a = γ) (hb : χ b = γ),
      (if selectedClusterPair x n pad γ.val c.val a.val b.val then (1 : ℤ) else 0) =
        M.Ψ a b (z ⟨a, ha⟩) (z ⟨b, hb⟩) := by
    intro a b hab ha hb
    have he := selectedClusterPair_eq x n N k y hP pad hpad γ c ⟨a, ha⟩ ⟨b, hb⟩ hab.ne
    change selectedClusterPair x n pad γ.val c.val a.val b.val = true ↔ M.Ψ a b (z ⟨a, ha⟩) (z ⟨b, hb⟩) = 1 at he
    have hbin := h01.interactions a b hab (z ⟨a, ha⟩) (z ⟨b, hb⟩)
    cases hw : selectedClusterPair x n pad γ.val c.val a.val b.val <;> rcases hbin with hbin | hbin <;> simp_all
  exact weightCheck_eq g P (by positivity) χ γ M z
    (fun i => decide (selectedGroupDigit x n pad i = γ.val))
    (fun i => selectedClusterWeight x n pad N k y γ.val c.val i)
    (fun i j => decide (i < j ∧ selectedGroupDigit x n pad i = γ.val ∧ selectedGroupDigit x n pad j = γ.val))
    (fun i j => selectedClusterPair x n pad γ.val c.val i j) hm hn hu hp

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedTopWeight_polyT

#print axioms IntegerFactoring.HalfAppendix.selectedTopWeight_eq
