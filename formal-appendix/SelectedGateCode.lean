/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time extraction of a multiplication-gate code from a selected
mixed-radix top-group code. The gate's internal fields and summands remain open.
-/
import SelectedStateDigit
import ClusterStateDigit

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedGateCode (x n : ℕ) (pad : List Unit) (γ code cluster gate : ℕ) : ℕ :=
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  let B := topEntryBits (sOf n) P
  clusterStateDigit (sOf n) g B cluster (selectedStateDigit x n pad γ code cluster) gate

theorem selectedGateCode_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n γ code : α → ℕ}
    (hn : Unary n) (hγ : PolyT γ) (hc : PolyT code) :
    PolyT fun p : ((α × List Unit) × ℕ) × ℕ =>
      selectedGateCode x (n p.1.1.1) p.1.1.2 (γ p.1.1.1) (code p.1.1.1) p.1.2 p.2 := by
  let β := (α × List Unit) × ℕ
  let n' : β → ℕ := fun p => n p.1.1
  let P : β → ℕ := fun p => K x (selectedLevel x (n' p))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let B : β → ℕ := fun p => topEntryBits (sOf (n' p)) (P p)
  let c : β → ℕ := fun p => selectedStateDigit x (n' p) p.1.2 (γ p.1.1) (code p.1.1) p.2
  have hn' : Unary n' := Unary.comp hn PolyT.id'.fst.fst
  have hP : Unary P := selected_K_unary x hx hn'
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hB : Unary B := selected_topEntryBits_unary x hx hn'
  have hc' : PolyT c := selectedStateDigit_polyT x hx hn hγ hc
  exact clusterStateDigit_polyT (n := n') (g := g) (B := B) (γ := fun p : β => p.2) (code := c)
    hn' hg.poly hB PolyT.id'.snd hc'

theorem selectedGateCode_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    let χ0 := fun a => clMG (sOf n) g (by positivity) ((MG.enum (pA (sOf n))).symm a)
    let κ0 := mulCodingZ (sOf n) N k y
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (cluster : {a : Fin g // χ a = γ})
      (gate : {a : Fin (mulNodes (sOf n)) // χ0 a = cluster.val}),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedGateCode x n pad γ.val c.val cluster.val.val gate.val.val =
        ((κ0.e gate.val) (z cluster gate)).val := by
  dsimp only
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping (sOf n) N k y P hP
  let κ := clusterCodingZ (sOf n) g N k y (by positivity)
  let χ0 := fun a => clMG (sOf n) g (by positivity) ((MG.enum (pA (sOf n))).symm a)
  let κ0 := mulCodingZ (sOf n) N k y
  intro c cluster gate
  let z := ((κ.coarsen χ).e γ).symm c
  let d := (κ.e cluster.val) (z cluster)
  have htop := selectedStateDigit_eq x n N k y hP pad hpad γ c cluster
  have hgate := clusterStateDigit_eq (sOf n) g (topEntryBits (sOf n) P) N k y (by positivity)
    (selected_Wc_le_bits (sOf n) P hP) cluster.val d gate
  change selectedStateDigit x n pad γ.val c.val cluster.val.val = d.val at htop
  change clusterStateDigit (sOf n) g (topEntryBits (sOf n) P) cluster.val.val d.val gate.val.val =
    ((κ0.e gate.val) (((κ0.coarsen χ0).e cluster.val).symm d gate)).val at hgate
  have hd : ((κ0.coarsen χ0).e cluster.val).symm d = z cluster := by
    change (κ.e cluster.val).symm ((κ.e cluster.val) (z cluster)) = z cluster
    exact Equiv.symm_apply_apply _ _
  rw [hd] at hgate
  unfold selectedGateCode
  dsimp only
  rw [htop]
  exact hgate

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedGateCode_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedGateCode_eq
