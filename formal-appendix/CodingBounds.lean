/-
SPDX-License-Identifier: Apache-2.0

A member radix is bounded by its enclosing mixed-radix state space.
-/
import ExplicitCoding

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem coding_member_size_le {R : Type} [CommRing R] {g K : ℕ} {N : Net R g}
    (κ : Coding N) (χ : Fin g → Fin K) (γ : Fin K) (a : Fin g) (ha : χ a = γ) :
    κ.s a ≤ (κ.coarsen χ).s γ := by
  let rad := Coding.sg κ χ γ
  have hpos : ∀ i, 0 < rad i := Coding.sg_pos κ χ γ
  have hp : 0 < pp rad a.val := pp_pos (fun i _ => hpos i)
  have ht : 0 < pp rad g := pp_pos (fun i _ => hpos i)
  have hd := pp_dvd rad (Nat.le_of_lt a.isLt)
  have heq : rad a.val = κ.s a := Coding.sg_mem κ χ γ ha
  have hdiv : pp rad a.val * κ.s a ∣ pp rad g := by
    rcases hd with hd | hd
    · simpa [heq] using hd
    · omega
  have hfactor : κ.s a ≤ pp rad a.val * κ.s a := Nat.le_mul_of_pos_left _ hp
  have htotal := Nat.le_of_dvd ht hdiv
  change κ.s a ≤ if h : γ.val < K then pp (Coding.sg κ χ ⟨γ.val, h⟩) g else 1
  rw [dif_pos γ.isLt]
  exact hfactor.trans htotal

theorem selected_gate_size_bound (s N k y P : ℕ) (hP : 2 ≤ P)
    (γ : Fin P) :
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping s N k y P hP
    let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
    ∀ (cluster : Fin g) (gate : Fin (mulNodes s)), χ cluster = γ → χ0 gate = cluster →
      mulSz s gate.val ≤ 2 ^ topEntryBits s P := by
  dsimp only
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping s N k y P hP
  let χ0 := fun a => clMG s g (by positivity) ((MG.enum (pA s)).symm a)
  let κ0 := mulCodingZ s N k y
  let κ := clusterCodingZ s g N k y (by positivity)
  intro cluster gate hcluster hgate
  have h1 : κ0.s gate ≤ κ.s cluster := coding_member_size_le κ0 χ0 cluster gate hgate
  have h2 : κ.s cluster ≤ (κ.coarsen χ).s γ := coding_member_size_le κ χ γ cluster hcluster
  have h3 : (κ.coarsen χ).s γ ≤ 2 ^ topEntryBits s P := selected_explicit_code_bound s N k y P hP γ
  exact h1.trans (h2.trans h3)

theorem selected_gate_kind_size_bound (s N k y P : ℕ) (hP : 2 ≤ P)
    (γ : Fin P) (cluster : Fin (3 * P ^ 8)) (gate : MG (pA s))
    (hc : selectedTopGrouping s N k y P hP cluster = γ)
    (hg : clMG s (3 * P ^ 8) (by positivity) gate = cluster) :
    gsz s gate ≤ 2 ^ topEntryBits s P := by
  have h := selected_gate_size_bound s N k y P hP γ cluster (MG.enum (pA s) gate) hc (by simpa using hg)
  have hs := mulSz_fin s (show Fin (mulNodes s) from MG.enum (pA s) gate)
  rw [Equiv.symm_apply_apply] at hs
  rw [hs] at h
  exact h

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.coding_member_size_le
#print axioms IntegerFactoring.HalfAppendix.selected_gate_size_bound
#print axioms IntegerFactoring.HalfAppendix.selected_gate_kind_size_bound
