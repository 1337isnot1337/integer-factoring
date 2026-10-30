/-
SPDX-License-Identifier: Apache-2.0

Executable extraction of a cluster-state code from a selected top-group code.
Individual multiplication-gate settings and row summands remain to be decoded.
-/
import SelectedRadix

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def selectedStateDigit (x n : ℕ) (pad : List Unit) (γ code i : ℕ) : ℕ :=
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  let B := topEntryBits (sOf n) P
  cappedDigit (2 ^ B) g code i (selectedRadix x n pad γ)

theorem selectedStateDigit_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n γ code : α → ℕ}
    (hn : Unary n) (hγ : PolyT γ) (hc : PolyT code) :
    PolyT fun p : (α × List Unit) × ℕ => selectedStateDigit x (n p.1.1) p.1.2 (γ p.1.1) (code p.1.1) p.2 := by
  let β := α × List Unit
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let B : β → ℕ := fun p => topEntryBits (sOf (n p.1)) (P p)
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hB : Unary B := Unary.comp (selected_topEntryBits_unary x hx hn) PolyT.id'.fst
  have hcode : PolyT fun p : β => code p.1 := PolyT.comp hc PolyT.id'.fst
  have hr : PolyT fun p : β × ℕ => selectedRadix x (n p.1.1) p.1.2 (γ p.1.1) p.2 :=
    selectedRadix_polyT x hx hn hγ
  exact cappedDigit_polyT (B := fun p : β => 2 ^ B p) (g := g) (code := fun p : β => code p.1)
    (radix := fun p i => selectedRadix x (n p.1) p.2 (γ p.1) i) (PolyT.two_pow hB) hg hcode hr

theorem selectedStateDigit_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (a : {a : Fin g // χ a = γ}),
      selectedStateDigit x n pad γ.val c.val a.val.val =
        ((κ.e a.val) (((κ.coarsen χ).e γ).symm c a)).val := by
  dsimp only
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  let χ := selectedTopGrouping (sOf n) N k y P hP
  let κ := clusterCodingZ (sOf n) g N k y (by positivity)
  intro c a
  have hrad : selectedRadix x n pad γ.val = Coding.sg κ χ γ := by
    funext i
    exact selectedRadix_eq x n N k y hP pad hpad γ i
  have hpos : ∀ i, 0 < Coding.sg κ χ γ i := Coding.sg_pos κ χ γ
  have hmono : Monotone (pp (Coding.sg κ χ γ)) := by
    apply monotone_nat_of_le_succ
    intro i
    exact Nat.le_mul_of_pos_right _ (hpos i)
  have hfull : pp (Coding.sg κ χ γ) g = (κ.coarsen χ).s γ := by
    change pp (Coding.sg κ χ γ) g = if h : γ.val < P then pp (Coding.sg κ χ ⟨γ.val, h⟩) g else 1
    rw [dif_pos γ.isLt]
  have hfullBound : (κ.coarsen χ).s γ ≤ 2 ^ topEntryBits (sOf n) P :=
    selected_explicit_code_bound (sOf n) N k y P hP γ
  have hprefix : pp (Coding.sg κ χ γ) a.val.val ≤ 2 ^ topEntryBits (sOf n) P :=
    (hmono (Nat.le_of_lt a.val.isLt)).trans (by rw [hfull]; exact hfullBound)
  unfold selectedStateDigit
  dsimp only
  rw [cappedDigit_eq _ _ _ _ _ (Nat.le_of_lt a.val.isLt) (by
    intro i _
    rw [hrad]
    exact hpos i) (by rw [hrad]; exact hprefix), hrad]
  exact (Coding.coarsen_digit κ χ γ c a).symm

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedStateDigit_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedStateDigit_eq
