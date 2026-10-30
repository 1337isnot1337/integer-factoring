/-
SPDX-License-Identifier: Apache-2.0

Executable selected-group mixed-radix sizes, conditional on sufficient pad.
-/
import ClusterSize
import SelectedGroupingPoly

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def selectedRadix (x n : ℕ) (pad : List Unit) (γ i : ℕ) : ℕ :=
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  if i < g ∧ selectedGroupDigit x n pad i = γ then
    cappedClusterSize (sOf n) g (topEntryBits (sOf n) P) i
  else 1

theorem selected_Wc_le_bits (s P : ℕ) (hP : 2 ≤ P) :
    Wc s (3 * P ^ 8) ≤ topEntryBits s P := by
  have hp : P ≤ P ^ 8 := Nat.le_self_pow (by norm_num) _
  have hd : 1 ≤ 2 * (3 * P ^ 8) / P := by
    apply (Nat.le_div_iff_mul_le (by omega : 0 < P)).2
    nlinarith
  unfold topEntryBits
  nlinarith

theorem selectedRadix_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n γ : α → ℕ}
    (hn : Unary n) (hγ : PolyT γ) :
    PolyT fun p : (α × List Unit) × ℕ => selectedRadix x (n p.1.1) p.1.2 (γ p.1.1) p.2 := by
  let β := α × List Unit
  let P : β → ℕ := fun p => K x (selectedLevel x (n p.1))
  let g : β → ℕ := fun p => 3 * P p ^ 8
  let B : β → ℕ := fun p => topEntryBits (sOf (n p.1)) (P p)
  have hP : Unary P := Unary.comp (selected_K_unary x hx hn) PolyT.id'.fst
  have hg : Unary g := Unary.mul (Unary.const 3) (unary_pow_const hP 8)
  have hB : Unary B := Unary.comp (selected_topEntryBits_unary x hx hn) PolyT.id'.fst
  have hn' : Unary fun p : β × ℕ => n p.1.1 := Unary.comp hn PolyT.id'.fst.fst
  have hg' : PolyT fun p : β × ℕ => g p.1 := PolyT.comp hg.poly PolyT.id'.fst
  have hB' : Unary fun p : β × ℕ => B p.1 := Unary.comp hB PolyT.id'.fst
  have hc : PolyT fun p : β × ℕ => selectedGroupDigit x (n p.1.1) p.1.2 p.2 := selectedGroupDigit_polyT x hx hn
  have hγ' : PolyT fun p : β × ℕ => γ p.1.1 := PolyT.comp hγ PolyT.id'.fst.fst
  have hs : PolyT fun p : β × ℕ => cappedClusterSize (sOf (n p.1.1)) (g p.1) (B p.1) p.2 :=
    cappedClusterSize_polyT hn' hg' hB' PolyT.id'.snd
  have hlt := PolyT.decide_lt (PolyT.id'.snd : PolyT fun p : β × ℕ => p.2) hg'
  have heq := PolyT.decide_eq hc hγ'
  exact PolyT.ite' (PolyT.decide_and hlt heq) hs (PolyT.const 1)

theorem selectedRadix_eq (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ : Fin (K x (selectedLevel x n))) (i : ℕ) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    selectedRadix x n pad γ.val i = Coding.sg (clusterCodingZ (sOf n) g N k y (by positivity))
      (selectedTopGrouping (sOf n) N k y P hP) γ i := by
  dsimp only
  let P := K x (selectedLevel x n)
  let g := 3 * P ^ 8
  by_cases hi : i < g
  · let j : Fin g := ⟨i, hi⟩
    have hd : selectedGroupDigit x n pad i = (selectedTopGrouping (sOf n) N k y P hP j).val :=
      selectedGroupDigit_eq x n N k y hP pad hpad j
    have hs : cappedClusterSize (sOf n) g (topEntryBits (sOf n) P) i =
        (clusterCodingZ (sOf n) g N k y (by positivity)).s j :=
      cappedClusterSize_eq (sOf n) g (topEntryBits (sOf n) P) N k y (by positivity) j (selected_Wc_le_bits (sOf n) P hP)
    simp [selectedRadix, Coding.sg, P, g, hi, hd, hs, Fin.ext_iff, j]
  · simp [selectedRadix, Coding.sg, P, g, hi]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedRadix_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedRadix_eq
