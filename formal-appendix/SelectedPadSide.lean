/-
SPDX-License-Identifier: Apache-2.0

The padded factor-side check used by a selected clique row.
-/
import SelectedTopSide
import SelectedTopEdge

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedPadSide (x n : ℕ) (pad : List Unit) (γ code δ q : ℕ) : Bool :=
  let P := K x (selectedLevel x n)
  let L := topWidth (sOf n) P + 1
  let w := selectedTopEdge x n pad γ δ
  decide (w ≤ L) && selectedTopLeft x n pad γ code δ q && decide (q / 2 ^ w % 2 ^ (L - w) = 0)

set_option maxHeartbeats 2000000 in
theorem selectedPadSide_polyT (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α]
    {n γ code δ q : α → ℕ} (hn : Unary n) (hγ : PolyT γ) (hc : PolyT code)
    (hδ : PolyT δ) (hq : PolyT q) :
    PolyT fun p : α × List Unit => selectedPadSide x (n p.1) p.2 (γ p.1) (code p.1) (δ p.1) (q p.1) := by
  let β := α × List Unit
  let L : β → ℕ := fun p => topWidth (sOf (n p.1)) (K x (selectedLevel x (n p.1))) + 1
  let w : β → ℕ := fun p => selectedTopEdge x (n p.1) p.2 (γ p.1) (δ p.1)
  have hL : Unary L := Unary.comp (selected_topWidth_unary x hx hn) PolyT.id'.fst
  have hw : Unary w := selectedTopEdge_unary x hx hn hγ hδ
  have hsubp : PolyT fun p : β => L p - w p := PolyT.sub hL.poly hw.poly
  have hsub : Unary fun p : β => L p - w p := Unary.of_le hL hsubp (fun p => Nat.sub_le _ _)
  have hq' : PolyT fun p : β => q p.1 := PolyT.comp hq PolyT.id'.fst
  have hrem : PolyT fun p : β => q p.1 / 2 ^ w p % 2 ^ (L p - w p) :=
    PolyT.mod (PolyT.div hq' (PolyT.two_pow hw)) (PolyT.two_pow hsub)
  exact PolyT.and' (PolyT.and' (PolyT.decide_le hw.poly hL.poly) (selectedTopLeft_polyT x hx hn hγ hc hδ hq))
    (PolyT.decide_eq hrem (PolyT.const 0))


set_option maxHeartbeats 2000000 in
theorem selectedPadSide_eq_left (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ δ : Fin (K x (selectedLevel x n))) (hne : γ ≠ δ) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let L := topWidth (sOf n) P + 1
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    let M := (clusterNet (sOf n) ℤ N k y g (by positivity)).coarsen χ
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (q : ℕ),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedPadSide x n pad γ.val c.val δ.val q = true ↔ (M.padF L γ δ).l z (bdec L q) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let L := topWidth s P + 1
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let M := (clusterNet s ℤ N k y g (by positivity)).coarsen χ
  intro c q
  let z := ((κ.coarsen χ).e γ).symm c
  have he := selectedTopLeft_eq x n N k y hP pad hpad γ δ hne c q
  change selectedTopLeft x n pad γ.val c.val δ.val q = true ↔ (M.F γ δ).l z (bdec _ q) = 1 at he
  have hw := selectedTopEdge_eq x n N k y hP pad hpad γ δ
  change selectedTopEdge x n pad γ.val δ.val = M.c γ δ at hw
  unfold selectedPadSide
  dsimp only
  rw [hw]
  simp only [Bool.and_eq_true, decide_eq_true_eq]
  change ((M.c γ δ ≤ L ∧ selectedTopLeft x n pad γ.val c.val δ.val q = true) ∧
    q / 2 ^ M.c γ δ % 2 ^ (L - M.c γ δ) = 0) ↔ (M.padF L γ δ).l z (bdec L q) = 1
  by_cases hf : M.c γ δ ≤ L
  · simp only [hf, true_and, Net.padF, dif_pos hf]
    simp only [dif_pos True.intro]
    rw [Fac.pad_l_bdec]
    by_cases hz : q / 2 ^ M.c γ δ % 2 ^ (L - M.c γ δ) = 0
    · simp [hz, he]
    · simp [hz]
  · simp [hf, Net.padF, Fac.zero]

set_option maxHeartbeats 2000000 in
theorem selectedPadSide_eq_right (x n N k y : ℕ) (hP : 2 ≤ K x (selectedLevel x n))
    (pad : List Unit) (hpad : (K x (selectedLevel x n)) ^ (3 * (K x (selectedLevel x n)) ^ 8) ≤ pad.length)
    (γ δ : Fin (K x (selectedLevel x n))) (hne : γ ≠ δ) :
    let P := K x (selectedLevel x n)
    let g := 3 * P ^ 8
    let L := topWidth (sOf n) P + 1
    let χ := selectedTopGrouping (sOf n) N k y P hP
    let κ := clusterCodingZ (sOf n) g N k y (by positivity)
    let M := (clusterNet (sOf n) ℤ N k y g (by positivity)).coarsen χ
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (q : ℕ),
      let z := ((κ.coarsen χ).e γ).symm c
      selectedPadSide x n pad γ.val c.val δ.val q = true ↔ (M.padF L δ γ).r z (bdec L q) = 1 := by
  dsimp only
  let P := K x (selectedLevel x n)
  let s := sOf n
  let g := 3 * P ^ 8
  let L := topWidth s P + 1
  let χ := selectedTopGrouping s N k y P hP
  let κ := clusterCodingZ s g N k y (by positivity)
  let M := (clusterNet s ℤ N k y g (by positivity)).coarsen χ
  intro c q
  let z := ((κ.coarsen χ).e γ).symm c
  have he := selectedTopRight_eq x n N k y hP pad hpad γ δ hne c q
  change selectedTopLeft x n pad γ.val c.val δ.val q = true ↔ (M.F δ γ).r z (bdec _ q) = 1 at he
  have hw := selectedTopEdge_eq x n N k y hP pad hpad δ γ
  rw [selectedTopEdge_symm] at hw
  change selectedTopEdge x n pad γ.val δ.val = M.c δ γ at hw
  unfold selectedPadSide
  dsimp only
  rw [hw]
  simp only [Bool.and_eq_true, decide_eq_true_eq]
  change ((M.c δ γ ≤ L ∧ selectedTopLeft x n pad γ.val c.val δ.val q = true) ∧
    q / 2 ^ M.c δ γ % 2 ^ (L - M.c δ γ) = 0) ↔ (M.padF L δ γ).r z (bdec L q) = 1
  by_cases hf : M.c δ γ ≤ L
  · simp only [hf, true_and, Net.padF, dif_pos hf]
    simp only [dif_pos True.intro]
    rw [Fac.pad_r_bdec]
    by_cases hz : q / 2 ^ M.c δ γ % 2 ^ (L - M.c δ γ) = 0
    · simp [hz, he]
    · simp [hz]
  · simp [hf, Net.padF, Fac.zero]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedPadSide_polyT

#print axioms IntegerFactoring.HalfAppendix.selectedPadSide_eq_left
#print axioms IntegerFactoring.HalfAppendix.selectedPadSide_eq_right
