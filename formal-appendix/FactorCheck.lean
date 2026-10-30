/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time numeric checking of one coarsened factor, given polynomial-time
numeric widths, membership and side checks with a unary width bound.
-/
import PairOffset

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def factorCheck (g q : ℕ) (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool)
    (side : ℕ → ℕ → ℕ → Bool) : Bool :=
  (List.range (g * g)).all fun t =>
    let a := t / g % g
    let b := t % g
    !used a b || side a b (q / 2 ^ pairOffset g t width used)

theorem factorCheck_polyT {α : Type} [Enc α] {g W : α → ℕ} {q : α → ℕ}
    {width : α → ℕ → ℕ → ℕ} {used : α → ℕ → ℕ → Bool} {side : α → ℕ → ℕ → ℕ → Bool}
    (hg : Unary g) (hW : Unary W) (hq : PolyT q)
    (hw : PolyT fun p : (α × ℕ) × ℕ => width p.1.1 p.1.2 p.2)
    (hu : PolyT fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
    (hs : PolyT fun p : ((α × ℕ) × ℕ) × ℕ => side p.1.1.1 p.1.1.2 p.1.2 p.2)
    (hbound : ∀ a i j, width a i j ≤ W a) :
    PolyT fun a => factorCheck (g a) (q a) (width a) (used a) (side a) := by
  let β := α × ℕ
  have hg' : Unary fun p : β => g p.1 := Unary.comp hg PolyT.id'.fst
  have hW' : Unary fun p : β => W p.1 := Unary.comp hW PolyT.id'.fst
  have hw' : PolyT fun p : (β × ℕ) × ℕ => width p.1.1.1 p.1.2 p.2 :=
    PolyT.comp (f := fun p : (α × ℕ) × ℕ => width p.1.1 p.1.2 p.2)
      (g := fun p : (β × ℕ) × ℕ => ((p.1.1.1, p.1.2), p.2)) hw
      (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst.fst PolyT.id'.fst.snd) PolyT.id'.snd)
  have hu' : PolyT fun p : (β × ℕ) × ℕ => used p.1.1.1 p.1.2 p.2 :=
    PolyT.comp (f := fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
      (g := fun p : (β × ℕ) × ℕ => ((p.1.1.1, p.1.2), p.2)) hu
      (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst.fst PolyT.id'.fst.snd) PolyT.id'.snd)
  have hoff : Unary fun p : β => pairOffset (g p.1) p.2 (width p.1) (used p.1) :=
    pairOffset_unary (g := fun p : β => g p.1) (W := fun p : β => W p.1) (k := fun p : β => p.2)
      (width := fun p i j => width p.1 i j) (used := fun p i j => used p.1 i j)
      hg' hW' PolyT.id'.snd hw' hu' (fun p i j => hbound p.1 i j)
  have hgp : PolyT fun p : β => g p.1 := hg'.poly
  have ha : PolyT fun p : β => p.2 / g p.1 % g p.1 := PolyT.mod (PolyT.div PolyT.id'.snd hgp) hgp
  have hb : PolyT fun p : β => p.2 % g p.1 := PolyT.mod PolyT.id'.snd hgp
  have hq' : PolyT fun p : β => q p.1 := PolyT.comp hq PolyT.id'.fst
  have hshift : PolyT fun p : β => q p.1 / 2 ^ pairOffset (g p.1) p.2 (width p.1) (used p.1) :=
    PolyT.div hq' (PolyT.two_pow hoff)
  have hp : PolyT fun p : β => (((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1),
      q p.1 / 2 ^ pairOffset (g p.1) p.2 (width p.1) (used p.1)) :=
    PolyT.prodMk (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst ha) hb) hshift
  have hs' : PolyT fun p : β => side p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1)
      (q p.1 / 2 ^ pairOffset (g p.1) p.2 (width p.1) (used p.1)) :=
    PolyT.comp (f := fun p : ((α × ℕ) × ℕ) × ℕ => side p.1.1.1 p.1.1.2 p.1.2 p.2)
      (g := fun p : β => (((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1),
        q p.1 / 2 ^ pairOffset (g p.1) p.2 (width p.1) (used p.1))) hs hp
  have huargs : PolyT fun p : β => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1) :=
    PolyT.prodMk (PolyT.prodMk PolyT.id'.fst ha) hb
  have hu'' : PolyT fun p : β => used p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) :=
    PolyT.comp (f := fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
      (g := fun p : β => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1)) hu huargs
  have htest := PolyT.or' (PolyT.not' hu'') hs'
  have hr := PolyT.range (Unary.mul hg hg)
  exact PolyT.all (f := fun a t =>
    !used a (t / g a % g a) (t % g a) ||
      side a (t / g a % g a) (t % g a) (q a / 2 ^ pairOffset (g a) t (width a) (used a))) htest hr

theorem factorCheck_true_iff (g q : ℕ) (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool)
    (side : ℕ → ℕ → ℕ → Bool) : factorCheck g q width used side = true ↔
      ∀ t < g * g, used (t / g % g) (t % g) = true →
        side (t / g % g) (t % g) (q / 2 ^ pairOffset g t width used) = true := by
  simp only [factorCheck, List.all_eq_true, List.mem_range, Bool.or_eq_true, Bool.not_eq_true]
  constructor
  · intro h t ht hu
    rcases h t ht with hfalse | htrue
    · simp_all
    · exact htrue
  · intro h t ht
    cases hu : used (t / g % g) (t % g) with
    | false => exact Or.inl rfl
    | true => exact Or.inr (h t ht hu)

private theorem binary_prod_one {ι : Type} [DecidableEq ι] (s : Finset ι) (f : ι → ℤ)
    (hf : ∀ i ∈ s, f i = 0 ∨ f i = 1) :
    (∏ i ∈ s, f i) = 1 ↔ ∀ i ∈ s, f i = 1 := by
  constructor
  · intro h i hi
    rcases hf i hi with h0 | h1
    · have hz : (∏ i ∈ s, f i) = 0 := Finset.prod_eq_zero hi h0
      omega
    · exact h1
  · intro h
    exact Finset.prod_eq_one h

theorem factorCheck_eq_left (g K : ℕ) (hg : 0 < g) (χ : Fin g → Fin K)
    (u v : Fin K) (N : Net ℤ g) (z : (N.coarsen χ).S u) (q : ℕ)
    (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool) (side : ℕ → ℕ → ℕ → Bool)
    (hw : ∀ a b : Fin g, width a.val b.val = N.c a b)
    (hu : ∀ a b : Fin g, used a.val b.val = decide (Net.UsedP χ u v (a, b)))
    (hs : ∀ (a b : Fin g) (q' : ℕ),
      (if side a.val b.val q' then (1 : ℤ) else 0) = (N.orient χ u v (a, b)).l z (bdec _ q')) :
    factorCheck g q width used side = true ↔ ((N.coarsen χ).F u v).l z (bdec _ q) = 1 := by
  rw [factorCheck_true_iff, Net.coarsen_l_code N hg χ u v z q]
  have hbinary : ∀ t ∈ Finset.range (g * g),
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).l z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 0 ∨
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).l z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 1 := by
    intro t ht
    by_cases h : Net.UsedP χ u v (Net.pairOf g hg t)
    · rw [if_pos h, ← hs (Net.pairOf g hg t).1 (Net.pairOf g hg t).2]
      split <;> simp
    · simp [h]
  rw [binary_prod_one _ _ hbinary]
  constructor
  · intro h t ht
    have hlt := Finset.mem_range.mp ht
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt hlt) χ u v N width used hw hu
    have hu' := hu a b
    have hs' := hs a b (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)
    by_cases hU : Net.UsedP χ u v (a, b)
    · have hcheck := h t hlt
      change used a.val b.val = true → side a.val b.val (q / 2 ^ pairOffset g t width used) = true at hcheck
      rw [ho] at hcheck
      have hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) = true :=
        hcheck (by simpa [hU] using hu')
      change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      rw [if_pos hU, ← hs', hside]
      rfl
    · change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      simp [hU]
  · intro h t ht hused
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt ht) χ u v N width used hw hu
    have hu' := hu a b
    have hU : Net.UsedP χ u v (a, b) := by
      have : decide (Net.UsedP χ u v (a, b)) = true := by simpa [a, b, Net.pairOf] using hu'.symm.trans hused
      exact of_decide_eq_true this
    have hprod := h t (Finset.mem_range.mpr ht)
    change (if Net.UsedP χ u v (a, b) then _ else 1) = 1 at hprod
    rw [if_pos hU, ← hs a b] at hprod
    change side a.val b.val (q / 2 ^ pairOffset g t width used) = true
    rw [ho]
    cases hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) <;> simp_all

theorem factorCheck_eq_left_used (g K : ℕ) (hg : 0 < g) (χ : Fin g → Fin K)
    (u v : Fin K) (N : Net ℤ g) (z : (N.coarsen χ).S u) (q : ℕ)
    (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool) (side : ℕ → ℕ → ℕ → Bool)
    (hw : ∀ a b : Fin g, width a.val b.val = N.c a b)
    (hu : ∀ a b : Fin g, used a.val b.val = decide (Net.UsedP χ u v (a, b)))
    (hs : ∀ (a b : Fin g), Net.UsedP χ u v (a, b) → ∀ (q' : ℕ),
      (if side a.val b.val q' then (1 : ℤ) else 0) = (N.orient χ u v (a, b)).l z (bdec _ q')) :
    factorCheck g q width used side = true ↔ ((N.coarsen χ).F u v).l z (bdec _ q) = 1 := by
  rw [factorCheck_true_iff, Net.coarsen_l_code N hg χ u v z q]
  have hbinary : ∀ t ∈ Finset.range (g * g),
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).l z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 0 ∨
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).l z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 1 := by
    intro t ht
    by_cases h : Net.UsedP χ u v (Net.pairOf g hg t)
    · rw [if_pos h, ← hs (Net.pairOf g hg t).1 (Net.pairOf g hg t).2 h]
      split <;> simp
    · simp [h]
  rw [binary_prod_one _ _ hbinary]
  constructor
  · intro h t ht
    have hlt := Finset.mem_range.mp ht
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt hlt) χ u v N width used hw hu
    have hu' := hu a b
    by_cases hU : Net.UsedP χ u v (a, b)
    · have hs' := hs a b hU (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)
      have hcheck := h t hlt
      change used a.val b.val = true → side a.val b.val (q / 2 ^ pairOffset g t width used) = true at hcheck
      rw [ho] at hcheck
      have hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) = true :=
        hcheck (by simpa [hU] using hu')
      change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      rw [if_pos hU, ← hs', hside]
      rfl
    · change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      simp [hU]
  · intro h t ht hused
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt ht) χ u v N width used hw hu
    have hu' := hu a b
    have hU : Net.UsedP χ u v (a, b) := by
      have : decide (Net.UsedP χ u v (a, b)) = true := by simpa [a, b, Net.pairOf] using hu'.symm.trans hused
      exact of_decide_eq_true this
    have hprod := h t (Finset.mem_range.mpr ht)
    change (if Net.UsedP χ u v (a, b) then _ else 1) = 1 at hprod
    rw [if_pos hU, ← hs a b hU] at hprod
    change side a.val b.val (q / 2 ^ pairOffset g t width used) = true
    rw [ho]
    cases hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) <;> simp_all


theorem factorCheck_eq_right (g K : ℕ) (hg : 0 < g) (χ : Fin g → Fin K)
    (u v : Fin K) (N : Net ℤ g) (z : (N.coarsen χ).S v) (q : ℕ)
    (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool) (side : ℕ → ℕ → ℕ → Bool)
    (hw : ∀ a b : Fin g, width a.val b.val = N.c a b)
    (hu : ∀ a b : Fin g, used a.val b.val = decide (Net.UsedP χ u v (a, b)))
    (hs : ∀ (a b : Fin g) (q' : ℕ),
      (if side a.val b.val q' then (1 : ℤ) else 0) = (N.orient χ u v (a, b)).r z (bdec _ q')) :
    factorCheck g q width used side = true ↔ ((N.coarsen χ).F u v).r z (bdec _ q) = 1 := by
  rw [factorCheck_true_iff, Net.coarsen_r_code N hg χ u v z q]
  have hbinary : ∀ t ∈ Finset.range (g * g),
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).r z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 0 ∨
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).r z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 1 := by
    intro t ht
    by_cases h : Net.UsedP χ u v (Net.pairOf g hg t)
    · rw [if_pos h, ← hs (Net.pairOf g hg t).1 (Net.pairOf g hg t).2]
      split <;> simp
    · simp [h]
  rw [binary_prod_one _ _ hbinary]
  constructor
  · intro h t ht
    have hlt := Finset.mem_range.mp ht
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt hlt) χ u v N width used hw hu
    have hu' := hu a b
    have hs' := hs a b (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)
    by_cases hU : Net.UsedP χ u v (a, b)
    · have hcheck := h t hlt
      change used a.val b.val = true → side a.val b.val (q / 2 ^ pairOffset g t width used) = true at hcheck
      rw [ho] at hcheck
      have hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) = true :=
        hcheck (by simpa [hU] using hu')
      change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      rw [if_pos hU, ← hs', hside]
      rfl
    · change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      simp [hU]
  · intro h t ht hused
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt ht) χ u v N width used hw hu
    have hu' := hu a b
    have hU : Net.UsedP χ u v (a, b) := by
      have : decide (Net.UsedP χ u v (a, b)) = true := by simpa [a, b, Net.pairOf] using hu'.symm.trans hused
      exact of_decide_eq_true this
    have hprod := h t (Finset.mem_range.mpr ht)
    change (if Net.UsedP χ u v (a, b) then _ else 1) = 1 at hprod
    rw [if_pos hU, ← hs a b] at hprod
    change side a.val b.val (q / 2 ^ pairOffset g t width used) = true
    rw [ho]
    cases hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) <;> simp_all


theorem factorCheck_eq_right_used (g K : ℕ) (hg : 0 < g) (χ : Fin g → Fin K)
    (u v : Fin K) (N : Net ℤ g) (z : (N.coarsen χ).S v) (q : ℕ)
    (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool) (side : ℕ → ℕ → ℕ → Bool)
    (hw : ∀ a b : Fin g, width a.val b.val = N.c a b)
    (hu : ∀ a b : Fin g, used a.val b.val = decide (Net.UsedP χ u v (a, b)))
    (hs : ∀ (a b : Fin g), Net.UsedP χ u v (a, b) → ∀ (q' : ℕ),
      (if side a.val b.val q' then (1 : ℤ) else 0) = (N.orient χ u v (a, b)).r z (bdec _ q')) :
    factorCheck g q width used side = true ↔ ((N.coarsen χ).F u v).r z (bdec _ q) = 1 := by
  rw [factorCheck_true_iff, Net.coarsen_r_code N hg χ u v z q]
  have hbinary : ∀ t ∈ Finset.range (g * g),
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).r z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 0 ∨
      (if Net.UsedP χ u v (Net.pairOf g hg t) then
        (N.orient χ u v (Net.pairOf g hg t)).r z
          (bdec _ (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)) else 1) = 1 := by
    intro t ht
    by_cases h : Net.UsedP χ u v (Net.pairOf g hg t)
    · rw [if_pos h, ← hs (Net.pairOf g hg t).1 (Net.pairOf g hg t).2 h]
      split <;> simp
    · simp [h]
  rw [binary_prod_one _ _ hbinary]
  constructor
  · intro h t ht
    have hlt := Finset.mem_range.mp ht
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt hlt) χ u v N width used hw hu
    have hu' := hu a b
    by_cases hU : Net.UsedP χ u v (a, b)
    · have hs' := hs a b hU (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t)
      have hcheck := h t hlt
      change used a.val b.val = true → side a.val b.val (q / 2 ^ pairOffset g t width used) = true at hcheck
      rw [ho] at hcheck
      have hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) = true :=
        hcheck (by simpa [hU] using hu')
      change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      rw [if_pos hU, ← hs', hside]
      rfl
    · change (if Net.UsedP χ u v (a, b) then _ else 1) = 1
      simp [hU]
  · intro h t ht hused
    let a := (Net.pairOf g hg t).1
    let b := (Net.pairOf g hg t).2
    have ho := pairOffset_eq g K t hg (Nat.le_of_lt ht) χ u v N width used hw hu
    have hu' := hu a b
    have hU : Net.UsedP χ u v (a, b) := by
      have : decide (Net.UsedP χ u v (a, b)) = true := by simpa [a, b, Net.pairOf] using hu'.symm.trans hused
      exact of_decide_eq_true this
    have hprod := h t (Finset.mem_range.mpr ht)
    change (if Net.UsedP χ u v (a, b) then _ else 1) = 1 at hprod
    rw [if_pos hU, ← hs a b hU] at hprod
    change side a.val b.val (q / 2 ^ pairOffset g t width used) = true
    rw [ho]
    cases hside : side a.val b.val (q / 2 ^ offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) t) <;> simp_all

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.factorCheck_polyT
#print axioms IntegerFactoring.HalfAppendix.factorCheck_eq_left
#print axioms IntegerFactoring.HalfAppendix.factorCheck_eq_right

#print axioms IntegerFactoring.HalfAppendix.factorCheck_eq_right_used

#print axioms IntegerFactoring.HalfAppendix.factorCheck_eq_left_used
