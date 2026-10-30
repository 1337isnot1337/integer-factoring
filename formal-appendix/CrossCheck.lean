/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time checks of the cross-pair factors in a coarsened network.
-/
import WeightCheck

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def crossCheck (g : ℕ) (used pair : ℕ → ℕ → Bool) : Bool :=
  weightCheck g (fun _ => false) (fun _ => true) used pair

theorem crossCheck_polyT {α : Type} [Enc α] {g : α → ℕ}
    {used pair : α → ℕ → ℕ → Bool} (hg : Unary g)
    (hu : PolyT fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
    (hp : PolyT fun p : (α × ℕ) × ℕ => pair p.1.1 p.1.2 p.2) :
    PolyT fun a => crossCheck (g a) (used a) (pair a) := by
  exact weightCheck_polyT hg (PolyT.const false) (PolyT.const true) hu hp

private theorem list_prod_binary_iff {α : Type} (l : List α) (f : α → ℤ)
    (hf : ∀ a ∈ l, f a = 0 ∨ f a = 1) :
    (l.map f).prod = 1 ↔ ∀ a ∈ l, f a = 1 := by
  induction l with
  | nil => simp
  | cons a l ih =>
    have ha := hf a (by simp)
    have hl : ∀ b ∈ l, f b = 0 ∨ f b = 1 := by
      intro b hb
      exact hf b (by simp [hb])
    have htail : (l.map f).prod = 0 ∨ (l.map f).prod = 1 := by
      by_cases h : (l.map f).prod = 1
      · exact Or.inr h
      · have hnot : ¬ ∀ b ∈ l, f b = 1 := by simpa [ih hl] using h
        push Not at hnot
        obtain ⟨b, hb, hne⟩ := hnot
        have hz : f b = 0 := (hl b hb).resolve_right hne
        have : 0 ∈ l.map f := by exact List.mem_map.mpr ⟨b, hb, hz⟩
        exact Or.inl (List.prod_eq_zero this)
    simp only [List.map_cons, List.prod_cons, List.mem_cons, forall_eq_or_imp]
    rcases ha with ha | ha <;> rcases htail with ht | ht <;> simp_all [ih hl]

theorem crossCheck_eq (g K : ℕ) (hg : 0 < g) (χ : Fin g → Fin K)
    (u v : Fin K) (N : Net ℤ g) (z : (N.coarsen χ).S u) (w : (N.coarsen χ).S v)
    (used pair : ℕ → ℕ → Bool)
    (hu : ∀ a b : Fin g, used a.val b.val = decide (Net.UsedP χ u v (a,b)))
    (hp : ∀ a b : Fin g, Net.UsedP χ u v (a,b) →
      (if pair a.val b.val then (1 : ℤ) else 0) = (N.orient χ u v (a,b)).eval z w) :
    crossCheck g used pair = true ↔ (N.coarsen χ).Ψ u v z w = 1 := by
  rw [crossCheck, weightCheck_true_iff g hg]
  simp only [Bool.false_eq_true, IsEmpty.forall_iff]
  have htr : (∀ i < g, True) := by simp
  rw [show (∀ i < g, True) ↔ True from iff_true_intro htr]
  simp only [true_and]
  change (∀ a < g, ∀ b < g, used a b = true → pair a b = true) ↔ _
  rw [Net.Ψ, Net.coarsen_F_eval]
  let l := Net.crossList χ u v
  let f := fun p : Fin g × Fin g => (N.orient χ u v p).eval z w
  have hmem : ∀ p : Fin g × Fin g, p ∈ l ↔ Net.UsedP χ u v p := by
    intro p
    simp [l, Net.crossList, Net.UsedP, mem_pairsList]
  have hbinary : ∀ p ∈ l, f p = 0 ∨ f p = 1 := by
    intro p hpl
    have h := hp p.1 p.2 ((hmem p).mp hpl)
    cases hpair : pair p.1.val p.2.val <;> simp_all [f]
  rw [list_prod_binary_iff l f hbinary]
  constructor
  · intro h p hpl
    have huse := (hmem p).mp hpl
    have hv : used p.1.val p.2.val = true := by rw [hu]; simp [huse]
    have hv' := h p.1.val p.1.isLt p.2.val p.2.isLt hv
    have he := hp p.1 p.2 huse
    simpa [f, hv'] using he.symm
  · intro h a ha b hb huse
    let p : Fin g × Fin g := (⟨a, ha⟩, ⟨b, hb⟩)
    have hpuse : Net.UsedP χ u v p := by
      have he := hu p.1 p.2
      simp only [p] at he
      rw [huse] at he
      exact of_decide_eq_true he.symm
    have hv := h p ((hmem p).mpr hpuse)
    have he := hp p.1 p.2 hpuse
    change (if pair a b then (1 : ℤ) else 0) = f p at he
    rw [hv] at he
    cases hpair : pair a b <;> simp_all

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.crossCheck_polyT
#print axioms IntegerFactoring.HalfAppendix.crossCheck_eq
