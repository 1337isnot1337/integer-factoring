/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time numeric checks for the node and internal-pair factors of a
coarsened network. The supplied callbacks carry the underlying semantics.
-/
import FactorCheck

open Groundwork Factoring.Aux Finset

namespace IntegerFactoring.HalfAppendix

def weightCheck (g : ℕ) (member node : ℕ → Bool) (used pair : ℕ → ℕ → Bool) : Bool :=
  (List.range g).all (fun i => !member i || node i) &&
  (List.range (g * g)).all (fun t => !used (t / g % g) (t % g) || pair (t / g % g) (t % g))

theorem weightCheck_polyT {α : Type} [Enc α] {g : α → ℕ}
    {member node : α → ℕ → Bool} {used pair : α → ℕ → ℕ → Bool}
    (hg : Unary g) (hm : PolyT fun p : α × ℕ => member p.1 p.2)
    (hn : PolyT fun p : α × ℕ => node p.1 p.2)
    (hu : PolyT fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
    (hp : PolyT fun p : (α × ℕ) × ℕ => pair p.1.1 p.1.2 p.2) :
    PolyT fun a => weightCheck (g a) (member a) (node a) (used a) (pair a) := by
  have hnode := PolyT.or' (PolyT.not' hm) hn
  have hnodes := PolyT.all (f := fun a i => !member a i || node a i) hnode (PolyT.range hg)
  have hg' : PolyT fun p : α × ℕ => g p.1 := PolyT.comp hg.poly PolyT.id'.fst
  have ha : PolyT fun p : α × ℕ => p.2 / g p.1 % g p.1 := PolyT.mod (PolyT.div PolyT.id'.snd hg') hg'
  have hb : PolyT fun p : α × ℕ => p.2 % g p.1 := PolyT.mod PolyT.id'.snd hg'
  have hargs : PolyT fun p : α × ℕ => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1) :=
    PolyT.prodMk (PolyT.prodMk PolyT.id'.fst ha) hb
  have hu' : PolyT fun p : α × ℕ => used p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) :=
    PolyT.comp (f := fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
      (g := fun p : α × ℕ => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1)) hu hargs
  have hp' : PolyT fun p : α × ℕ => pair p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) :=
    PolyT.comp (f := fun p : (α × ℕ) × ℕ => pair p.1.1 p.1.2 p.2)
      (g := fun p : α × ℕ => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1)) hp hargs
  have hpair := PolyT.or' (PolyT.not' hu') hp'
  have hpairs := PolyT.all (f := fun a t => !used a (t / g a % g a) (t % g a) || pair a (t / g a % g a) (t % g a))
    hpair (PolyT.range (Unary.mul hg hg))
  exact PolyT.and' hnodes hpairs

theorem weightCheck_true_iff (g : ℕ) (hg : 0 < g) (member node : ℕ → Bool) (used pair : ℕ → ℕ → Bool) :
    weightCheck g member node used pair = true ↔
      (∀ i < g, member i = true → node i = true) ∧
      (∀ a < g, ∀ b < g, used a b = true → pair a b = true) := by
  simp only [weightCheck, Bool.and_eq_true, List.all_eq_true, List.mem_range, Bool.or_eq_true, Bool.not_eq_true]
  constructor
  · rintro ⟨hnode, hpair⟩
    constructor
    · intro i hi hm
      rcases hnode i hi with hfalse | htrue
      · simp_all
      · exact htrue
    · intro a ha b hb hu
      let t := a * g + b
      have ht : t < g * g := by dsimp [t]; nlinarith
      have ht1 : t / g % g = a := by
        dsimp [t]
        rw [show a * g + b = b + g * a by ring]
        rw [Nat.add_mul_div_left _ _ hg, Nat.div_eq_of_lt hb, zero_add, Nat.mod_eq_of_lt ha]
      have ht2 : t % g = b := by dsimp [t]; simp [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hb]
      have h := hpair t ht
      rw [ht1, ht2] at h
      rcases h with hfalse | htrue
      · simp_all
      · exact htrue
  · rintro ⟨hnode, hpair⟩
    constructor
    · intro i hi
      cases hm : member i with
      | false => exact Or.inl rfl
      | true => exact Or.inr (hnode i hi hm)
    · intro t ht
      have ha : t / g % g < g := Nat.mod_lt _ hg
      have hb : t % g < g := Nat.mod_lt _ hg
      cases hu : used (t / g % g) (t % g) with
      | false => exact Or.inl rfl
      | true => exact Or.inr (hpair _ ha _ hb hu)

private theorem prod_binary (s : Finset ι) (f : ι → ℤ) (hf : ∀ i ∈ s, f i = 0 ∨ f i = 1) :
    (∏ i ∈ s, f i) = 0 ∨ (∏ i ∈ s, f i) = 1 := by
  classical
  by_cases h : ∀ i ∈ s, f i = 1
  · exact Or.inr (Finset.prod_eq_one h)
  · push_neg at h
    obtain ⟨i, hi, hne⟩ := h
    rcases hf i hi with hz | ho
    · exact Or.inl (Finset.prod_eq_zero hi hz)
    · exact False.elim (hne ho)

theorem weightCheck_eq (g K : ℕ) (hg : 0 < g) (χ : Fin g → Fin K) (γ : Fin K)
    (N : Net ℤ g) (z : (N.coarsen χ).S γ) (member node : ℕ → Bool) (used pair : ℕ → ℕ → Bool)
    (hm : ∀ a : Fin g, member a.val = decide (χ a = γ))
    (hn : ∀ (a : Fin g) (h : χ a = γ), (if node a.val then (1 : ℤ) else 0) = N.W a (z ⟨a, h⟩))
    (hu : ∀ a b : Fin g, used a.val b.val = decide (a < b ∧ χ a = γ ∧ χ b = γ))
    (hp : ∀ (a b : Fin g) (ha : χ a = γ) (hb : χ b = γ),
      (if pair a.val b.val then (1 : ℤ) else 0) = N.Ψ a b (z ⟨a, ha⟩) (z ⟨b, hb⟩)) :
    weightCheck g member node used pair = true ↔ (N.coarsen χ).W γ z = 1 := by
  rw [weightCheck_true_iff g hg, Net.coarsen_W_code]
  let f : Fin g → ℤ := fun a => if h : χ a = γ then N.W a (z ⟨a, h⟩) else 1
  let p : Fin g × Fin g → ℤ := fun ab =>
    if h : χ ab.1 = γ ∧ χ ab.2 = γ then N.Ψ ab.1 ab.2 (z ⟨ab.1, h.1⟩) (z ⟨ab.2, h.2⟩) else 1
  change _ ↔ (∏ a, f a) * (∏ ab ∈ pairs g, p ab) = 1
  have hf : ∀ a ∈ (Finset.univ : Finset (Fin g)), f a = 0 ∨ f a = 1 := by
    intro a _
    by_cases h : χ a = γ
    · have hw := hn a h
      cases hnode : node a.val <;> simp_all [f, h]
    · simp [f, h]
  have hp' : ∀ ab ∈ pairs g, p ab = 0 ∨ p ab = 1 := by
    intro ab _
    by_cases h : χ ab.1 = γ ∧ χ ab.2 = γ
    · have hw := hp ab.1 ab.2 h.1 h.2
      cases hpair : pair ab.1.val ab.2.val <;> simp_all [p, h]
    · simp [p, h]
  have hfb := prod_binary (Finset.univ : Finset (Fin g)) f hf
  have hpb := prod_binary (pairs g) p hp'
  have hmul : ((∏ a, f a) * (∏ ab ∈ pairs g, p ab) = 1) ↔
      ((∏ a, f a) = 1 ∧ (∏ ab ∈ pairs g, p ab) = 1) := by
    rcases hfb with h0 | h1 <;> rcases hpb with h2 | h3 <;> simp_all
  rw [hmul]
  constructor
  · rintro ⟨hnode, hpair⟩
    constructor
    · apply Finset.prod_eq_one
      intro a _
      by_cases h : χ a = γ
      · have hm' : member a.val = true := by simpa [h] using hm a
        have hv := hnode a.val a.isLt hm'
        simp [f, h, ← hn a h, hv]
      · simp [f, h]
    · apply Finset.prod_eq_one
      intro ab hab
      have hlt := mem_pairs.mp hab
      by_cases h : χ ab.1 = γ ∧ χ ab.2 = γ
      · have hu' : used ab.1.val ab.2.val = true := by simpa [hlt, h] using hu ab.1 ab.2
        have hv := hpair ab.1.val ab.1.isLt ab.2.val ab.2.isLt hu'
        simp [p, h, ← hp ab.1 ab.2 h.1 h.2, hv]
      · simp [p, h]
  · rintro ⟨hnode, hpair⟩
    constructor
    · intro i hi hm'
      let a : Fin g := ⟨i, hi⟩
      have h : χ a = γ := by
        have hd : decide (χ a = γ) = true := (hm a).symm.trans (by simpa [a] using hm')
        exact of_decide_eq_true hd
      have hfa : f a = 1 := by
        rcases hf a (Finset.mem_univ a) with hz | ho
        · have hzprod : (∏ a, f a) = 0 := Finset.prod_eq_zero (Finset.mem_univ a) hz
          omega
        · exact ho
      have heq := hn a h
      have hw : N.W a (z ⟨a, h⟩) = 1 := by simpa [f, h] using hfa
      have : (if node i then (1 : ℤ) else 0) = 1 := by simpa [a] using heq.trans hw
      cases hni : node i <;> simp_all
    · intro i hi j hj hu'
      let a : Fin g := ⟨i, hi⟩
      let b : Fin g := ⟨j, hj⟩
      have hcond : a < b ∧ χ a = γ ∧ χ b = γ := by
        have hd : decide (a < b ∧ χ a = γ ∧ χ b = γ) = true :=
          (hu a b).symm.trans (by simpa [a, b] using hu')
        exact of_decide_eq_true hd
      have hab : (a, b) ∈ pairs g := mem_pairs.mpr hcond.1
      have hpab : p (a, b) = 1 := by
        rcases hp' (a, b) hab with hz | ho
        · have hzprod : (∏ ab ∈ pairs g, p ab) = 0 := Finset.prod_eq_zero hab hz
          omega
        · exact ho
      have heq := hp a b hcond.2.1 hcond.2.2
      have hw : N.Ψ a b (z ⟨a, hcond.2.1⟩) (z ⟨b, hcond.2.2⟩) = 1 := by
        simpa [p, hcond] using hpab
      have : (if pair i j then (1 : ℤ) else 0) = 1 := by simpa [a, b] using heq.trans hw
      cases hpij : pair i j <;> simp_all

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.weightCheck_polyT
#print axioms IntegerFactoring.HalfAppendix.weightCheck_true_iff
#print axioms IntegerFactoring.HalfAppendix.weightCheck_eq
