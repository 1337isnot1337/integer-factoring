/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time capped products for mixed-radix prefix decoding. This is an
arithmetic component; it does not yet implement the factoring state decoder.
-/
import TopParamsPoly
import FactoringAux.EngineCode

open Groundwork

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def cappedProduct (B : ℕ) (f : ℕ → ℕ) (n : ℕ) : ℕ :=
  (List.range n).foldl (fun acc i => min (acc * f i) B) (min 1 B)

theorem cappedFold_le (B : ℕ) (f : ℕ → ℕ) : ∀ (ks : List ℕ) (a : ℕ), a ≤ B →
    ks.foldl (fun acc i => min (acc * f i) B) a ≤ B
  | [], a, ha => ha
  | i :: ks, a, _ => cappedFold_le B f ks _ (Nat.min_le_right _ _)

theorem cappedProduct_polyT {α : Type} [Enc α] {B n : α → ℕ} {f : α → ℕ → ℕ}
    (hB : PolyT B) (hn : Unary n) (hf : PolyT fun p : α × ℕ => f p.1 p.2) :
    PolyT fun a => cappedProduct (B a) (f a) (n a) := by
  have hlist : PolyT fun a => List.range (n a) := PolyT.range hn
  have hb' : PolyT fun p : α × ℕ × ℕ => B p.1 := PolyT.comp hB PolyT.id'.fst
  have hf' : PolyT fun p : α × ℕ × ℕ => f p.1 p.2.2 :=
    PolyT.comp hf (PolyT.prodMk PolyT.id'.fst PolyT.id'.snd.snd)
  have hstep : PolyT fun p : α × ℕ × ℕ => min (p.2.1 * f p.1 p.2.2) (B p.1) := by fun_prop
  have hinit : PolyT fun a => min 1 (B a) := by fun_prop
  obtain ⟨_, t, ht, _, _, hlen⟩ := hB.elim
  have hfold : PolyT fun a => (List.range (n a)).foldl (fun acc i => min (acc * f a i) (B a)) (min 1 (B a)) := by
    refine PolyT.foldl hlist hinit hstep ?_
    refine ⟨t, ht, ?_⟩
    intro a k
    let ks := (List.range (n a)).take k
    have hb := cappedFold_le (B a) (f a) ks (min 1 (B a)) (Nat.min_le_right _ _)
    have hs := Nat.size_le_size hb
    have hl := hlen a
    change len (ks.foldl _ _) ≤ t (len a)
    rw [len_nat] at hl ⊢
    omega
  exact hfold

theorem min_mul_cap (B a b : ℕ) (hb : 1 ≤ b) : min (min a B * b) B = min (a * b) B := by
  by_cases h : a ≤ B
  · simp [Nat.min_eq_left h]
  · have h1 : B ≤ a * b := by nlinarith
    have h2 : B ≤ B * b := by nlinarith
    simp [Nat.min_eq_right (Nat.le_of_not_ge h), Nat.min_eq_right h1, Nat.min_eq_right h2]

theorem cappedFold_eq_min (B : ℕ) (f : ℕ → ℕ) : ∀ (ks : List ℕ) (a : ℕ),
    (∀ i ∈ ks, 1 ≤ f i) →
    ks.foldl (fun acc i => min (acc * f i) B) (min a B) = min (a * (ks.map f).prod) B
  | [], a, _ => by simp
  | i :: ks, a, h => by
      simp only [List.foldl_cons, List.map_cons, List.prod_cons]
      rw [min_mul_cap B a (f i) (h i (by simp))]
      rw [cappedFold_eq_min B f ks (a * f i) (by intro j hj; exact h j (by simp [hj]))]
      congr 1
      ring

theorem cappedProduct_eq (B n : ℕ) (f : ℕ → ℕ)
    (hpos : ∀ i < n, 1 ≤ f i) (hbound : ((List.range n).map f).prod ≤ B) :
    cappedProduct B f n = ((List.range n).map f).prod := by
  unfold cappedProduct
  rw [cappedFold_eq_min B f (List.range n) 1 (by simpa using hpos)]
  simpa [Nat.min_eq_left hbound]

theorem range_prod_eq_pp (f : ℕ → ℕ) : ∀ n, ((List.range n).map f).prod = Factoring.Aux.pp f n
  | 0 => rfl
  | n + 1 => by
      rw [List.range_succ, List.map_append, List.prod_append, range_prod_eq_pp f n]
      simp [Factoring.Aux.pp]

def cappedDigit (B g code i : ℕ) (radix : ℕ → ℕ) : ℕ :=
  code / cappedProduct B radix (min i g) % radix i

theorem cappedDigit_polyT {α : Type} [Enc α] {B code : α → ℕ} {g : α → ℕ} {radix : α → ℕ → ℕ}
    (hB : PolyT B) (hg : Unary g) (hc : PolyT code)
    (hr : PolyT fun p : α × ℕ => radix p.1 p.2) :
    PolyT fun p : α × ℕ => cappedDigit (B p.1) (g p.1) (code p.1) p.2 (radix p.1) := by
  have hg' : Unary fun p : α × ℕ => g p.1 := Unary.comp hg PolyT.id'.fst
  have hidx : PolyT fun p : α × ℕ => min p.2 (g p.1) := PolyT.min' PolyT.id'.snd (Unary.poly hg')
  have hidxu : Unary fun p : α × ℕ => min p.2 (g p.1) := Unary.of_le hg' hidx (fun p => Nat.min_le_right _ _)
  have hB' : PolyT fun p : α × ℕ => B p.1 := PolyT.comp hB PolyT.id'.fst
  have hr' : PolyT fun p : (α × ℕ) × ℕ => radix p.1.1 p.2 :=
    PolyT.comp (f := fun q : α × ℕ => radix q.1 q.2)
      (g := fun p : (α × ℕ) × ℕ => (p.1.1, p.2)) hr (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hprod : PolyT fun p : α × ℕ => cappedProduct (B p.1) (radix p.1) (min p.2 (g p.1)) :=
    cappedProduct_polyT (B := fun p : α × ℕ => B p.1)
      (n := fun p : α × ℕ => min p.2 (g p.1))
      (f := fun p : α × ℕ => radix p.1) hB' hidxu hr'
  have hcode : PolyT fun p : α × ℕ => code p.1 := PolyT.comp hc PolyT.id'.fst
  exact PolyT.mod (PolyT.div hcode hprod) hr

theorem cappedDigit_eq (B g code i : ℕ) (radix : ℕ → ℕ) (hi : i ≤ g)
    (hpos : ∀ j < i, 1 ≤ radix j) (hbound : Factoring.Aux.pp radix i ≤ B) :
    cappedDigit B g code i radix = Factoring.Aux.dgx radix code i := by
  unfold cappedDigit Factoring.Aux.dgx
  rw [Nat.min_eq_left hi, cappedProduct_eq B i radix hpos]
  · rw [range_prod_eq_pp]
  · rw [range_prod_eq_pp]; exact hbound

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.cappedProduct_polyT
#print axioms IntegerFactoring.HalfAppendix.cappedProduct_eq
#print axioms IntegerFactoring.HalfAppendix.cappedDigit_polyT
#print axioms IntegerFactoring.HalfAppendix.cappedDigit_eq
