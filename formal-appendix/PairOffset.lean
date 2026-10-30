/-
SPDX-License-Identifier: Apache-2.0

Uniform numeric prefix offsets for coarsened-network message blocks.
-/
import CrossQueryCorrect
import ToTuring.Poly.NetFun

open Groundwork Factoring.Aux Finset

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def pairOffset (g k : ℕ) (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool) : ℕ :=
  ((List.range (min k (g * g))).map fun t =>
    if used (t / g % g) (t % g) then width (t / g % g) (t % g) else 0).sum

theorem pairOffset_unary {α : Type} [Enc α] {g W k : α → ℕ}
    {width : α → ℕ → ℕ → ℕ} {used : α → ℕ → ℕ → Bool}
    (hg : Unary g) (hW : Unary W) (hk : PolyT k)
    (hw : PolyT fun p : (α × ℕ) × ℕ => width p.1.1 p.1.2 p.2)
    (hu : PolyT fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
    (hbound : ∀ a i j, width a i j ≤ W a) :
    Unary fun a => pairOffset (g a) (k a) (width a) (used a) := by
  have hg2 := Unary.mul hg hg
  have hn : Unary fun a => min (k a) (g a * g a) :=
    Unary.of_le hg2 (PolyT.min' hk hg2.poly) (fun _ => Nat.min_le_right _ _)
  have hlist : PolyT fun a => List.range (min (k a) (g a * g a)) := PolyT.range hn
  have hg' : PolyT fun p : α × ℕ => g p.1 := PolyT.comp hg.poly PolyT.id'.fst
  have ha : PolyT fun p : α × ℕ => p.2 / g p.1 % g p.1 := PolyT.mod (PolyT.div PolyT.id'.snd hg') hg'
  have hb : PolyT fun p : α × ℕ => p.2 % g p.1 := PolyT.mod PolyT.id'.snd hg'
  have hpair : PolyT fun p : α × ℕ => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1) :=
    PolyT.prodMk (PolyT.prodMk PolyT.id'.fst ha) hb
  have hval : PolyT fun p : α × ℕ => width p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) :=
    PolyT.comp (f := fun p : (α × ℕ) × ℕ => width p.1.1 p.1.2 p.2)
      (g := fun p : α × ℕ => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1)) hw hpair
  have hused : PolyT fun p : α × ℕ => used p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) :=
    PolyT.comp (f := fun p : (α × ℕ) × ℕ => used p.1.1 p.1.2 p.2)
      (g := fun p : α × ℕ => ((p.1, p.2 / g p.1 % g p.1), p.2 % g p.1)) hu hpair
  have hterm : PolyT fun p : α × ℕ =>
      if used p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) then width p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) else 0 :=
    (PolyT.cond' hused hval (PolyT.const 0)).congr (fun p => by cases used p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) <;> rfl)
  have hWu : Unary fun p : α × ℕ => W p.1 := Unary.comp hW PolyT.id'.fst
  have htermU : Unary fun p : α × ℕ =>
      if used p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) then width p.1 (p.2 / g p.1 % g p.1) (p.2 % g p.1) else 0 := by
    apply Unary.of_le hWu hterm
    intro p
    split <;> simp_all [hbound]
  exact Groundwork.Prog.NetE.unary_sum_map (f := fun a t =>
    if used a (t / g a % g a) (t % g a) then width a (t / g a % g a) (t % g a) else 0) htermU hlist

theorem pairOffset_eq (g K k : ℕ) (hg : 0 < g) (hk : k ≤ g * g)
    (χ : Fin g → Fin K) (u v : Fin K) (N : Net ℤ g) (width : ℕ → ℕ → ℕ) (used : ℕ → ℕ → Bool)
    (hw : ∀ a b : Fin g, width a.val b.val = N.c a b)
    (hu : ∀ a b : Fin g, used a.val b.val = decide (Net.UsedP χ u v (a, b))) :
    pairOffset g k width used = offR (fun p => N.c p.1 p.2) (Net.UsedP χ u v) (Net.pairOf g hg) k := by
  unfold pairOffset offR
  rw [Nat.min_eq_left hk, ← List.sum_toFinset _ List.nodup_range, List.toFinset_range]
  apply Finset.sum_congr rfl
  intro t ht
  let a : Fin g := ⟨t / g % g, Nat.mod_lt _ hg⟩
  let b : Fin g := ⟨t % g, Nat.mod_lt _ hg⟩
  have hused := hu a b
  have hwidth := hw a b
  change (if used a.val b.val then width a.val b.val else 0) = _
  rw [hused, hwidth]
  simp only [decide_eq_true_eq]
  rfl

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.pairOffset_unary
#print axioms IntegerFactoring.HalfAppendix.pairOffset_eq
