/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time exhaustive grouping-code search when its search range is
provided in unary and the grouping predicate has a polynomial-time test.
-/
import FactoringAux.EngineGroup
import PaddedDescent.F.CertNet
import ToTuring.Lib.List2
import ToTuring.Lib.Nat2

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def firstGoodFold (ok : ℕ → Bool) (S n : ℕ) : ℕ :=
  (List.range n).foldl (fun s x => if s = S && ok x then x else s) S

theorem firstGoodFold_eq (ok : ℕ → Bool) (S : ℕ) : ∀ n,
    firstGoodFold ok S n = firstGood (fun x => ok x = true) S n
  | 0 => rfl
  | n + 1 => by
      have ih := firstGoodFold_eq ok S n
      unfold firstGoodFold at ih ⊢
      rw [List.range_succ, List.foldl_append]
      simp only [List.foldl_cons, List.foldl_nil]
      change (if (decide (_ = S) && ok n) then n else _) =
        fgStep (fun x => ok x = true) S n (firstGood (fun x => ok x = true) S n)
      rw [ih]
      simp [fgStep]

theorem firstGoodFold_bound (ok : ℕ → Bool) (S : ℕ) : ∀ (xs : List ℕ) (s : ℕ),
    s ≤ S → (∀ x ∈ xs, x ≤ S) →
      xs.foldl (fun s x => if s = S && ok x then x else s) s ≤ S
  | [], s, hs, _ => hs
  | x :: xs, s, hs, hx => by
      apply firstGoodFold_bound ok S xs _ ?_ (by intro y hy; exact hx y (by simp [hy]))
      dsimp
      split <;> simp_all

theorem firstGood_polyT {α : Type} [Enc α] {S : α → ℕ} {ok : α → ℕ → Bool}
    (hS : Unary S) (hok : PolyT fun p : α × ℕ => ok p.1 p.2) :
    PolyT fun x => firstGood (fun n => ok x n = true) (S x) (S x) := by
  have hSp : PolyT S := Unary.poly hS
  have hlist : PolyT fun x => List.range (S x) := PolyT.range hS
  have hstep : PolyT fun p : α × ℕ × ℕ =>
      if p.2.1 = S p.1 && ok p.1 p.2.2 then p.2.2 else p.2.1 := by fun_prop
  obtain ⟨_, t, ht, _, _, hlen⟩ := hS.elim
  have hfold : PolyT fun x => (List.range (S x)).foldl
      (fun s n => if s = S x && ok x n then n else s) (S x) := by
    refine PolyT.foldl hlist ?_ hstep ?_
    · exact hSp
    · refine ⟨t, ht, ?_⟩
      intro x k
      let xs := (List.range (S x)).take k
      have hb := firstGoodFold_bound (ok x) (S x) xs (S x) (le_refl _) (by
        intro z hz
        have : z ∈ List.range (S x) := List.mem_of_mem_take hz
        exact Nat.le_of_lt (List.mem_range.mp this))
      have hs := Nat.size_le_size hb
      have hsize : Nat.size (S x) ≤ S x := Nat.size_le.2 Nat.lt_two_pow_self
      have hlen' := hlen x
      rw [len_units, List.length_replicate] at hlen'
      change len (xs.foldl _ _) ≤ t (len x)
      rw [len_nat]
      omega
  apply hfold.congr
  intro x
  exact firstGoodFold_eq (ok x) (S x) (S x)

/-- A supplied unit pad can be truncated to an exponentially sized search
    range in polynomial time when its exponent is unary. -/
theorem bounded_exp_pad_polyT {α : Type} [Enc α] {e : α → ℕ}
    (he : Unary e) : PolyT fun p : α × List Unit => p.2.take (2 ^ e p.1) := by
  have he' : Unary fun p : α × List Unit => e p.1 := Unary.comp he (by fun_prop)
  exact PolyT.take (PolyT.two_pow he') (by fun_prop)

theorem bounded_exp_pad_eq (e : ℕ) (pad : List Unit) (h : 2 ^ e ≤ pad.length) :
    pad.take (2 ^ e) = List.replicate (2 ^ e) () := by
  have hp : pad = List.replicate pad.length () := List.eq_replicate_length.mpr (by intro x _; cases x; rfl)
  conv_lhs => rw [hp]
  simp [h]

theorem bounded_exp_unary {α : Type} [Enc α] {e : α → ℕ} (he : Unary e) :
    Unary fun p : α × List Unit => min (2 ^ e p.1) p.2.length := by
  have h := Unary.of_length (bounded_exp_pad_polyT he)
  exact h.congr (fun p => by simp [List.length_take, Nat.min_comm])

theorem padded_firstGood_polyT {α : Type} [Enc α] {e : α → ℕ}
    {ok : (α × List Unit) → ℕ → Bool} (he : Unary e)
    (hok : PolyT fun p : (α × List Unit) × ℕ => ok p.1 p.2) :
    PolyT fun p : α × List Unit =>
      firstGood (fun n => ok p n = true) (min (2 ^ e p.1) p.2.length)
        (min (2 ^ e p.1) p.2.length) :=
  firstGood_polyT (bounded_exp_unary he) hok

/-- The candidate count can be any polynomial-time binary number: the unit
    pad supplies the unary budget for the exhaustive scan. -/
theorem padded_firstGood_of_count_polyT {α : Type} [Enc α] {S : α → ℕ}
    {ok : (α × List Unit) → ℕ → Bool} (hS : PolyT S)
    (hok : PolyT fun p : (α × List Unit) × ℕ => ok p.1 p.2) :
    PolyT fun p : α × List Unit =>
      firstGood (fun n => ok p n = true) (min (S p.1) p.2.length)
        (min (S p.1) p.2.length) := by
  have hc : PolyT fun p : α × List Unit => S p.1 := PolyT.comp hS (by fun_prop)
  have hp : PolyT fun p : α × List Unit => p.2.take (S p.1) := PolyT.take hc (by fun_prop)
  have hu : Unary fun p : α × List Unit => min (S p.1) p.2.length := by
    exact (Unary.of_length hp).congr (fun p => by simp [List.length_take, Nat.min_comm])
  exact firstGood_polyT hu hok

theorem padded_firstGood_of_count_eq {α : Type} (S : α → ℕ)
    (ok : (α × List Unit) → ℕ → Bool) (p : α × List Unit)
    (h : S p.1 ≤ p.2.length) :
    firstGood (fun n => ok p n = true) (min (S p.1) p.2.length)
        (min (S p.1) p.2.length) =
      firstGood (fun n => ok p n = true) (S p.1) (S p.1) := by
  rw [Nat.min_eq_left h]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.firstGood_polyT
#print axioms IntegerFactoring.HalfAppendix.bounded_exp_pad_polyT
#print axioms IntegerFactoring.HalfAppendix.bounded_exp_pad_eq
#print axioms IntegerFactoring.HalfAppendix.padded_firstGood_polyT
#print axioms IntegerFactoring.HalfAppendix.padded_firstGood_of_count_polyT
#print axioms IntegerFactoring.HalfAppendix.padded_firstGood_of_count_eq
