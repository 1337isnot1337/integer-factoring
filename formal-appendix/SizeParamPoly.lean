/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time computation of the reference multiplication-network size
parameter from a unary input width.
-/
import LevelPoly
import FactoringAux.EngineTop

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def cappedSizeStep (n s : ℕ) : ℕ :=
  if 2 ^ (4 * min s n + 2) < n then s + 1 else s

def sizeFold (n j : ℕ) : ℕ :=
  (List.range j).foldl (fun s _ => cappedSizeStep n s) 0

theorem iterSize_le (n : ℕ) : ∀ j, iterF (sStep n) 0 j ≤ j
  | 0 => le_rfl
  | j + 1 => by
      have ih := iterSize_le n j
      simp only [iterF, sStep]
      split <;> omega

theorem sizeFold_eq (n : ℕ) : ∀ j, j ≤ n → sizeFold n j = iterF (sStep n) 0 j
  | 0, _ => rfl
  | j + 1, hj => by
      have ih := sizeFold_eq n j (by omega)
      have hs := iterSize_le n j
      unfold sizeFold at ih ⊢
      rw [List.range_succ, List.foldl_append]
      simp only [List.foldl_cons, List.foldl_nil, iterF]
      rw [ih]
      simp [cappedSizeStep, sStep, Nat.min_eq_left (show iterF (sStep n) 0 j ≤ n by omega)]

theorem sizeFold_bound (n : ℕ) : ∀ (ks : List ℕ) (s : ℕ),
    ks.foldl (fun s _ => cappedSizeStep n s) s ≤ s + ks.length
  | [], s => by simp
  | k :: ks, s => by
      have h := sizeFold_bound n ks (cappedSizeStep n s)
      have hs : cappedSizeStep n s ≤ s + 1 := by unfold cappedSizeStep; split <;> omega
      simp only [List.foldl_cons, List.length_cons]
      omega

theorem sOf_polyT {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    PolyT fun a => sOf (n a) := by
  have hnp : PolyT n := Unary.poly hn
  have hlist : PolyT fun a => List.range (n a) := PolyT.range hn
  have hcap : Unary fun p : α × ℕ => min p.2 (n p.1) := by
    have hb : Unary fun p : α × ℕ => n p.1 := Unary.comp hn (by fun_prop)
    have hp : PolyT fun p : α × ℕ => min p.2 (n p.1) := by fun_prop
    exact Unary.of_le hb hp (fun p => Nat.min_le_right _ _)
  have hexp : Unary fun p : α × ℕ => 4 * min p.2 (n p.1) + 2 := by fun_prop
  have hpow : PolyT fun p : α × ℕ => 2 ^ (4 * min p.2 (n p.1) + 2) := PolyT.two_pow hexp
  have hstep : PolyT fun p : α × ℕ × ℕ => cappedSizeStep (n p.1) p.2.1 := by
    unfold cappedSizeStep
    fun_prop
  obtain ⟨_, t, ht, _, _, hlen⟩ := hn.elim
  have hf : PolyT fun a => (List.range (n a)).foldl (fun s _ => cappedSizeStep (n a) s) 0 := by
    refine PolyT.foldl hlist (PolyT.const 0) hstep ?_
    refine ⟨t, ht, ?_⟩
    intro a k
    let ks := (List.range (n a)).take k
    have hb := sizeFold_bound (n a) ks 0
    have hks : ks.length ≤ n a := by simp [ks]
    have hs : Nat.size (ks.foldl (fun s _ => cappedSizeStep (n a) s) 0) ≤ n a := by
      exact (Nat.size_le.2 Nat.lt_two_pow_self).trans (by omega)
    have hlen' := hlen a
    rw [len_units, List.length_replicate] at hlen'
    change len (ks.foldl _ _) ≤ t (len a)
    rw [len_nat]
    omega
  apply hf.congr
  intro a
  exact sizeFold_eq (n a) (n a) le_rfl

theorem sOf_unary {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    Unary fun a => sOf (n a) := by
  exact Unary.of_le hn (sOf_polyT hn) (fun a => iterSize_le (n a) (n a))

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.sOf_polyT
#print axioms IntegerFactoring.HalfAppendix.sOf_unary
