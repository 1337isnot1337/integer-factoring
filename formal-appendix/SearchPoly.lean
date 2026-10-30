/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time prefix search conditional on a polynomial-time count query
and unary input width. It does not certify the factoring-specific query.
-/
import FactoringAux.Search
import PaddedDescent.F.CertNet
import ToTuring.Lib.List2

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def searchFold (c : ℕ → ℕ → ℕ) (n p : ℕ) : ℕ :=
  (List.range n).reverse.foldl (fun p k => searchStep c k p) p

theorem searchFold_eq (c : ℕ → ℕ → ℕ) : ∀ n p, searchFold c n p = searchGo c n p
  | 0, p => rfl
  | n + 1, p => by
      simp only [searchFold, List.range_succ, List.reverse_append, List.reverse_singleton, List.singleton_append,
        List.foldl_cons, searchGo]
      exact searchFold_eq c n (searchStep c n p)

theorem fold_prefix_bound (c : ℕ → ℕ → ℕ) :
    ∀ (ks : List ℕ) (p m : ℕ), p < 2 ^ m →
      ks.foldl (fun p k => searchStep c k p) p < 2 ^ (m + ks.length)
  | [], p, m, h => by simpa using h
  | k :: ks, p, m, h => by
      have hs : searchStep c k p < 2 ^ (m + 1) := by
        unfold searchStep
        rw [pow_succ]
        split_ifs <;> omega
      have ih := fold_prefix_bound c ks (searchStep c k p) (m + 1) hs
      simpa [List.foldl_cons, Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using ih

/-- A uniform query fragment yields a uniform search fragment when the input
width is available in unary. The query fragment itself is a hypothesis. -/
theorem search_polyT {α : Type} [Enc α] {width : α → ℕ} {c : α → ℕ → ℕ → ℕ}
    (hw : Unary width) (hc : PolyT fun p : α × ℕ × ℕ => c p.1 p.2.1 p.2.2) :
    PolyT fun x => search (c x) (width x) := by
  have hlist : PolyT fun x => (List.range (width x)).reverse := PolyT.reverse (PolyT.range hw)
  have hstep : PolyT fun p : α × ℕ × ℕ => searchStep (c p.1) p.2.2 p.2.1 := by
    unfold searchStep
    fun_prop
  obtain ⟨_, t, ht, _, _, hlen⟩ := hw.elim
  have hfold : PolyT fun x => ((List.range (width x)).reverse).foldl
      (fun p k => searchStep (c x) k p) 0 := by
    refine PolyT.foldl (l := fun x => (List.range (width x)).reverse) (a := fun _ => 0)
      (g := fun x p k => searchStep (c x) k p) hlist (PolyT.const 0) hstep ?_
    refine ⟨t, ht, ?_⟩
    intro x k
    let ks := ((List.range (width x)).reverse).take k
    have hs : ks.length ≤ width x := by simp [ks]
    have hb := fold_prefix_bound (c x) ks 0 0 (by norm_num)
    have hb' : ks.foldl (fun p k => searchStep (c x) k p) 0 < 2 ^ width x := by
      exact hb.trans_le (Nat.pow_le_pow_right (by norm_num) (by simpa using hs))
    have hsize : Nat.size (ks.foldl (fun p k => searchStep (c x) k p) 0) ≤ width x := Nat.size_le.2 hb'
    have hlen' := hlen x
    rw [len_units, List.length_replicate] at hlen'
    change len (ks.foldl (fun p k => searchStep (c x) k p) 0) ≤ t (len x)
    rw [len_nat]
    omega
  apply hfold.congr
  intro x
  exact (searchFold_eq (c x) (width x) 0).trans rfl

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.searchFold_eq
#print axioms IntegerFactoring.HalfAppendix.search_polyT
