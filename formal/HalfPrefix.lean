/-
SPDX-License-Identifier: Apache-2.0

The Appendix B top-network value counts all divisors. The factoring search
in the pinned reference counts only divisors at least two. This lemma
accounts for the sole excluded divisor; it does not construct either machine.
-/
import FactoringAux.Search
import Mathlib

namespace IntegerFactoring.HalfPrefix

open Finset Factoring.Aux

theorem all_divisors_count (N k x : ℕ) (hN : 2 ≤ N) :
    ((N.divisors.filter fun d => d / 2 ^ k = x).card) =
      cntPre N k x + if 1 / 2 ^ k = x then 1 else 0 := by
  let s := N.divisors.filter fun d => d / 2 ^ k = x
  have hlarge : (s.filter fun d => 2 ≤ d).card = cntPre N k x := by
    simp [s, cntPre, Finset.filter_filter, and_comm]
  have hsmall : s.filter (fun d => ¬ 2 ≤ d) =
      if 1 / 2 ^ k = x then {1} else ∅ := by
    ext d
    by_cases h : 1 / 2 ^ k = x
    · simp only [h, ite_true, mem_filter, mem_singleton, s]
      constructor
      · rintro ⟨⟨hd, _⟩, hlt⟩
        have hp := Nat.pos_of_mem_divisors hd
        omega
      · intro hd
        subst d
        simp [show N ≠ 0 by omega, h]
    · simp only [h, ite_false, mem_filter, Finset.notMem_empty, iff_false, s]
      rintro ⟨⟨hd, hprefix⟩, hlt⟩
      have hp := Nat.pos_of_mem_divisors hd
      have : d = 1 := by omega
      subst d
      exact h hprefix
  have hpart := Finset.card_filter_add_card_filter_not (s := s) (fun d => 2 ≤ d)
  rw [hlarge, hsmall] at hpart
  by_cases h : 1 / 2 ^ k = x <;> simp [h, s] at hpart ⊢ <;> omega

theorem count_without_one (N k x : ℕ) (hN : 2 ≤ N) :
    cntPre N k x = (N.divisors.filter fun d => d / 2 ^ k = x).card -
      (if 1 / 2 ^ k = x then 1 else 0) := by
  have h := all_divisors_count N k x hN
  omega

end IntegerFactoring.HalfPrefix

#print axioms IntegerFactoring.HalfPrefix.all_divisors_count
#print axioms IntegerFactoring.HalfPrefix.count_without_one
