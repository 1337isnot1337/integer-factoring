/-
SPDX-License-Identifier: Apache-2.0

A sufficient padding exponent for scanning every top-grouping code. The
predicate and the rest of the executable construction are not certified here.
-/
import SelectedPadCost

open Levels.Gen IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

/-- The number of grouping codes is bounded by an explicit power of two. -/
theorem grouping_count_le_pad (P : ℕ) :
    P ^ (3 * P ^ 8) ≤ 2 ^ (3 * P ^ 9) := by
  calc
    P ^ (3 * P ^ 8) ≤ (2 ^ P) ^ (3 * P ^ 8) :=
      Nat.pow_le_pow_left (Nat.lt_two_pow_self.le) _
    _ = 2 ^ (3 * P ^ 9) := by rw [← pow_mul]; congr 1; ring

/-- Under the selected level's hypotheses, the extra exhaustive-search pad
    has exponent at most a fixed multiple of n/r. -/
theorem grouping_pad_cost (P n r : ℕ) (hP : 2 ≤ P)
    (hn : P ^ 40 ≤ n) (hr : r ≤ P) :
    r * (3 * P ^ 9) ≤ 3 * n := by
  calc
    r * (3 * P ^ 9) ≤ P * (3 * P ^ 9) := Nat.mul_le_mul_right _ hr
    _ = 3 * P ^ 10 := by ring
    _ ≤ 3 * P ^ 40 := Nat.mul_le_mul_left _ (Nat.pow_le_pow_right (by omega : 0 < P) (by omega : 10 ≤ 40))
    _ ≤ 3 * n := Nat.mul_le_mul_left _ hn

theorem scheduled_grouping_pad_cost (x l n : ℕ) (hx : 1 ≤ x)
    (hn : (K x l) ^ 40 ≤ n) :
    2 ^ (l / 2) * (3 * (K x l) ^ 9) ≤ 3 * n :=
  grouping_pad_cost (K x l) n (2 ^ (l / 2)) (K_two x l hx) hn (half_level_le_K x l hx)

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.grouping_count_le_pad
#print axioms IntegerFactoring.HalfAppendix.scheduled_grouping_pad_cost
