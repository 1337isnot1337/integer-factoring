/-
SPDX-License-Identifier: Apache-2.0

An explicit finite search for the maximal Appendix B level. This file proves
the selection inequalities, not a machine-time bound for evaluating it.
-/
import HalfSchedule

open Levels.Gen

namespace IntegerFactoring.HalfAppendix

def selectedLevel (x n : ℕ) : ℕ :=
  Nat.findGreatest (fun l => (K x l) ^ 40 ≤ n) (Nat.log 2 n + 1)

theorem selectedLevel_spec (x n : ℕ) (hx : 1 ≤ x) (h0 : (K x 0) ^ 40 ≤ n) :
    (K x (selectedLevel x n)) ^ 40 ≤ n ∧ n < (K x (selectedLevel x n + 1)) ^ 40 := by
  let l := selectedLevel x n
  have hselected : (K x l) ^ 40 ≤ n := by
    exact Nat.findGreatest_spec (P := fun d => (K x d) ^ 40 ≤ n) (Nat.zero_le _) h0
  have hlpow : 2 ^ l ≤ n :=
    (two_pow_le_K x l hx).trans ((Nat.le_self_pow (by norm_num : 40 ≠ 0) _).trans hselected)
  have hlt : l < Nat.log 2 n + 1 := by
    by_contra h
    have hpow := Nat.pow_le_pow_right (by norm_num : 0 < 2)
      (show Nat.log 2 n + 1 ≤ l by omega)
    have hbig := Nat.lt_pow_succ_log_self (by norm_num : 1 < 2) n
    omega
  have hnext : n < (K x (l + 1)) ^ 40 := by
    by_contra h
    have hc : (K x (l + 1)) ^ 40 ≤ n := by omega
    have hmax : l + 1 ≤ l := by
      exact Nat.le_findGreatest (show l + 1 ≤ Nat.log 2 n + 1 by omega) hc
    omega
  exact ⟨hselected, hnext⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedLevel_spec
