/-
SPDX-License-Identifier: Apache-2.0

An integer rank-seven factorization of a signed 8 by 8 complement block.
Eight disjoint XOR blocks give the 56-channel construction in
research/crossing-compression.md. This file proves the block factorization
and the scalar crossing identity; it does not formalize the full recursion.
-/
module

public import FactoringAux.Cross
public import Mathlib.Tactic.Ring
public import Mathlib.Tactic.Linarith

@[expose] public section

namespace Factoring.Improvements

open Finset Factoring.Aux

def complementLeft (s : Fin 8 → ℤ) (x : Fin 8) (a : Fin 7) : ℤ :=
  (if x = a.succ then 1 else 0) - (if x = 0 then s a.succ else 0)

def complementRight (s : Fin 8 → ℤ) (a : Fin 7) (v : Fin 8) : ℤ :=
  8 * (if a.succ = v then 1 else 0) - s a.succ * s v

/-- All coefficients are integers. The seven channels are indexed by `Fin 7`. -/
theorem signed_complement_factorization (s : Fin 8 → ℤ) (h0 : s 0 = 1)
    (hs : ∀ a, s a * s a = 1) (x v : Fin 8) :
    ∑ a : Fin 7, complementLeft s x a * complementRight s a v =
      8 * (if x = v then 1 else 0) - s x * s v := by
  have hzero (a : Fin 7) : (0 : Fin 8) ≠ a.succ := (Fin.succ_ne_zero a).symm
  refine Fin.cases ?_ (fun x => ?_) x
  · refine Fin.cases ?_ (fun v => ?_) v
    · simp [complementLeft, complementRight, hzero, Fin.succ_ne_zero, h0, neg_mul, hs]
    · have hterm (a : Fin 7) :
          (-s a.succ) * (8 * (if a = v then 1 else 0) - s a.succ * s v.succ) =
            -8 * (if a = v then s a.succ else 0) + s v.succ := by
        have h : s a.succ * s a.succ * s v.succ = s v.succ := by
          rw [hs, one_mul]
        split_ifs <;> nlinarith
      simp only [complementLeft, complementRight, hzero,
        Fin.succ_inj, ite_false, ite_true, zero_sub]
      simp_rw [hterm]
      simp [Finset.sum_add_distrib, h0]
      ring
  · simp [complementLeft, complementRight, Fin.succ_ne_zero, Fin.succ_inj,
      ite_mul, eq_comm]

/-- The new residual is an integer linear combination of the original two
crossing tensors. No inverse is introduced in any coefficient ring. -/
theorem crossing_complement_identity (x y u v : Blk) :
    8 * (if x = v ∧ y = u then (1 : ℤ) else 0) - Rx x y u v =
      3 * Rx x y u v + 4 * Cx x y u v := by
  have h := crossing x y u v
  split_ifs at h ⊢ <;> linarith

/-- Exact rational envelope used by the recurrence certificate. -/
theorem seven_pow_26_lt_two_pow_73 : (7 : ℕ) ^ 26 < 2 ^ 73 := by
  norm_num

#print axioms signed_complement_factorization
#print axioms crossing_complement_identity
#print axioms seven_pow_26_lt_two_pow_73

end Factoring.Improvements
