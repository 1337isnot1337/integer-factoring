/-
SPDX-License-Identifier: Apache-2.0

Modular positivity for a scaled count: the multiplier need not be invertible.
This uses the upstream small-moduli theorem and does not construct the evaluator
or prove its multiplier and count bounds.
-/
module

public import FactoringAux.Moduli
public import Mathlib.Data.Nat.ModEq

@[expose] public section

namespace Factoring.Improvements

open Factoring.Aux

/-- The explicit bit budget in the candidate supplies the smallness hypothesis
of the upstream modular zero test. -/
theorem modular_budget_small (M V : ℕ) (hM : 4 ≤ M)
    (hV : V ≤ 2 ^ (M - 4)) : (2 * M + 1) * V < 4 ^ M := by
  have htwo : ∀ n : ℕ, n + 1 ≤ 2 ^ n := by
    intro n
    induction n with
    | zero => simp
    | succ n ih =>
        rw [pow_succ]
        omega
  have hfactor : 2 * M + 1 ≤ 2 ^ (M + 1) := by
    have h := Nat.mul_le_mul_left 2 (htwo M)
    rw [pow_succ]
    omega
  calc
    (2 * M + 1) * V ≤ 2 ^ (M + 1) * 2 ^ (M - 4) := Nat.mul_le_mul hfactor hV
    _ = 2 ^ (M + 1 + (M - 4)) := (pow_add 2 (M + 1) (M - 4)).symm
    _ < 2 ^ (M + M) := Nat.pow_lt_pow_right (by norm_num) (by omega)
    _ = 4 ^ M := by
      have he : M + M = 2 * M := by omega
      rw [he, pow_mul]
      norm_num

/-- Removing a known baseline from a scaled count can be checked by comparing
residues, including when the multiplier is not invertible modulo the modulus. -/
theorem scaled_positive_iff_residue_difference (H q b M : ℕ)
    (hH : 0 < H) (hM : 0 < M)
    (hsmall : (2 * M + 1) * (H * q) < 4 ^ M) :
    (∃ m, 2 ≤ m ∧ m ≤ 2 * M ∧ (H * (q + b)) % m ≠ (H * b) % m) ↔ 0 < q := by
  constructor
  · rintro ⟨m, _, _, hne⟩
    by_contra hq
    have : q = 0 := Nat.eq_zero_of_not_pos hq
    simp [this] at hne
  · intro hq
    by_contra hnone
    push Not at hnone
    have hdvd : ∀ m, 1 ≤ m → m ≤ 2 * M → m ∣ H * q := by
      intro m hm hmM
      by_cases hm1 : m = 1
      · simp [hm1]
      · have hm2 : 2 ≤ m := by omega
        have heq := hnone m hm2 hmM
        have hcongr : H * q + H * b ≡ 0 + H * b [MOD m] := by
          show (H * q + H * b) % m = (0 + H * b) % m
          simpa [Nat.mul_add] using heq
        have hzero : H * q ≡ 0 [MOD m] :=
          Nat.ModEq.add_right_cancel (Nat.ModEq.refl (H * b)) hcongr
        exact Nat.dvd_of_mod_eq_zero (by simpa [Nat.ModEq] using hzero)
    have hz := eq_zero_of_all_dvd hM hdvd hsmall
    have hp : 0 < H * q := Nat.mul_pos hH hq
    omega

/-- A multiplier and count with known bit budgets admit a finite residue test.
The bound on the multiplier is an assumption to be supplied by the evaluator. -/
theorem scaled_positive_from_bit_budgets (H q b n S : ℕ)
    (hH : 0 < H) (hHbound : H ≤ 2 ^ S) (hqbound : q ≤ 2 ^ n) :
    (∃ m, 2 ≤ m ∧ m ≤ 2 * (n + S + 4) ∧
      (H * (q + b)) % m ≠ (H * b) % m) ↔ 0 < q := by
  let M := n + S + 4
  have hM : 4 ≤ M := by dsimp [M]; omega
  have hV : H * q ≤ 2 ^ (M - 4) := by
    calc
      H * q ≤ 2 ^ S * 2 ^ n := Nat.mul_le_mul hHbound hqbound
      _ = 2 ^ (M - 4) := by
        rw [← pow_add]
        congr 1
        dsimp [M]
        omega
  exact scaled_positive_iff_residue_difference H q b M hH (by omega)
    (modular_budget_small M (H * q) hM hV)

#print axioms scaled_positive_iff_residue_difference
#print axioms modular_budget_small
#print axioms scaled_positive_from_bit_budgets

end Factoring.Improvements
