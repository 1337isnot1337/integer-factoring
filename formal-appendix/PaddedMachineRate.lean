/-
SPDX-License-Identifier: Apache-2.0

Exponentiate the conditional finite-alphabet machine's polynomial overhead.
The concrete factoring algorithm remains a hypothesis.
-/
import PaddedWordMachine

open Groundwork PaddedDescent

namespace IntegerFactoring.HalfAppendix

theorem padded_machine_time_exp (q A n t E W : ℕ) (ht : 1 ≤ t)
    (hn : n ≤ 2 ^ (A * t)) (hE : E ≤ A * t) (hW : W ≤ 2 ^ (A * t)) :
    q * (E + 2) * (3 * n + 2 * 2 ^ E + W + 5) ^ q ≤
      2 ^ ((q + A + 2 + (A + 4) * q) * t) := by
  let H := 2 ^ (A * t)
  have hH : 1 ≤ H := Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
  have hpow : 2 ^ E ≤ H := Nat.pow_le_pow_right (by norm_num) hE
  have hE' : E ≤ H := Nat.lt_two_pow_self.le.trans hpow
  have harg : 3 * n + 2 * 2 ^ E + W + 5 ≤ 16 * H := by omega
  have he2 : E + 2 ≤ 4 * H := by omega
  have hq : q ≤ 2 ^ (q * t) :=
    Nat.lt_two_pow_self.le.trans (Nat.pow_le_pow_right (by norm_num) (by nlinarith))
  have h4 : 4 * H ≤ 2 ^ ((A + 2) * t) := by
    calc
      4 * H = 2 ^ (A * t + 2) := by dsimp [H]; rw [pow_add]; norm_num; ring
      _ ≤ 2 ^ ((A + 2) * t) := Nat.pow_le_pow_right (by norm_num) (by nlinarith)
  have h16 : 16 * H ≤ 2 ^ ((A + 4) * t) := by
    calc
      16 * H = 2 ^ (A * t + 4) := by dsimp [H]; rw [pow_add]; norm_num; ring
      _ ≤ 2 ^ ((A + 4) * t) := Nat.pow_le_pow_right (by norm_num) (by nlinarith)
  calc
    q * (E + 2) * (3 * n + 2 * 2 ^ E + W + 5) ^ q ≤
        2 ^ (q * t) * 2 ^ ((A + 2) * t) * (2 ^ ((A + 4) * t)) ^ q :=
      Nat.mul_le_mul (Nat.mul_le_mul hq (he2.trans h4)) (Nat.pow_le_pow_left (harg.trans h16) q)
    _ = 2 ^ ((q + A + 2 + (A + 4) * q) * t) := by
      rw [← pow_mul, ← pow_add, ← pow_add]
      congr 1
      ring

/-- The reference factoring-in-time statement follows once a correct padded
    polynomial-time algorithm and the stated schedule bounds are supplied. -/
theorem padded_factoring_exp {e : List Bool → ℕ} {F : List Bool × List Unit → List Bool}
    (A : ℕ) (t : ℕ → ℕ) (he : PolyT e) (hF : PolyT F)
    (ht : ∀ n, 1 ≤ t n) (hn : ∀ n, n ≤ 2 ^ (A * t n))
    (hE : ∀ N, e (Computability.encodeNat N) ≤ A * t (Computability.encodeNat N).length)
    (hW : ∀ N, size (F (Computability.encodeNat N, pad (e (Computability.encodeNat N)))) ≤
      2 ^ (A * t (Computability.encodeNat N).length))
    (hc : ∀ N, F (Computability.encodeNat N, pad (e (Computability.encodeNat N))) =
      Computability.encodeNat (Factoring.factor N)) :
    ∃ D, Factoring.FactoringInTime (fun n => 2 ^ (D * t n)) := by
  obtain ⟨M, q, hfin, htime⟩ := padded_factoring_machine
    (fun n => A * t n) (fun n => 2 ^ (A * t n)) he hF hE hW hc
  refine ⟨q + A + 2 + (A + 4) * q, M, hfin, ?_⟩
  intro n
  exact (htime n).trans (padded_machine_time_exp q A n (t n) (A * t n)
    (2 ^ (A * t n)) (ht n) (hn n) le_rfl le_rfl)

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.padded_machine_time_exp
#print axioms IntegerFactoring.HalfAppendix.padded_factoring_exp
