/-
SPDX-License-Identifier: Apache-2.0

Charge a polynomial-time Appendix fragment against an exponential input
length. This is not a runtime proof for constructing that input or for the
complete factoring machine.
-/
import WrittenPoly

open Groundwork Groundwork.TM

namespace IntegerFactoring.HalfAppendix

theorem polyT_fragment_exp {α β : Type} [Enc α] [Enc β] (f : α → β)
    (hf : PolyT f) (A : ℕ) :
    ∃ (P : Fr) (D : ℕ), ∀ (x : α) (q : ℕ), 1 ≤ q → len x ≤ 2 ^ (A * q) →
      ∀ (M S T N : List Sym), ∃ t ≤ 2 ^ (D * q),
        P.Runs none (st (wd x ++ M) S T N) t (st (wd (f x) ++ M) S T N) := by
  obtain ⟨P, c, hrun, _⟩ := hf
  refine ⟨P, c * (A + 2), ?_⟩
  intro x q hq hlen M S T N
  obtain ⟨t, ht, hrun'⟩ := hrun x M S T N
  have hone : 1 ≤ 2 ^ (A * q) := Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
  have hlen' : len x + 1 ≤ 2 ^ ((A + 1) * q) := by
    calc
      len x + 1 ≤ 2 ^ (A * q) + 1 := by omega
      _ ≤ 2 ^ (A * q) + 2 ^ (A * q) := by omega
      _ = 2 ^ (A * q + 1) := by rw [pow_add]; norm_num; ring
      _ ≤ 2 ^ ((A + 1) * q) := Nat.pow_le_pow_right (by norm_num) (by nlinarith)
  have hc : c ≤ 2 ^ (c * q) := by
    exact Nat.lt_two_pow_self.le.trans (Nat.pow_le_pow_right (by norm_num) (by nlinarith))
  have htime : c * (len x + 1) ^ c ≤ 2 ^ (c * (A + 2) * q) := by
    calc
      c * (len x + 1) ^ c ≤ 2 ^ (c * q) * (2 ^ ((A + 1) * q)) ^ c :=
        Nat.mul_le_mul hc (Nat.pow_le_pow_left hlen' c)
      _ = 2 ^ (c * q + ((A + 1) * q) * c) := by rw [← pow_mul, ← pow_add]
      _ = 2 ^ (c * (A + 2) * q) := congrArg (fun n : ℕ => 2 ^ n) (by ring)
  exact ⟨t, ht.trans htime, hrun'⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.polyT_fragment_exp
