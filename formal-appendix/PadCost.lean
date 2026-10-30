/-
SPDX-License-Identifier: Apache-2.0

Arithmetic cost envelope for the explicit evaluator pad. The hypotheses
must still be established for the selected schedule and top constructor.
-/
import HalfPad

open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem pad_exponent_mul_le (P L B C r n A D : ℕ)
    (hrP : r ≤ P) (hP4 : P ^ 4 ≤ n)
    (hwidth : P ^ 2 * L ≤ A * n)
    (hlog : Nat.log 2 (P * (L + B)) + 1 ≤ D * P) :
    r * padExponent P L B C r ≤ (A + 1 + C * D) * n := by
  have hrn : r ≤ n := hrP.trans ((Nat.le_self_pow (by norm_num : 4 ≠ 0) P).trans hP4)
  have hquot : r * (P ^ 2 * L / r) ≤ A * n :=
    (Nat.mul_div_le (P ^ 2 * L) r).trans hwidth
  have hterm : r * (C * P ^ 2 * (Nat.log 2 (P * (L + B)) + 1)) ≤ (C * D) * n := by
    calc
      r * (C * P ^ 2 * (Nat.log 2 (P * (L + B)) + 1)) ≤ r * (C * P ^ 2 * (D * P)) := by gcongr
      _ = C * D * r * P ^ 3 := by ring
      _ = (C * D) * (r * P ^ 3) := by ring
      _ ≤ (C * D) * (P * P ^ 3) := Nat.mul_le_mul_left _ (Nat.mul_le_mul_right _ hrP)
      _ = C * D * P ^ 4 := by ring
      _ ≤ C * D * n := by gcongr
  calc
    r * padExponent P L B C r = r * (P ^ 2 * L / r) + r +
        r * (C * P ^ 2 * (Nat.log 2 (P * (L + B)) + 1)) := by simp [padExponent]; ring
    _ ≤ A * n + n + (C * D) * n := Nat.add_le_add (Nat.add_le_add hquot hrn) hterm
    _ = (A + 1 + C * D) * n := by ring

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.pad_exponent_mul_le
