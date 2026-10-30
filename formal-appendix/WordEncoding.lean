/-
SPDX-License-Identifier: Apache-2.0

The reference challenge's binary encoding agrees with the executable word code.
-/
import WordCode
import Challenge

open Groundwork Computability

namespace IntegerFactoring.HalfAppendix

theorem encodePosNum_eq_bits (p : PosNum) : encodePosNum p = Nat.bits (p : ℕ) := by
  induction p with
  | one => simp [encodePosNum, Nat.one_bits]
  | bit0 p ih =>
    rw [encodePosNum, PosNum.cast_bit0]
    rw [show (p : ℕ) + p = 2 * (p : ℕ) by omega, Nat.bit0_bits _ (ne_of_gt (PosNum.to_nat_pos p)), ih]
  | bit1 p ih =>
    rw [encodePosNum, PosNum.cast_bit1]
    rw [show (p : ℕ) + p + 1 = 2 * (p : ℕ) + 1 by omega, Nat.bit1_bits, ih]

theorem encodeNat_eq_bits (n : ℕ) : encodeNat n = Nat.bits n := by
  have h := Num.to_of_nat n
  cases hn : (n : Num) with
  | zero =>
    have : n = 0 := by simpa [hn] using h.symm
    subst n
    rfl
  | pos p =>
    have he : (p : ℕ) = n := by simpa [hn] using h
    unfold encodeNat
    rw [hn]
    exact (encodePosNum_eq_bits p).trans (congrArg Nat.bits he)

theorem wordCode_encodeNat (n : ℕ) : PolyT.bval (encodeNat n) = n := by
  rw [encodeNat_eq_bits, PolyT.bval_bits]

theorem encodeNat_length (n : ℕ) : (encodeNat n).length = n.size := by
  rw [encodeNat_eq_bits, Nat.size_eq_bits_len]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.encodeNat_eq_bits
#print axioms IntegerFactoring.HalfAppendix.wordCode_encodeNat
#print axioms IntegerFactoring.HalfAppendix.encodeNat_length
