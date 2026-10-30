/-
SPDX-License-Identifier: Apache-2.0

Numeric coding of the reference's little-endian Boolean words.
-/
import SelectedRadix
import StepMachines.Math.CutWord

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem wordCode_polyT {α : Type} [Enc α] {w : α → List Bool} (hw : PolyT w) :
    PolyT fun a => PolyT.bval (w a) := PolyT.bval' hw

theorem testBit_wordCode (w : List Bool) (i : ℕ) : (PolyT.bval w).testBit i = w.getD i false := by
  induction w generalizing i with
  | nil => simp [PolyT.bval]
  | cons b w ih =>
    cases i with
    | zero => simp [PolyT.bval]
    | succ i =>
      change (Nat.bit b (PolyT.bval w)).testBit (i + 1) = w.getD i false
      rw [Nat.testBit_bit_succ, ih]

theorem bitsOf_wordCode (L : ℕ) (w : List Bool) :
    bdec L (PolyT.bval w) = Groundwork.Grid.bitsOf L w := by
  funext i
  exact testBit_wordCode w i.val

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.wordCode_polyT
#print axioms IntegerFactoring.HalfAppendix.bitsOf_wordCode
