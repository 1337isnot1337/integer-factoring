/-
SPDX-License-Identifier: Apache-2.0

Zero-one semantics of the reference multiplication network. This is a
precursor to bounding the entries after grouping.
-/
import BinaryFac
import FactoringAux.MulCount

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

structure Net01 {g : ℕ} (N : Net ℤ g) : Prop where
  weights : ∀ a z, Is01 (N.W a z)
  factors : ∀ a b, Fac01 (N.F a b)
  interactions : ∀ a b, a < b → ∀ z z', Is01 (N.Ψ a b z z')

theorem multiplicationNet_01 (s N k y : ℕ) : Net01 (mulNet s ℤ N k y) := by
  constructor
  · intro a z
    change Is01 (mW s N k y ((MG.enum (pA s)).symm a) z)
    exact mW_01 s N k y _ z
  · intro a b
    exact numTest_01 _ _ _
  · intro a b hab z z'
    have hm (i j : MG (pA s)) (u : mS s i) (v : mS s j) :
        Is01 ((Fac.numTest (mC s i j) (mV s i j) (mV s j i) : Fac ℤ _ _ _).eval u v) := by
      rw [Fac.eval_numTest]
      split <;> simp [Is01]
    exact hm ((MG.enum (pA s)).symm a) ((MG.enum (pA s)).symm b) z z'

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.multiplicationNet_01
