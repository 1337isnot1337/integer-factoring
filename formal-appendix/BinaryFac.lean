/-
SPDX-License-Identifier: Apache-2.0

Zero-one bounds for the reference's factorization operations.
-/
import FactoringAux.Gates
import Mathlib

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def Is01 (z : ℤ) : Prop := z = 0 ∨ z = 1

theorem is01_mul {x y : ℤ} (hx : Is01 x) (hy : Is01 y) : Is01 (x * y) := by
  rcases hx with rfl | rfl <;> rcases hy with rfl | rfl <;> simp [Is01]

def Fac01 {X Y : Type} {c : ℕ} (f : Fac ℤ X Y c) : Prop :=
  (∀ x t, Is01 (f.l x t)) ∧ (∀ y t, Is01 (f.r y t))

theorem fac_cast_01 {X Y : Type} {c d : ℕ} (h : c = d) (f : Fac ℤ X Y c)
    (hf : Fac01 f) : Fac01 (f.cast h) := by
  cases h
  exact hf

theorem fac_mul_01 {X Y : Type} {c d : ℕ} (f : Fac ℤ X Y c) (g : Fac ℤ X Y d)
    (hf : Fac01 f) (hg : Fac01 g) : Fac01 (f.mul g) := by
  constructor
  · intro x t; exact is01_mul (hf.1 _ _) (hg.1 _ _)
  · intro y t; exact is01_mul (hf.2 _ _) (hg.2 _ _)

theorem fac_prodList_01 {X Y ι : Type} (c : ι → ℕ) (F : ∀ i, Fac ℤ X Y (c i))
    (hF : ∀ i, Fac01 (F i)) (xs : List ι) : Fac01 (Fac.prodList c F xs) := by
  induction xs with
  | nil =>
      apply fac_cast_01
      simp [Fac01, Fac.one, Is01]
  | cons i xs ih =>
      apply fac_cast_01
      exact fac_mul_01 (F i) (Fac.prodList c F xs) (hF i) ih

theorem fac_pad_01 {X Y : Type} {c : ℕ} (f : Fac ℤ X Y c) (L : ℕ) (h : c ≤ L)
    (hf : Fac01 f) : Fac01 (f.pad L h) := by
  unfold Fac.pad
  apply fac_cast_01
  apply fac_mul_01 f (Fac.unit (L - c)) hf
  constructor <;> intro x t <;> simp [Fac.unit, Is01] <;> tauto

theorem numTest_01 {X Y : Type} (c : ℕ) (v : X → ℕ) (w : Y → ℕ) :
    Fac01 (Fac.numTest c v w : Fac ℤ X Y c) := by
  constructor <;> intro x t <;> simp [Fac.numTest, Is01] <;> tauto

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.fac_prodList_01
#print axioms IntegerFactoring.HalfAppendix.fac_pad_01
#print axioms IntegerFactoring.HalfAppendix.numTest_01
