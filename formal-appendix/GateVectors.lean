/-
SPDX-License-Identifier: Apache-2.0

Numeric input-digit and residue-vector accessors for multiplication gates.
-/
import GateFields

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def inputVector (s B code i : ℕ) : ℕ :=
  if i < ph s then fixedRadixDigit B (ph s) (2 ^ pL s) code i else 0

def residueVector (s B code i : ℕ) : ℕ :=
  if i < pA s then fixedRadixDigit B (pA s) (pM s) code i else 0

theorem inputVector_polyT {α : Type} [Enc α] {n B code : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hc : PolyT code) :
    PolyT fun p : α × ℕ => inputVector (sOf (n p.1)) (B p.1) (code p.1) p.2 := by
  have hh : Unary fun a => ph (sOf (n a)) := selected_S_unary hn
  have hL := selected_pL_unary hn
  have hd : PolyT fun p : α × ℕ => fixedRadixDigit (B p.1) (ph (sOf (n p.1)))
      (2 ^ pL (sOf (n p.1))) (code p.1) p.2 :=
    fixedRadixDigit_polyT (B := B) (count := fun a => ph (sOf (n a)))
      (base := fun a => 2 ^ pL (sOf (n a))) (code := code) hB hh (PolyT.two_pow hL) hc
  have hh' : PolyT fun p : α × ℕ => ph (sOf (n p.1)) := PolyT.comp (f := fun a => ph (sOf (n a))) hh.poly PolyT.id'.fst
  exact PolyT.ite' (PolyT.decide_lt PolyT.id'.snd hh') hd (PolyT.const 0)

theorem residueVector_polyT {α : Type} [Enc α] {n B code : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hc : PolyT code) :
    PolyT fun p : α × ℕ => residueVector (sOf (n p.1)) (B p.1) (code p.1) p.2 := by
  have hA := selected_pA_unary hn
  have hm := selected_pM_polyT hn
  have hd : PolyT fun p : α × ℕ => fixedRadixDigit (B p.1) (pA (sOf (n p.1)))
      (pM (sOf (n p.1))) (code p.1) p.2 :=
    fixedRadixDigit_polyT (B := B) (count := fun a => pA (sOf (n a)))
      (base := fun a => pM (sOf (n a))) (code := code) hB hA hm hc
  have hA' : PolyT fun p : α × ℕ => pA (sOf (n p.1)) := PolyT.comp (f := fun a => pA (sOf (n a))) hA.poly PolyT.id'.fst
  exact PolyT.ite' (PolyT.decide_lt PolyT.id'.snd hA') hd (PolyT.const 0)

theorem inputVector_eq (s B : ℕ) (w : Fin (ph s) → Fin (2 ^ pL s)) (hsize : pD s ≤ 2 ^ B) (i : ℕ) :
    (inputVector s B (dCode s w).val i : ZMod (pM s)) = digZ s w i := by
  by_cases hi : i < ph s
  · simp [inputVector, digZ, ext0, hi, fixed_input_digit_of_size s B w ⟨i, hi⟩ hsize]
  · simp [inputVector, digZ, ext0, hi]

theorem residueVector_eq (s B : ℕ) (z : Fin (pA s) → ZMod (pM s)) (hsize : pZ s ≤ 2 ^ B) (i : ℕ) :
    (residueVector s B (zCode s z).val i : ZMod (pM s)) = ext0 z i := by
  by_cases hi : i < pA s
  · simp [residueVector, ext0, hi, fixed_residue_digit_of_size s B z ⟨i, hi⟩ hsize]
  · simp [residueVector, ext0, hi]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.inputVector_polyT
#print axioms IntegerFactoring.HalfAppendix.residueVector_polyT
#print axioms IntegerFactoring.HalfAppendix.inputVector_eq
#print axioms IntegerFactoring.HalfAppendix.residueVector_eq
