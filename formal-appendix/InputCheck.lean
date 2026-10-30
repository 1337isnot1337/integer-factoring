/-
SPDX-License-Identifier: Apache-2.0

Numeric input-gate weight check. Gate-message transforms and the complete row
summand remain separate obligations.
-/
import GateFields

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def inputCheck (A h L B k y r code : ℕ) : Bool :=
  let base := 2 ^ L
  let bitsOK := (List.range h).all fun i => (List.range L).all fun t =>
    let pos := L * (r + A * i) + t
    decide (k ≤ pos → (fixedRadixDigit B h base code i).testBit t = y.testBit (pos - k))
  bitsOK && decide (r = 0 → y < 2 ^ (L * L - k))

theorem inputCheck_polyT {α : Type} [Enc α] {A h L B : α → ℕ} {k y r code : α → ℕ}
    (hA : Unary A) (hh : Unary h) (hL : Unary L) (hB : Unary B)
    (hk : PolyT k) (hy : PolyT y) (hr : PolyT r) (hc : PolyT code) :
    PolyT fun a => inputCheck (A a) (h a) (L a) (B a) (k a) (y a) (r a) (code a) := by
  let β := (α × ℕ) × ℕ
  have hbase := PolyT.two_pow hL
  have hd : PolyT fun p : α × ℕ => fixedRadixDigit (B p.1) (h p.1) (2 ^ L p.1) (code p.1) p.2 :=
    fixedRadixDigit_polyT (B := B) (count := h) (base := fun a => 2 ^ L a) (code := code) hB hh hbase hc
  have hd' : PolyT fun p : β => fixedRadixDigit (B p.1.1) (h p.1.1) (2 ^ L p.1.1) (code p.1.1) p.1.2 :=
    PolyT.comp hd PolyT.id'.fst
  have hA' : PolyT fun p : β => A p.1.1 := PolyT.comp hA.poly PolyT.id'.fst.fst
  have hL' : PolyT fun p : β => L p.1.1 := PolyT.comp hL.poly PolyT.id'.fst.fst
  have hk' : PolyT fun p : β => k p.1.1 := PolyT.comp hk PolyT.id'.fst.fst
  have hy' : PolyT fun p : β => y p.1.1 := PolyT.comp hy PolyT.id'.fst.fst
  have hr' : PolyT fun p : β => r p.1.1 := PolyT.comp hr PolyT.id'.fst.fst
  have hi' : PolyT fun p : β => p.1.2 := PolyT.id'.fst.snd
  have ht' : PolyT fun p : β => p.2 := PolyT.id'.snd
  have hpos : PolyT fun p : β => L p.1.1 * (r p.1.1 + A p.1.1 * p.1.2) + p.2 :=
    PolyT.add (PolyT.mul hL' (PolyT.add hr' (PolyT.mul hA' hi'))) ht'
  have hleft := PolyT.testBit hd' ht'
  have hright := PolyT.testBit hy' (PolyT.sub hpos hk')
  have himp := PolyT.or' (PolyT.not' (PolyT.decide_le hk' hpos)) (PolyT.decide_eq hleft hright)
  have hpred : PolyT fun p : β => decide (k p.1.1 ≤ L p.1.1 * (r p.1.1 + A p.1.1 * p.1.2) + p.2 →
      (fixedRadixDigit (B p.1.1) (h p.1.1) (2 ^ L p.1.1) (code p.1.1) p.1.2).testBit p.2 =
        (y p.1.1).testBit (L p.1.1 * (r p.1.1 + A p.1.1 * p.1.2) + p.2 - k p.1.1)) := by
    exact himp.congr fun p => by simp
  have hinnerRange : PolyT fun p : α × ℕ => List.range (L p.1) := PolyT.range (Unary.comp hL PolyT.id'.fst)
  have hinner := PolyT.all (f := fun p t => decide (k p.1 ≤ L p.1 * (r p.1 + A p.1 * p.2) + t →
    (fixedRadixDigit (B p.1) (h p.1) (2 ^ L p.1) (code p.1) p.2).testBit t =
      (y p.1).testBit (L p.1 * (r p.1 + A p.1 * p.2) + t - k p.1))) hpred hinnerRange
  have hbits := PolyT.all (f := fun a i => (List.range (L a)).all fun t =>
    decide (k a ≤ L a * (r a + A a * i) + t →
      (fixedRadixDigit (B a) (h a) (2 ^ L a) (code a) i).testBit t =
        (y a).testBit (L a * (r a + A a * i) + t - k a))) hinner (PolyT.range hh)
  have hL2 := Unary.mul hL hL
  have hexpp : PolyT fun a => L a * L a - k a := PolyT.sub hL2.poly hk
  have hexp : Unary fun a => L a * L a - k a := Unary.of_le hL2 hexpp (fun a => Nat.sub_le _ _)
  have hbound := PolyT.decide_lt hy (PolyT.two_pow hexp)
  have hguard := PolyT.or' (PolyT.not' (PolyT.decide_eq hr (PolyT.const 0))) hbound
  have hfinal := PolyT.and' hbits hguard
  exact hfinal.congr fun a => by simp [inputCheck]

theorem inputCheck_eq (s B N k y : ℕ) (r : Fin (pA s)) (w : Fin (ph s) → Fin (2 ^ pL s))
    (hsize : pD s ≤ 2 ^ B) :
    inputCheck (pA s) (ph s) (pL s) B k y r.val (dCode s w).val =
      decide (mW s N k y (.faP r) w = (1 : ℤ)) := by
  have hb : ∀ (i : ℕ) (hi : i < ph s), fixedRadixDigit B (ph s) (2 ^ pL s) (dCode s w).val i = (w ⟨i, hi⟩).val := by
    intro i hi
    exact fixed_input_digit_of_size s B w ⟨i, hi⟩ hsize
  apply Bool.eq_iff_iff.mpr
  simp only [inputCheck, Bool.and_eq_true, List.all_eq_true, List.mem_range, decide_eq_true_eq]
  simp only [mW]
  rw [ite_eq_left_iff]
  norm_num
  unfold preOK
  intro _
  constructor
  · intro h
    intro i t hk
    have h' := h i.val i.isLt t.val t.isLt hk
    simpa [hb i.val i.isLt] using h'
  · intro h
    intro i hi t ht hk
    have h' := h ⟨i, hi⟩ ⟨t, ht⟩ hk
    simpa [hb i hi] using h'

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.inputCheck_polyT
#print axioms IntegerFactoring.HalfAppendix.inputCheck_eq
