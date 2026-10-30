/-
SPDX-License-Identifier: Apache-2.0

Numeric field extraction for the reference multiplication-gate encoding.
The factorial gate checks and complete clique-row summand remain open.
-/
import CappedProducts
import GateSize

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def fixedRadixDigit (B count base code i : ℕ) : ℕ :=
  cappedDigit (2 ^ B) count code i (fun _ => base)

theorem fixedRadixDigit_polyT {α : Type} [Enc α] {B count : α → ℕ} {base code : α → ℕ}
    (hB : Unary B) (hcount : Unary count) (hb : PolyT base) (hc : PolyT code) :
    PolyT fun p : α × ℕ => fixedRadixDigit (B p.1) (count p.1) (base p.1) (code p.1) p.2 := by
  have hr : PolyT fun p : α × ℕ => base p.1 := PolyT.comp hb PolyT.id'.fst
  exact cappedDigit_polyT (B := fun a => 2 ^ B a) (g := count) (code := code)
    (radix := fun a _ => base a) (PolyT.two_pow hB) hcount hc hr

theorem fixedRadixDigit_eq (B count base code i : ℕ) (hi : i ≤ count) (hb : 1 ≤ base)
    (hbound : pp (fun _ => base) i ≤ 2 ^ B) :
    fixedRadixDigit B count base code i = code / pp (fun _ => base) i % base := by
  exact cappedDigit_eq (2 ^ B) count code i (fun _ => base) hi (fun _ _ => hb) hbound

theorem fixed_input_digit (s B : ℕ) (w : Fin (ph s) → Fin (2 ^ pL s)) (i : Fin (ph s))
    (hbound : pp (fun _ => 2 ^ pL s) i.val ≤ 2 ^ B) :
    fixedRadixDigit B (ph s) (2 ^ pL s) (dCode s w).val i.val = (w i).val := by
  have hb : 1 ≤ 2 ^ pL s := Nat.one_le_two_pow
  rw [fixedRadixDigit_eq B (ph s) (2 ^ pL s) (dCode s w).val i.val (Nat.le_of_lt i.isLt) hb hbound]
  exact (dCode_dig s w i).symm

theorem fixed_residue_digit (s B : ℕ) (z : Fin (pA s) → ZMod (pM s)) (i : Fin (pA s))
    (hbound : pp (fun _ => pM s) i.val ≤ 2 ^ B) :
    fixedRadixDigit B (pA s) (pM s) (zCode s z).val i.val = (z i).val := by
  have hb : 1 ≤ pM s := Nat.le_add_left 1 _
  rw [fixedRadixDigit_eq B (pA s) (pM s) (zCode s z).val i.val (Nat.le_of_lt i.isLt) hb hbound]
  exact (zCode_dig s z i).symm

theorem constant_pp_prefix_le (base i n : ℕ) (hb : 1 ≤ base) (hi : i ≤ n) :
    pp (fun _ => base) i ≤ pp (fun _ => base) n := by
  apply monotone_nat_of_le_succ (f := pp (fun _ => base)) ?_ hi
  intro j
  change pp (fun _ => base) j ≤ pp (fun _ => base) j * base
  exact Nat.le_mul_of_pos_right _ hb

theorem fixed_input_digit_of_size (s B : ℕ) (w : Fin (ph s) → Fin (2 ^ pL s)) (i : Fin (ph s))
    (hsize : pD s ≤ 2 ^ B) :
    fixedRadixDigit B (ph s) (2 ^ pL s) (dCode s w).val i.val = (w i).val := by
  apply fixed_input_digit
  exact (constant_pp_prefix_le (2 ^ pL s) i.val (ph s) Nat.one_le_two_pow (Nat.le_of_lt i.isLt)).trans hsize

theorem fixed_residue_digit_of_size (s B : ℕ) (z : Fin (pA s) → ZMod (pM s)) (i : Fin (pA s))
    (hsize : pZ s ≤ 2 ^ B) :
    fixedRadixDigit B (pA s) (pM s) (zCode s z).val i.val = (z i).val := by
  apply fixed_residue_digit
  exact (constant_pp_prefix_le (pM s) i.val (pA s) (Nat.le_add_left 1 _) (Nat.le_of_lt i.isLt)).trans hsize

def midLeft (Z code : ℕ) : ℕ := code / Z % Z
def midRight (Z code : ℕ) : ℕ := code % Z
def carryResidue (M F code : ℕ) : ℕ := code / 2 ^ F % M
def carryIncoming (F code : ℕ) : ℕ := code % 2 ^ F

theorem gateSimpleFields_polyT {α : Type} [Enc α] {n code : α → ℕ}
    (hn : Unary n) (hc : PolyT code) :
    PolyT (fun a => (midLeft (pZ (sOf (n a))) (code a), midRight (pZ (sOf (n a))) (code a),
      carryResidue (pM (sOf (n a))) (fw (sOf (n a))) (code a), carryIncoming (fw (sOf (n a))) (code a))) := by
  have hz := selected_pZ_polyT hn
  have hm := selected_pM_polyT hn
  have hf := PolyT.two_pow (selected_fw_unary hn)
  have hleft := PolyT.mod (PolyT.div hc hz) hz
  have hright := PolyT.mod hc hz
  have hcarry := PolyT.mod (PolyT.div hc hf) hm
  have hin := PolyT.mod hc hf
  exact PolyT.prodMk hleft (PolyT.prodMk hright (PolyT.prodMk hcarry hin))

theorem gateSimpleFields_eq (s : ℕ) (k : Fin (pA s)) (mid : mS s (.mid k))
    (j : Fin (pA s * pA s)) (carry : mS s (.carry j)) :
    midLeft (pZ s) (gcode s (.mid k) mid).val = (zCode s mid.1).val ∧
    midRight (pZ s) (gcode s (.mid k) mid).val = (zCode s mid.2).val ∧
    carryResidue (pM s) (fw s) (gcode s (.carry j) carry).val = carry.1.val ∧
    carryIncoming (fw s) (gcode s (.carry j) carry).val = carry.2.val := by
  exact ⟨mid_dig1 s k mid, mid_dig2 s k mid, carry_y s j carry, carry_c s j carry⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.fixedRadixDigit_polyT
#print axioms IntegerFactoring.HalfAppendix.fixed_input_digit
#print axioms IntegerFactoring.HalfAppendix.fixed_residue_digit
#print axioms IntegerFactoring.HalfAppendix.fixed_input_digit_of_size
#print axioms IntegerFactoring.HalfAppendix.fixed_residue_digit_of_size
#print axioms IntegerFactoring.HalfAppendix.gateSimpleFields_polyT
#print axioms IntegerFactoring.HalfAppendix.gateSimpleFields_eq
