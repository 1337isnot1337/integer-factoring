/-
SPDX-License-Identifier: Apache-2.0

Numeric checks for the two sides of a multiplication-gate factor.
-/
import GatePair

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

private theorem bnum_bdec_remainder (c q : ℕ) : bnum (bdec c q) = q % 2 ^ c := by
  have e : bdec c (q % 2 ^ c) = bdec c q := by
    funext i
    simp [bdec, Nat.testBit_mod_two_pow, i.2]
  rw [← e, bnum_bdec (Nat.mod_lt _ (by positivity))]

def gateSide (s B i j code q : ℕ) : Bool :=
  decide (gateMessage s B i j code = q % 2 ^ gateEdge (pA s) (fw s) i j)

theorem gateSide_polyT {α : Type} [Enc α] {n B i j code q : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hi : PolyT i) (hj : PolyT j) (hc : PolyT code) (hq : PolyT q) :
    PolyT fun a => gateSide (sOf (n a)) (B a) (i a) (j a) (code a) (q a) := by
  have hA := (selected_pA_unary hn).poly
  have hF := selected_fw_unary hn
  have he := gateEdge_polyT hA hF.poly hi hj
  have he' : Unary fun a => gateEdge (pA (sOf (n a))) (fw (sOf (n a))) (i a) (j a) :=
    Unary.of_le hF he (fun a => by unfold gateEdge; split <;> simp)
  exact PolyT.decide_eq (gateMessage_polyT hn hB hi hj hc) (PolyT.mod hq (PolyT.two_pow he'))

theorem gateSide_eq_left (s B N k y : ℕ) (i j : Fin (mulNodes s))
    (z : (mulNet s ℤ N k y).S i) (q : ℕ) (hD : pD s ≤ 2 ^ B) (hZ : pZ s ≤ 2 ^ B) :
    (if gateSide s B i.val j.val ((mulCodingZ s N k y).e i z).val q then (1 : ℤ) else 0) =
      ((mulNet s ℤ N k y).F i j).l z (bdec _ q) := by
  let a := (MG.enum (pA s)).symm i
  let b := (MG.enum (pA s)).symm j
  have hi : (MG.enum (pA s) a).val = i.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) i)
  have hj : (MG.enum (pA s) b).val = j.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) j)
  have hc : ((mulCodingZ s N k y).e i z).val = (gcode s a z).val := rfl
  have hm := gateMessage_eq s B a b z hD hZ
  have he := gateEdge_eq s a b
  rw [hi, hj, ← hc] at hm
  rw [hi, hj] at he
  unfold gateSide
  rw [hm, he]
  change (if decide (mV s a b z = q % 2 ^ mC s a b) then (1 : ℤ) else 0) =
    if mV s a b z = bnum (bdec (mC s a b) q) then 1 else 0
  rw [bnum_bdec_remainder]
  simp

theorem gateSide_eq_right (s B N k y : ℕ) (i j : Fin (mulNodes s))
    (z : (mulNet s ℤ N k y).S j) (q : ℕ) (hD : pD s ≤ 2 ^ B) (hZ : pZ s ≤ 2 ^ B) :
    (if gateSide s B j.val i.val ((mulCodingZ s N k y).e j z).val q then (1 : ℤ) else 0) =
      ((mulNet s ℤ N k y).F i j).r z (bdec _ q) := by
  let a := (MG.enum (pA s)).symm i
  let b := (MG.enum (pA s)).symm j
  have hi : (MG.enum (pA s) a).val = i.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) i)
  have hj : (MG.enum (pA s) b).val = j.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) j)
  have hc : ((mulCodingZ s N k y).e j z).val = (gcode s b z).val := rfl
  have hm := gateMessage_eq s B b a z hD hZ
  have he := gateEdge_eq s b a
  rw [hj, hi, ← hc] at hm
  rw [hj, hi, mC_symm] at he
  unfold gateSide
  rw [hm, he]
  change (if decide (mV s b a z = q % 2 ^ mC s a b) then (1 : ℤ) else 0) =
    if mV s b a z = bnum (bdec (mC s a b) q) then 1 else 0
  rw [bnum_bdec_remainder]
  simp

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.gateSide_polyT
#print axioms IntegerFactoring.HalfAppendix.gateSide_eq_left
#print axioms IntegerFactoring.HalfAppendix.gateSide_eq_right
