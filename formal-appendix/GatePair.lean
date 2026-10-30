/-
SPDX-License-Identifier: Apache-2.0

Numeric gate-pair agreement for the reference multiplication network.
-/
import GateMessage
import ExplicitCoding
import FactoringAux.MulCount

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def gatePair (s B i j ci cj : ℕ) : Bool :=
  decide (gateMessage s B i j ci = gateMessage s B j i cj)

theorem gatePair_polyT {α : Type} [Enc α] {n B i j ci cj : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hi : PolyT i) (hj : PolyT j) (hci : PolyT ci) (hcj : PolyT cj) :
    PolyT fun a => gatePair (sOf (n a)) (B a) (i a) (j a) (ci a) (cj a) := by
  have hleft := gateMessage_polyT (n := n) (B := B) (i := i) (j := j) (code := ci) hn hB hi hj hci
  have hright := gateMessage_polyT (n := n) (B := B) (i := j) (j := i) (code := cj) hn hB hj hi hcj
  exact PolyT.decide_eq hleft hright

theorem gatePair_eq (s B N k y : ℕ) (i j : Fin (mulNodes s))
    (z : (mulNet s ℤ N k y).S i) (w : (mulNet s ℤ N k y).S j)
    (hD : pD s ≤ 2 ^ B) (hZ : pZ s ≤ 2 ^ B) :
    (if gatePair s B i.val j.val ((mulCodingZ s N k y).e i z).val ((mulCodingZ s N k y).e j w).val then (1 : ℤ) else 0) =
      (mulNet s ℤ N k y).Ψ i j z w := by
  let a := (MG.enum (pA s)).symm i
  let b := (MG.enum (pA s)).symm j
  have hidxi : (MG.enum (pA s) a).val = i.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) i)
  have hidxj : (MG.enum (pA s) b).val = j.val := congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) j)
  have hci : ((mulCodingZ s N k y).e i z).val = (gcode s a z).val := rfl
  have hcj : ((mulCodingZ s N k y).e j w).val = (gcode s b w).val := rfl
  have hi : gateMessage s B i.val j.val ((mulCodingZ s N k y).e i z).val = mV s a b z := by
    have h := gateMessage_eq s B a b z hD hZ
    rw [hidxi, hidxj, ← hci] at h
    exact h
  have hj : gateMessage s B j.val i.val ((mulCodingZ s N k y).e j w).val = mV s b a w := by
    have h := gateMessage_eq s B b a w hD hZ
    rw [hidxj, hidxi, ← hcj] at h
    exact h
  have hbound : mV s b a w < 2 ^ mC s a b := by
    have h := mV_lt s b a w
    simpa [mC_symm] using h
  unfold gatePair
  simp only [decide_eq_true_eq]
  rw [hi, hj]
  change (if mV s a b z = mV s b a w then (1 : ℤ) else 0) =
    (Fac.numTest (mC s a b) (mV s a b) (mV s b a) : Fac ℤ _ _ _).eval z w
  have heval := Fac.eval_numTest (R := ℤ) (mC s a b) (mV s a b) (mV s b a)
    (show mS s a from z) (show mS s b from w)
  simp only [hbound, and_true] at heval
  exact heval.symm

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.gatePair_polyT
#print axioms IntegerFactoring.HalfAppendix.gatePair_eq
