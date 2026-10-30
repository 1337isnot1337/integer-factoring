/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time numeric gate-weight checks for the reference multiplication
network. Message checks and row-summand composition remain open.
-/
import InputCheck
import CarryCheck
import CodingBounds

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def gateWeight (s B N k y i code : ℕ) : Bool :=
  let A := pA s
  let layer := gateLayer A i
  let offset := gateOffset A i
  if layer = 0 then inputCheck A (ph s) (pL s) B k y offset code
  else if layer = 4 then carryCheck A (pL s) (pM s) (fw s) N offset code
  else true

theorem gateWeight_polyT {α : Type} [Enc α] {n B : α → ℕ} {N k y i code : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hN : PolyT N) (hk : PolyT k) (hy : PolyT y)
    (hi : PolyT i) (hc : PolyT code) : PolyT fun a => gateWeight (sOf (n a)) (B a) (N a) (k a) (y a) (i a) (code a) := by
  have hA := selected_pA_unary hn
  have hh : Unary fun a => ph (sOf (n a)) := selected_S_unary hn
  have hL := selected_pL_unary hn
  have hF := selected_fw_unary hn
  have hM := selected_pM_polyT hn
  have hl := gateLayer_polyT hA.poly hi
  have ho := gateOffset_polyT hA.poly hi
  have hin : PolyT fun a => inputCheck (pA (sOf (n a))) (ph (sOf (n a))) (pL (sOf (n a)))
      (B a) (k a) (y a) (gateOffset (pA (sOf (n a))) (i a)) (code a) :=
    inputCheck_polyT (A := fun a => pA (sOf (n a))) (h := fun a => ph (sOf (n a)))
      (L := fun a => pL (sOf (n a))) (B := B) (k := k) (y := y)
      (r := fun a => gateOffset (pA (sOf (n a))) (i a)) (code := code) hA hh hL hB hk hy ho hc
  have hcarry : PolyT fun a => carryCheck (pA (sOf (n a))) (pL (sOf (n a))) (pM (sOf (n a)))
      (fw (sOf (n a))) (N a) (gateOffset (pA (sOf (n a))) (i a)) (code a) :=
    carryCheck_polyT (A := fun a => pA (sOf (n a))) (L := fun a => pL (sOf (n a)))
      (F := fun a => fw (sOf (n a))) (M := fun a => pM (sOf (n a))) (N := N)
      (j := fun a => gateOffset (pA (sOf (n a))) (i a)) (code := code) hA hL hF hM hN ho hc
  exact PolyT.ite' (PolyT.decide_eq hl (PolyT.const 0)) hin
    (PolyT.ite' (PolyT.decide_eq hl (PolyT.const 4)) hcarry (PolyT.const true))

set_option backward.isDefEq.respectTransparency false in
theorem gateWeight_eq (s B N k y : ℕ) (gate : MG (pA s)) (z : mS s gate)
    (hD : ∀ r, gate = .faP r → pD s ≤ 2 ^ B) :
    gateWeight s B N k y (MG.enum (pA s) gate).val (gcode s gate z).val =
      decide (mW s N k y gate z = (1 : ℤ)) := by
  have h := gate_decode (pA s) gate
  cases gate with
  | faP r =>
      have hi := inputCheck_eq s B N k y r z (hD r rfl)
      simp only [gateWeight, h.1, h.2, gcode, if_pos]
      exact hi
  | faQ r => simp_all [gateWeight, mW]
  | mid r => simp_all [gateWeight, mW]
  | ib r => simp_all [gateWeight, mW]
  | carry r =>
      have hc := carryCheck_eq s N k y r z
      simpa [gateWeight, h] using hc

theorem selected_gateWeight_eq (s N k y P : ℕ) (hP : 2 ≤ P) (γ : Fin P)
    (cluster : Fin (3 * P ^ 8)) (gate : MG (pA s)) (z : mS s gate)
    (hc : selectedTopGrouping s N k y P hP cluster = γ)
    (hg : clMG s (3 * P ^ 8) (by positivity) gate = cluster) :
    gateWeight s (topEntryBits s P) N k y (MG.enum (pA s) gate).val (gcode s gate z).val =
      decide (mW s N k y gate z = (1 : ℤ)) := by
  apply gateWeight_eq
  intro r hr
  subst gate
  simpa [gsz] using selected_gate_kind_size_bound s N k y P hP γ cluster (.faP r) hc hg

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.gateWeight_polyT
#print axioms IntegerFactoring.HalfAppendix.gateWeight_eq
#print axioms IntegerFactoring.HalfAppendix.selected_gateWeight_eq
