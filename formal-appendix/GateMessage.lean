/-
SPDX-License-Identifier: Apache-2.0

Numeric reference gate-message computation. Full row-summand composition is
still a separate obligation.
-/
import GateTransforms

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def gateMessage (s B i j code : ℕ) : ℕ :=
  let A := pA s
  let ti := gateLayer A i
  let tj := gateLayer A j
  let ri := gateOffset A i
  let rj := gateOffset A j
  if (ti = 0 ∨ ti = 1) ∧ tj = 2 then faMessage s B ri rj code
  else if ti = 2 ∧ tj = 0 then residueVector s B (midLeft (pZ s) code) rj
  else if ti = 2 ∧ tj = 1 then residueVector s B (midRight (pZ s) code) rj
  else if ti = 2 ∧ tj = 3 then midMessage s B ri rj code
  else if ti = 3 ∧ tj = 2 then residueVector s B code rj
  else if ti = 3 ∧ tj = 4 ∧ rj % A = ri then ibMessage s B (rj / A) code
  else if ti = 4 ∧ tj = 3 ∧ ri % A = rj then carryResidue (pM s) (fw s) code
  else if ti = 4 ∧ tj = 4 ∧ rj = ri + 1 then (carryIncoming (fw s) code + carryResidue (pM s) (fw s) code) / 2 ^ pL s
  else if ti = 4 ∧ tj = 4 ∧ ri = rj + 1 then carryIncoming (fw s) code
  else 0

theorem gateMessage_polyT {α : Type} [Enc α] {n B i j code : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hi : PolyT i) (hj : PolyT j) (hc : PolyT code) :
    PolyT fun a => gateMessage (sOf (n a)) (B a) (i a) (j a) (code a) := by
  have hA := selected_pA_unary hn
  have hM := selected_pM_polyT hn
  have hZ := selected_pZ_polyT hn
  have hL := selected_pL_unary hn
  have hF := selected_fw_unary hn
  have hti := gateLayer_polyT hA.poly hi
  have htj := gateLayer_polyT hA.poly hj
  have hri := gateOffset_polyT hA.poly hi
  have hrj := gateOffset_polyT hA.poly hj
  have hfa := faMessage_polyT hn hB hri hrj hc
  have hmid := midMessage_polyT hn hB hri hrj hc
  have hleft : PolyT fun a => midLeft (pZ (sOf (n a))) (code a) := PolyT.mod (PolyT.div hc hZ) hZ
  have hright : PolyT fun a => midRight (pZ (sOf (n a))) (code a) := PolyT.mod hc hZ
  have hvl := residueVector_polyT hn hB hleft
  have hvr := residueVector_polyT hn hB hright
  have hv := residueVector_polyT hn hB hc
  have hmap : PolyT fun a => (a, gateOffset (pA (sOf (n a))) (j a)) := PolyT.prodMk PolyT.id' hrj
  have hleft' : PolyT fun a => residueVector (sOf (n a)) (B a) (midLeft (pZ (sOf (n a))) (code a)) (gateOffset (pA (sOf (n a))) (j a)) :=
    PolyT.comp (f := fun p : α × ℕ => residueVector (sOf (n p.1)) (B p.1) (midLeft (pZ (sOf (n p.1))) (code p.1)) p.2)
      (g := fun a => (a, gateOffset (pA (sOf (n a))) (j a))) hvl hmap
  have hright' : PolyT fun a => residueVector (sOf (n a)) (B a) (midRight (pZ (sOf (n a))) (code a)) (gateOffset (pA (sOf (n a))) (j a)) :=
    PolyT.comp (f := fun p : α × ℕ => residueVector (sOf (n p.1)) (B p.1) (midRight (pZ (sOf (n p.1))) (code p.1)) p.2)
      (g := fun a => (a, gateOffset (pA (sOf (n a))) (j a))) hvr hmap
  have hv' : PolyT fun a => residueVector (sOf (n a)) (B a) (code a) (gateOffset (pA (sOf (n a))) (j a)) :=
    PolyT.comp (f := fun p : α × ℕ => residueVector (sOf (n p.1)) (B p.1) (code p.1) p.2)
      (g := fun a => (a, gateOffset (pA (sOf (n a))) (j a))) hv hmap
  have hdiv := PolyT.div hrj hA.poly
  have hib := ibMessage_polyT hn hB hdiv hc
  have hFpow := PolyT.two_pow hF
  have hLpow := PolyT.two_pow hL
  have hres : PolyT fun a => carryResidue (pM (sOf (n a))) (fw (sOf (n a))) (code a) :=
    PolyT.mod (PolyT.div hc hFpow) hM
  have hinc : PolyT fun a => carryIncoming (fw (sOf (n a))) (code a) := PolyT.mod hc hFpow
  have hout := PolyT.div (PolyT.add hinc hres) hLpow
  have hi0 := PolyT.decide_eq hti (PolyT.const 0)
  have hi1 := PolyT.decide_eq hti (PolyT.const 1)
  have hi2 := PolyT.decide_eq hti (PolyT.const 2)
  have hi3 := PolyT.decide_eq hti (PolyT.const 3)
  have hi4 := PolyT.decide_eq hti (PolyT.const 4)
  have hj0 := PolyT.decide_eq htj (PolyT.const 0)
  have hj1 := PolyT.decide_eq htj (PolyT.const 1)
  have hj2 := PolyT.decide_eq htj (PolyT.const 2)
  have hj3 := PolyT.decide_eq htj (PolyT.const 3)
  have hj4 := PolyT.decide_eq htj (PolyT.const 4)
  have hc1 := PolyT.decide_and (PolyT.decide_or hi0 hi1) hj2
  have hc2 := PolyT.decide_and hi2 hj0
  have hc3 := PolyT.decide_and hi2 hj1
  have hc4 := PolyT.decide_and hi2 hj3
  have hc5 := PolyT.decide_and hi3 hj2
  have hc6 := PolyT.decide_and hi3 (PolyT.decide_and hj4 (PolyT.decide_eq (PolyT.mod hrj hA.poly) hri))
  have hc7 := PolyT.decide_and hi4 (PolyT.decide_and hj3 (PolyT.decide_eq (PolyT.mod hri hA.poly) hrj))
  have hc8 := PolyT.decide_and hi4 (PolyT.decide_and hj4 (PolyT.decide_eq hrj (PolyT.add hri (PolyT.const 1))))
  have hc9 := PolyT.decide_and hi4 (PolyT.decide_and hj4 (PolyT.decide_eq hri (PolyT.add hrj (PolyT.const 1))))
  exact PolyT.ite' hc1 hfa (PolyT.ite' hc2 hleft' (PolyT.ite' hc3 hright'
    (PolyT.ite' hc4 hmid (PolyT.ite' hc5 hv' (PolyT.ite' hc6 hib (PolyT.ite' hc7 hres
      (PolyT.ite' hc8 hout (PolyT.ite' hc9 hinc (PolyT.const 0)))))))))

set_option backward.isDefEq.respectTransparency false in
theorem gateMessage_eq (s B : ℕ) (source target : MG (pA s)) (z : mS s source)
    (hD : pD s ≤ 2 ^ B) (hZ : pZ s ≤ 2 ^ B) :
    gateMessage s B (MG.enum (pA s) source).val (MG.enum (pA s) target).val (gcode s source z).val =
      mV s source target z := by
  have hs := gate_decode (pA s) source
  have ht := gate_decode (pA s) target
  have hv : ∀ (v : Fin (pA s) → ZMod (pM s)) (r : Fin (pA s)),
      residueVector s B (zCode s v).val r.val = (v r).val := by
    intro v r
    simp [residueVector, r.isLt, fixed_residue_digit_of_size s B v r hZ]
  have hfaP : ∀ (r k : Fin (pA s)) (w : mS s (.faP r)),
      faMessage s B r.val k.val (gcode s (.faP r) w).val = (faV s r.val w k.val).val := by
    intro r k w
    exact faMessage_eq s B r k w hD
  have hfaQ : ∀ (r k : Fin (pA s)) (w : mS s (.faQ r)),
      faMessage s B r.val k.val (gcode s (.faQ r) w).val = (faV s r.val w k.val).val := by
    intro r k w
    exact faMessage_eq s B r k w hD
  have hmid : ∀ (k j : Fin (pA s)) (w : mS s (.mid k)),
      midMessage s B k.val j.val (gcode s (.mid k) w).val = (midV s k.val w j.val).val := by
    intro k j w
    exact midMessage_eq s B k j w hZ
  have hleft : ∀ (k r : Fin (pA s)) (w : mS s (.mid k)),
      residueVector s B (midLeft (pZ s) (gcode s (.mid k) w).val) r.val = (w.1 r).val := by
    intro k r w
    change residueVector s B ((gcode s (.mid k) w).val / pZ s % pZ s) r.val = _
    rw [mid_dig1]
    exact hv w.1 r
  have hright : ∀ (k r : Fin (pA s)) (w : mS s (.mid k)),
      residueVector s B (midRight (pZ s) (gcode s (.mid k) w).val) r.val = (w.2 r).val := by
    intro k r w
    change residueVector s B ((gcode s (.mid k) w).val % pZ s) r.val = _
    rw [mid_dig2]
    exact hv w.2 r
  have hibv : ∀ (j r : Fin (pA s)) (w : mS s (.ib j)),
      residueVector s B (gcode s (.ib j) w).val r.val = (w r).val := by
    intro j r w
    exact hv w r
  have hib : ∀ (j : Fin (pA s)) (w : mS s (.ib j)) (q : Fin (pA s * pA s)),
      ibMessage s B (q.val / pA s) (gcode s (.ib j) w).val = (ibV s w (q.val / pA s)).val := by
    intro j w q
    have hq : q.val / pA s < pA s := by
      apply (Nat.div_lt_iff_lt_mul (pA_pos s)).2
      exact q.isLt
    exact ibMessage_eq s B j w ⟨q.val / pA s, hq⟩ hZ
  have hres : ∀ (j : Fin (pA s * pA s)) (w : mS s (.carry j)),
      carryResidue (pM s) (fw s) (gcode s (.carry j) w).val = w.1.val := carry_y s
  have hinc : ∀ (j : Fin (pA s * pA s)) (w : mS s (.carry j)),
      carryIncoming (fw s) (gcode s (.carry j) w).val = w.2.val := carry_c s
  cases source <;> cases target <;>
    simp_all [gateMessage, mV, hfaP, hfaQ, hmid, hleft, hright, hibv, hib, hres, hinc]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.gateMessage_polyT
#print axioms IntegerFactoring.HalfAppendix.gateMessage_eq
