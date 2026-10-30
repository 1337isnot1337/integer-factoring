/-
SPDX-License-Identifier: Apache-2.0

Numeric polynomial-time check for the multiplication-network carry weight.
The input-gate weights and gate messages remain separate obligations.
-/
import GateFields

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def carryCheck (A L M F N j code : ℕ) : Bool :=
  let r := carryResidue M F code
  let c := carryIncoming F code
  let j' := min j (A * A)
  decide ((c + r) % 2 ^ L = dgt L N j' ∧ (j = 0 → c = 0) ∧
    (j + 1 = A * A → (c + r) / 2 ^ L = 0))

theorem carryCheck_polyT {α : Type} [Enc α] {A L F : α → ℕ} {M N j code : α → ℕ}
    (hA : Unary A) (hL : Unary L) (hF : Unary F)
    (hM : PolyT M) (hN : PolyT N) (hj : PolyT j) (hc : PolyT code) :
    PolyT fun a => carryCheck (A a) (L a) (M a) (F a) (N a) (j a) (code a) := by
  have hA2 := Unary.mul hA hA
  have hidxp := PolyT.min' hj hA2.poly
  have hidx : Unary fun a => min (j a) (A a * A a) := Unary.of_le hA2 hidxp (fun a => Nat.min_le_right _ _)
  have hexp := Unary.mul hL hidx
  have hpow := PolyT.two_pow hexp
  have hpowL := PolyT.two_pow hL
  have hpowF := PolyT.two_pow hF
  have hr : PolyT fun a => carryResidue (M a) (F a) (code a) := PolyT.mod (PolyT.div hc hpowF) hM
  have hcarry : PolyT fun a => carryIncoming (F a) (code a) := PolyT.mod hc hpowF
  have hdgt : PolyT fun a => dgt (L a) (N a) (min (j a) (A a * A a)) :=
    PolyT.mod (PolyT.div hN hpow) hpowL
  have hsum := PolyT.add hcarry hr
  have hfirst := PolyT.decide_eq (PolyT.mod hsum hpowL) hdgt
  have hsecond := PolyT.or' (PolyT.not' (PolyT.decide_eq hj (PolyT.const 0)))
    (PolyT.decide_eq hcarry (PolyT.const 0))
  have hthird := PolyT.or' (PolyT.not' (PolyT.decide_eq (PolyT.add hj (PolyT.const 1)) hA2.poly))
    (PolyT.decide_eq (PolyT.div hsum hpowL) (PolyT.const 0))
  have h := PolyT.and' hfirst (PolyT.and' hsecond hthird)
  apply h.congr
  intro a
  simp [carryCheck]

theorem carryCheck_eq (s N k y : ℕ) (j : Fin (pA s * pA s)) (z : mS s (.carry j)) :
    carryCheck (pA s) (pL s) (pM s) (fw s) N j.val (gcode s (.carry j) z).val =
      decide (mW s N k y (.carry j) z = (1 : ℤ)) := by
  have hres : carryResidue (pM s) (fw s) (gcode s (.carry j) z).val = z.1.val := carry_y s j z
  have hinc : carryIncoming (fw s) (gcode s (.carry j) z).val = z.2.val := carry_c s j z
  simp [carryCheck, hres, hinc, Nat.min_eq_left (Nat.le_of_lt j.isLt), mW]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.carryCheck_polyT
#print axioms IntegerFactoring.HalfAppendix.carryCheck_eq
