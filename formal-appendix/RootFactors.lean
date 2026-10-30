/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time numeric representatives of the reference transform roots.
-/
import GateFields

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def inverseRoot (A : ℕ) : ℕ := 16 ^ (A * A - 1)

def normalization (s : ℕ) : ℕ := 2 ^ ((2 * pm s - 1) * (2 * s + 2))

theorem inverseRoot_polyT {α : Type} [Enc α] {A : α → ℕ} (hA : Unary A) :
    PolyT fun a => inverseRoot (A a) := by
  have hA2 := Unary.mul hA hA
  have hexpp : PolyT fun a => A a * A a - 1 := PolyT.sub hA2.poly (PolyT.const 1)
  have hexp : Unary fun a => A a * A a - 1 := Unary.of_le hA2 hexpp (fun _ => Nat.sub_le _ _)
  exact PolyT.pow_unary (PolyT.const 16) hexp

theorem inverseRoot_eq (s : ℕ) : (inverseRoot (pA s) : ZMod (pM s)) = ωi s := by
  unfold inverseRoot ωi
  have hA : pA s * pA s = 2 ^ (2 * s + 1 + 1) := by
    unfold pA
    rw [← pow_add]
    congr 1
    omega
  rw [hA]
  push_cast
  rfl

theorem normalization_polyT {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    PolyT fun a => normalization (sOf (n a)) := by
  have hs := sOf_unary hn
  have hL := selected_pL_unary hn
  have hpm : Unary fun a => pm (sOf (n a)) := Unary.mul (Unary.const 4) hL
  have h2pm := Unary.mul (Unary.const 2) hpm
  have hsubp : PolyT fun a => 2 * pm (sOf (n a)) - 1 := PolyT.sub h2pm.poly (PolyT.const 1)
  have hsub : Unary fun a => 2 * pm (sOf (n a)) - 1 := Unary.of_le h2pm hsubp (fun _ => Nat.sub_le _ _)
  have htail : Unary fun a => 2 * sOf (n a) + 2 := Unary.add (Unary.mul (Unary.const 2) hs) (Unary.const 2)
  exact PolyT.two_pow (Unary.mul hsub htail)

theorem normalization_eq (s : ℕ) : (normalization s : ZMod (pM s)) = tInv s := by
  simp [normalization, tInv]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.inverseRoot_polyT
#print axioms IntegerFactoring.HalfAppendix.inverseRoot_eq
#print axioms IntegerFactoring.HalfAppendix.normalization_polyT
#print axioms IntegerFactoring.HalfAppendix.normalization_eq
