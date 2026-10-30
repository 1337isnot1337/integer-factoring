/-
SPDX-License-Identifier: Apache-2.0

Executable unary bounds for the multiplication-network base parameters.
-/
import SizeParamPoly
import TopArithmetic
import BoundedTop

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

theorem selected_S_unary {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    Unary fun a => 2 ^ sOf (n a) := by
  have hp : PolyT fun a => 2 ^ sOf (n a) := PolyT.two_pow (sOf_unary hn)
  have hb : Unary fun a => 4 * n a + 4 := by fun_prop
  apply Unary.of_le hb hp
  intro a
  have h := S_le_n (n a)
  have hself : 2 ^ sOf (n a) ≤ (2 ^ sOf (n a)) ^ 4 := Nat.le_self_pow (by norm_num) _
  omega

theorem selected_pA_unary {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    Unary fun a => pA (sOf (n a)) := by
  have hs := selected_S_unary hn
  have hp : Unary fun a => 2 * 2 ^ sOf (n a) := Unary.mul (Unary.const 2) hs
  apply hp.congr
  intro a
  change List.replicate (2 * 2 ^ sOf (n a)) () = List.replicate (2 ^ (sOf (n a) + 1)) ()
  congr 1
  rw [pow_succ]
  ring

theorem selected_pL_unary {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    Unary fun a => pL (sOf (n a)) := by
  have hs := selected_S_unary hn
  have hp : Unary fun a => 2 * (2 ^ sOf (n a) * 2 ^ sOf (n a)) :=
    Unary.mul (Unary.const 2) (Unary.mul hs hs)
  apply hp.congr
  intro a
  change List.replicate (2 * (2 ^ sOf (n a) * 2 ^ sOf (n a))) () =
    List.replicate (2 ^ (2 * sOf (n a) + 1)) ()
  congr 1
  rw [pow_succ, Nat.mul_comm 2 (sOf (n a)), pow_mul, pow_two]
  ring

theorem selected_fw_unary {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    Unary fun a => fw (sOf (n a)) := by
  have hp := selected_pL_unary hn
  have hf : Unary fun a => 4 * pL (sOf (n a)) + 1 := Unary.add (Unary.mul (Unary.const 4) hp) (Unary.const 1)
  exact hf.congr fun a => by rfl

@[fun_prop] theorem unary_div {α : Type} [Enc α] {f g : α → ℕ} (hf : Unary f) (hg : PolyT g) :
    Unary fun a => f a / g a := by
  exact Unary.of_le hf (PolyT.div hf.poly hg) (fun a => Nat.div_le_self _ _)

@[fun_prop] theorem unary_pow_const {α : Type} [Enc α] {f : α → ℕ} (hf : Unary f) (k : ℕ) :
    Unary fun a => f a ^ k := by
  induction k with
  | zero => simpa using (Unary.const 1 : Unary fun _ : α => 1)
  | succ k ih =>
      have h := Unary.mul ih hf
      exact h.congr fun a => by simp [pow_succ]

theorem selected_topWidth_unary (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : Unary fun a => topWidth (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))) + 1 := by
  have hK := selected_K_unary x hx hn
  have hA := selected_pA_unary hn
  have hfw := selected_fw_unary hn
  have hg := Unary.mul (Unary.const 3) (unary_pow_const hK 8)
  have ha_div := unary_div hA hg.poly
  have haa := Unary.mul hA hA
  have haa_div := unary_div haa hg.poly
  have hcl := Unary.add
    (Unary.mul (Unary.mul (Unary.mul (Unary.const 7) (Unary.add ha_div (Unary.const 2))) hA) hfw)
    (Unary.mul (Unary.mul (Unary.const 3) (Unary.add haa_div (Unary.const 2))) hfw)
  have htop := unary_div
    (Unary.add (Unary.mul (Unary.mul (Unary.mul hK (Unary.const 10)) haa) hfw) (Unary.mul hg hcl))
    (unary_pow_const hK 3).poly
  apply (Unary.add htop (Unary.const 1)).congr
  intro a
  simp [topWidth, clP, Nat.mul_assoc]

theorem selected_topEntryBits_unary (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : Unary fun a => topEntryBits (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))) := by
  have hK := selected_K_unary x hx hn
  have hA := selected_pA_unary hn
  have hfw := selected_fw_unary hn
  have hg := Unary.mul (Unary.const 3) (unary_pow_const hK 8)
  have ha_div := unary_div hA hg.poly
  have haa := Unary.mul hA hA
  have haa_div := unary_div haa hg.poly
  have hw := Unary.add
    (Unary.mul (Unary.mul (Unary.mul (Unary.const 8) (Unary.add ha_div (Unary.const 2))) hA) hfw)
    (Unary.mul (Unary.mul (Unary.const 2) (Unary.add haa_div (Unary.const 2))) hfw)
  have hgroup := unary_div (Unary.mul (Unary.const 2) hg) hK.poly
  exact Unary.add (Unary.mul hw hgroup) (Unary.const 1)

theorem selected_rate_unary (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : Unary fun a => 2 ^ (selectedLevel x (n a) / 2) := by
  have hl := selectedLevel_unary x hn
  have hhalf := unary_div hl (PolyT.const 2)
  have hp := PolyT.two_pow hhalf
  exact Unary.of_le (selected_K_unary x hx hn) hp
    (fun a => half_level_le_K x (selectedLevel x (n a)) hx)

theorem selected_tableExponent_unary (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : Unary fun a =>
      (Levels.Gen.K x (selectedLevel x (n a)) - 1) *
        (topWidth (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))) + 1) := by
  have hK := selected_K_unary x hx hn
  have hp : PolyT fun a => Levels.Gen.K x (selectedLevel x (n a)) - 1 := PolyT.sub hK.poly (PolyT.const 1)
  have hsub := Unary.of_le hK hp (fun a => Nat.sub_le _ _)
  exact Unary.mul hsub (selected_topWidth_unary x hx hn)

theorem selected_padExponent_unary (x C : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : Unary fun a => IntegerFactoring.HalfPad.padExponent
      (Levels.Gen.K x (selectedLevel x (n a)))
      (topWidth (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))) + 1)
      (topEntryBits (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a)))) C
      (2 ^ (selectedLevel x (n a) / 2)) := by
  have hK := selected_K_unary x hx hn
  have hL := selected_topWidth_unary x hx hn
  have hB := selected_topEntryBits_unary x hx hn
  have hr := selected_rate_unary x hx hn
  have hK2 := unary_pow_const hK 2
  have hQ := Unary.mul hK (Unary.add hL hB)
  have hlogp : PolyT fun a => Nat.log2 ((Levels.Gen.K x (selectedLevel x (n a))) *
      ((topWidth (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))) + 1) +
        topEntryBits (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))))) := PolyT.log2' hQ.poly
  have hlog : Unary fun a => Nat.log2 ((Levels.Gen.K x (selectedLevel x (n a))) *
      ((topWidth (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))) + 1) +
        topEntryBits (sOf (n a)) (Levels.Gen.K x (selectedLevel x (n a))))) := by
    apply Unary.of_le hQ hlogp
    intro a
    rw [Nat.log2_eq_log_two]
    exact Nat.log_le_self _ _
  have h := Unary.add (Unary.add (unary_div (Unary.mul hK2 hL) hr.poly) (Unary.const 1))
    (Unary.mul (Unary.mul (Unary.const C) hK2) (Unary.add hlog (Unary.const 1)))
  apply h.congr
  intro a
  simp [IntegerFactoring.HalfPad.padExponent, Nat.log2_eq_log_two]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_S_unary
#print axioms IntegerFactoring.HalfAppendix.selected_pA_unary
#print axioms IntegerFactoring.HalfAppendix.selected_pL_unary
#print axioms IntegerFactoring.HalfAppendix.selected_fw_unary
#print axioms IntegerFactoring.HalfAppendix.selected_topWidth_unary
#print axioms IntegerFactoring.HalfAppendix.selected_topEntryBits_unary
#print axioms IntegerFactoring.HalfAppendix.selected_tableExponent_unary
#print axioms IntegerFactoring.HalfAppendix.selected_padExponent_unary
