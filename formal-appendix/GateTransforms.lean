/-
SPDX-License-Identifier: Apache-2.0

Numeric implementations of the reference gate transform messages. Other
message cases and the full row summand remain separate obligations.
-/
import NumericTransform
import RootFactors
import GateVectors

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def faMessage (s B r k code : ℕ) : ℕ :=
  numericFa (pA s) (pM s) 16 r k (inputVector s B code)

def midMessage (s B k j code : ℕ) : ℕ :=
  numericMid (pA s) (pM s) 16 (inverseRoot (pA s)) k j
    (residueVector s B (midLeft (pZ s) code)) (residueVector s B (midRight (pZ s) code))

def ibMessage (s B j code : ℕ) : ℕ :=
  numericIb (pA s) (pM s) (inverseRoot (pA s)) (normalization s) j (residueVector s B code)

theorem faMessage_polyT {α : Type} [Enc α] {n B r k code : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hr : PolyT r) (hk : PolyT k) (hc : PolyT code) :
    PolyT fun a => faMessage (sOf (n a)) (B a) (r a) (k a) (code a) := by
  have hA := selected_pA_unary hn
  have hM := selected_pM_polyT hn
  have hv := inputVector_polyT hn hB hc
  exact numericFa_polyT (A := fun a => pA (sOf (n a))) (M := fun a => pM (sOf (n a)))
    (root := fun _ => 16) (r := r) (k := k) (f := fun a i => inputVector (sOf (n a)) (B a) (code a) i)
    hA hM (PolyT.const 16) hr hk hv

theorem midMessage_polyT {α : Type} [Enc α] {n B k j code : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hk : PolyT k) (hj : PolyT j) (hc : PolyT code) :
    PolyT fun a => midMessage (sOf (n a)) (B a) (k a) (j a) (code a) := by
  have hA := selected_pA_unary hn
  have hM := selected_pM_polyT hn
  have hZ := selected_pZ_polyT hn
  have hinv := inverseRoot_polyT hA
  have hl : PolyT fun a => midLeft (pZ (sOf (n a))) (code a) := PolyT.mod (PolyT.div hc hZ) hZ
  have hr : PolyT fun a => midRight (pZ (sOf (n a))) (code a) := PolyT.mod hc hZ
  have hvl := residueVector_polyT hn hB hl
  have hvr := residueVector_polyT hn hB hr
  exact numericMid_polyT (A := fun a => pA (sOf (n a))) (M := fun a => pM (sOf (n a)))
    (root := fun _ => 16) (inv := fun a => inverseRoot (pA (sOf (n a)))) (k := k) (j := j)
    (f := fun a i => residueVector (sOf (n a)) (B a) (midLeft (pZ (sOf (n a))) (code a)) i)
    (g := fun a i => residueVector (sOf (n a)) (B a) (midRight (pZ (sOf (n a))) (code a)) i)
    hA hM (PolyT.const 16) hinv hk hj hvl hvr

theorem ibMessage_polyT {α : Type} [Enc α] {n B j code : α → ℕ}
    (hn : Unary n) (hB : Unary B) (hj : PolyT j) (hc : PolyT code) :
    PolyT fun a => ibMessage (sOf (n a)) (B a) (j a) (code a) := by
  have hA := selected_pA_unary hn
  have hM := selected_pM_polyT hn
  have hinv := inverseRoot_polyT hA
  have ht := normalization_polyT hn
  have hv := residueVector_polyT hn hB hc
  exact numericIb_polyT (A := fun a => pA (sOf (n a))) (M := fun a => pM (sOf (n a)))
    (inv := fun a => inverseRoot (pA (sOf (n a)))) (tinv := fun a => normalization (sOf (n a)))
    (j := j) (f := fun a i => residueVector (sOf (n a)) (B a) (code a) i) hA hM hinv ht hj hv

theorem faMessage_eq (s B : ℕ) (r k : Fin (pA s)) (w : Fin (ph s) → Fin (2 ^ pL s))
    (hsize : pD s ≤ 2 ^ B) :
    faMessage s B r.val k.val (dCode s w).val = (faV s r.val w k.val).val := by
  unfold faMessage faV
  rw [numericFa_eq (pA s) (pM s) 16 r.val k.val (Nat.le_of_lt r.isLt) (Nat.le_of_lt k.isLt)
    (inputVector s B (dCode s w).val) (digZ s w) (fun i _ => inputVector_eq s B w hsize i)]
  rfl

theorem midMessage_eq (s B : ℕ) (k j : Fin (pA s)) (z : mS s (.mid k)) (hsize : pZ s ≤ 2 ^ B) :
    midMessage s B k.val j.val (gcode s (.mid k) z).val = (midV s k.val z j.val).val := by
  unfold midMessage midV
  have hl : midLeft (pZ s) (gcode s (.mid k) z).val = (zCode s z.1).val := mid_dig1 s k z
  have hr : midRight (pZ s) (gcode s (.mid k) z).val = (zCode s z.2).val := mid_dig2 s k z
  rw [hl, hr]
  rw [numericMid_eq (pA s) (pM s) 16 (inverseRoot (pA s)) k.val j.val
    (Nat.le_of_lt k.isLt) (Nat.le_of_lt j.isLt)
    (residueVector s B (zCode s z.1).val) (residueVector s B (zCode s z.2).val)
    (ext0 z.1) (ext0 z.2) (fun i _ => residueVector_eq s B z.1 hsize i)
    (fun i _ => residueVector_eq s B z.2 hsize i)]
  rw [inverseRoot_eq]
  rfl

theorem ibMessage_eq (s B : ℕ) (j : Fin (pA s)) (z : mS s (.ib j))
    (k : Fin (pA s)) (hsize : pZ s ≤ 2 ^ B) :
    ibMessage s B k.val (gcode s (.ib j) z).val = (ibV s z k.val).val := by
  change numericIb (pA s) (pM s) (inverseRoot (pA s)) (normalization s) k.val
      (residueVector s B (zCode s z).val) = (ibV s z k.val).val
  unfold ibV
  rw [numericIb_eq (pA s) (pM s) (inverseRoot (pA s)) (normalization s) k.val (Nat.le_of_lt k.isLt)
    (residueVector s B (zCode s z).val) (ext0 z) (fun i _ => residueVector_eq s B z hsize i)]
  rw [inverseRoot_eq, normalization_eq]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.faMessage_polyT
#print axioms IntegerFactoring.HalfAppendix.midMessage_polyT
#print axioms IntegerFactoring.HalfAppendix.ibMessage_polyT
#print axioms IntegerFactoring.HalfAppendix.faMessage_eq
#print axioms IntegerFactoring.HalfAppendix.midMessage_eq
#print axioms IntegerFactoring.HalfAppendix.ibMessage_eq
