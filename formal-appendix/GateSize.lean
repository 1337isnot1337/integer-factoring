/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time numeric radix sizes for the reference multiplication gates.
-/
import NumericGates
import FactoringAux.EngineMul

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

private theorem pp_const (a : ℕ) : ∀ n, pp (fun _ => a) n = a ^ n
  | 0 => rfl
  | n + 1 => by simp [pp, pp_const a n, pow_succ]

theorem selected_pM_polyT {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    PolyT fun a => pM (sOf (n a)) := by
  have hL := selected_pL_unary hn
  have hpm : Unary fun a => 4 * pL (sOf (n a)) := Unary.mul (Unary.const 4) hL
  have hp := PolyT.two_pow hpm
  exact PolyT.add hp (PolyT.const 1)

theorem selected_pZ_polyT {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    PolyT fun a => pZ (sOf (n a)) := by
  have h := PolyT.pow_unary (selected_pM_polyT hn) (selected_pA_unary hn)
  exact h.congr fun a => (pp_const (pM (sOf (n a))) (pA (sOf (n a)))).symm

theorem selected_pD_polyT {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    PolyT fun a => pD (sOf (n a)) := by
  have hbase := PolyT.two_pow (selected_pL_unary hn)
  have hh : Unary fun a => ph (sOf (n a)) := selected_S_unary hn
  have h := PolyT.pow_unary hbase hh
  exact h.congr fun a => (pp_const (2 ^ pL (sOf (n a))) (ph (sOf (n a)))).symm

def numericGateSize (A D Z M F i : ℕ) : ℕ :=
  if i < A + (A + (A + (A + A * A))) then
    let t := gateLayer A i
    if t = 0 ∨ t = 1 then D else if t = 2 then Z * Z else if t = 3 then Z else M * 2 ^ F
  else 1

theorem numericGateSize_polyT {α : Type} [Enc α] {A D Z M F i : α → ℕ}
    (hA : PolyT A) (hD : PolyT D) (hZ : PolyT Z) (hM : PolyT M) (hF : Unary F) (hi : PolyT i) :
    PolyT fun a => numericGateSize (A a) (D a) (Z a) (M a) (F a) (i a) := by
  have ht := gateLayer_polyT hA hi
  have hp := PolyT.two_pow hF
  unfold numericGateSize
  fun_prop

theorem numericGateSize_eq (s : ℕ) (i : Fin (mulNodes s)) :
    numericGateSize (pA s) (pD s) (pZ s) (pM s) (fw s) i.val = mulSz s i.val := by
  obtain ⟨j, rfl⟩ := (MG.enum (pA s)).surjective i
  have h := gate_decode (pA s) j
  have hi := (MG.enum (pA s) j).isLt
  have hsz : mulSz s (MG.enum (pA s) j).val = gsz s j := by
    unfold mulSz
    rw [dif_pos (by simpa [mulNodes] using hi)]
    simp
  rw [hsz]
  cases j <;> simp_all [numericGateSize, gsz, mulNodes]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_pM_polyT
#print axioms IntegerFactoring.HalfAppendix.selected_pZ_polyT
#print axioms IntegerFactoring.HalfAppendix.selected_pD_polyT
#print axioms IntegerFactoring.HalfAppendix.numericGateSize_polyT
#print axioms IntegerFactoring.HalfAppendix.numericGateSize_eq
