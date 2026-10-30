/-
SPDX-License-Identifier: Apache-2.0

Numeric decoding of multiplication-network gates and their cluster numbers.
-/
import TopParamsPoly
import FactoringAux.EngineTop

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def gateLayer (A i : ℕ) : ℕ :=
  if i < A then 0 else if i < 2 * A then 1 else if i < 3 * A then 2 else if i < 4 * A then 3 else 4

def gateOffset (A i : ℕ) : ℕ :=
  if i < A then i else if i < 2 * A then i - A else if i < 3 * A then i - 2 * A
  else if i < 4 * A then i - 3 * A else i - 4 * A

def gateCluster (A g i : ℕ) : ℕ :=
  if gateLayer A i < 4 then gateOffset A i * g / A else gateOffset A i * g / (A * A)

theorem gateLayer_polyT {α : Type} [Enc α] {A i : α → ℕ} (hA : PolyT A) (hi : PolyT i) :
    PolyT fun a => gateLayer (A a) (i a) := by
  unfold gateLayer
  fun_prop

theorem gateOffset_polyT {α : Type} [Enc α] {A i : α → ℕ} (hA : PolyT A) (hi : PolyT i) :
    PolyT fun a => gateOffset (A a) (i a) := by
  unfold gateOffset
  fun_prop

theorem gateCluster_polyT {α : Type} [Enc α] {A g i : α → ℕ}
    (hA : PolyT A) (hg : PolyT g) (hi : PolyT i) : PolyT fun a => gateCluster (A a) (g a) (i a) := by
  have hl := gateLayer_polyT hA hi
  have ho := gateOffset_polyT hA hi
  unfold gateCluster
  fun_prop

theorem gate_decode (A : ℕ) (i : MG A) :
    gateLayer A (MG.enum A i).val = (match i with | .faP _ => 0 | .faQ _ => 1 | .mid _ => 2 | .ib _ => 3 | .carry _ => 4) ∧
    gateOffset A (MG.enum A i).val = (match i with | .faP r => r.val | .faQ r => r.val | .mid r => r.val | .ib r => r.val | .carry r => r.val) := by
  cases i with
  | faP r => simp [gateLayer, gateOffset, enum_faP, r.isLt]
  | faQ r => simp only [enum_faQ, gateLayer, gateOffset]; have h := r.isLt; split_ifs <;> omega
  | mid r => simp only [enum_mid, gateLayer, gateOffset]; have h := r.isLt; split_ifs <;> omega
  | ib r => simp only [enum_ib, gateLayer, gateOffset]; have h := r.isLt; split_ifs <;> omega
  | carry r => simp only [enum_carry, gateLayer, gateOffset]; split_ifs <;> omega

theorem gateCluster_eq (s g : ℕ) (hg : 0 < g) (i : MG (pA s)) :
    gateCluster (pA s) g (MG.enum (pA s) i).val = (clMG s g hg i).val := by
  have h := gate_decode (pA s) i
  cases i <;> simp_all [gateCluster, clMG]

def gateLink (A i j : ℕ) : Prop :=
  let ti := gateLayer A i
  let tj := gateLayer A j
  let ri := gateOffset A i
  let rj := gateOffset A j
  ((ti = 0 ∨ ti = 1) ∧ tj = 2) ∨ ((tj = 0 ∨ tj = 1) ∧ ti = 2) ∨
  (ti = 2 ∧ tj = 3) ∨ (ti = 3 ∧ tj = 2) ∨
  (ti = 3 ∧ tj = 4 ∧ rj % A = ri) ∨ (ti = 4 ∧ tj = 3 ∧ ri % A = rj) ∨
  (ti = 4 ∧ tj = 4 ∧ (rj = ri + 1 ∨ ri = rj + 1))

instance (A i j : ℕ) : Decidable (gateLink A i j) := by unfold gateLink; infer_instance

def gateEdge (A width i j : ℕ) : ℕ := if gateLink A i j then width else 0

theorem gateEdge_polyT {α : Type} [Enc α] {A width i j : α → ℕ}
    (hA : PolyT A) (hw : PolyT width) (hi : PolyT i) (hj : PolyT j) :
    PolyT fun a => gateEdge (A a) (width a) (i a) (j a) := by
  have hli := gateLayer_polyT hA hi
  have hlj := gateLayer_polyT hA hj
  have hoi := gateOffset_polyT hA hi
  have hoj := gateOffset_polyT hA hj
  unfold gateEdge gateLink
  fun_prop

theorem gateEdge_eq (s : ℕ) (i j : MG (pA s)) :
    gateEdge (pA s) (fw s) (MG.enum (pA s) i).val (MG.enum (pA s) j).val = mC s i j := by
  have hi := gate_decode (pA s) i
  have hj := gate_decode (pA s) j
  cases i <;> cases j <;> simp_all [gateEdge, gateLink, mC, mLink]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.gateCluster_polyT
#print axioms IntegerFactoring.HalfAppendix.gateCluster_eq
#print axioms IntegerFactoring.HalfAppendix.gateEdge_polyT
#print axioms IntegerFactoring.HalfAppendix.gateEdge_eq
