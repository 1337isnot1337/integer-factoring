/-
SPDX-License-Identifier: Apache-2.0

Concrete padded factoring algorithm in the reference binary encoding.
-/
import SelectedPrefix
import WordEncoding
import SearchPoly
import ReferenceStatement

open Groundwork Groundwork.Prog PaddedDescent.F Factoring.Aux Levels.Gen
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

def selectedFactorNat (x C : ℕ) (run : U × U × PaddedDescent.F.PInst → ℤ)
    (w : List Bool) (pad : List Unit) : ℕ :=
  let N := PolyT.bval w
  if N = 0 then 2 else if N = 1 then 1 else
    search (fun k y => selectedPrefix x C run w.length pad N k y) w.length

def selectedFactor (x C : ℕ) (run : U × U × PaddedDescent.F.PInst → ℤ)
    (p : List Bool × List Unit) : List Bool := Nat.bits (selectedFactorNat x C run p.1 p.2)

def selectedExponent (x C : ℕ) (w : List Bool) : ℕ := masterExponent x C w.length

set_option maxHeartbeats 2000000 in
theorem selectedExponent_polyT (x C : ℕ) (hx : 1 ≤ x) : PolyT (selectedExponent x C) := by
  have hw : Unary fun w : List Bool => w.length := Unary.of_length PolyT.id'
  exact (masterExponent_unary x C hx hw).poly

set_option maxHeartbeats 2000000 in
theorem selectedFactor_polyT (x C : ℕ) (hx : 1 ≤ x) (run : U × U × PaddedDescent.F.PInst → ℤ) (hrun : PolyT run) :
    PolyT (selectedFactor x C run) := by
  let β := List Bool × List Unit
  let α := List Bool × ℕ × ℕ
  have hw : Unary fun a : α => a.1.length := Unary.of_length PolyT.id'.fst
  have hN : PolyT fun a : α => PolyT.bval a.1 := PolyT.bval' PolyT.id'.fst
  have hbase : PolyT fun p : α × List Unit => selectedPrefix x C run p.1.1.length p.2
      (PolyT.bval p.1.1) p.1.2.1 p.1.2.2 :=
    selectedPrefix_polyT x C hx run hrun (n := fun a : α => a.1.length) (N := fun a : α => PolyT.bval a.1)
      (k := fun a : α => a.2.1) (y := fun a : α => a.2.2) hw hN PolyT.id'.snd.fst PolyT.id'.snd.snd
  have hquery : PolyT fun p : β × ℕ × ℕ => selectedPrefix x C run p.1.1.length p.1.2
      (PolyT.bval p.1.1) p.2.1 p.2.2 :=
    PolyT.comp (f := fun p : α × List Unit => selectedPrefix x C run p.1.1.length p.2
      (PolyT.bval p.1.1) p.1.2.1 p.1.2.2)
      (g := fun p : β × ℕ × ℕ => ((p.1.1, p.2), p.1.2)) hbase
      (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd) PolyT.id'.fst.snd)
  have hwidth : Unary fun p : β => p.1.length := Unary.of_length PolyT.id'.fst
  have hsearch : PolyT fun p : β => search (fun k y => selectedPrefix x C run p.1.length p.2 (PolyT.bval p.1) k y) p.1.length :=
    search_polyT (width := fun p : β => p.1.length)
      (c := fun p k y => selectedPrefix x C run p.1.length p.2 (PolyT.bval p.1) k y) hwidth hquery
  have hN' : PolyT fun p : β => PolyT.bval p.1 := PolyT.bval' PolyT.id'.fst
  have hfactor : PolyT fun p : β => if PolyT.bval p.1 = 0 then 2 else if PolyT.bval p.1 = 1 then 1 else
      search (fun k y => selectedPrefix x C run p.1.length p.2 (PolyT.bval p.1) k y) p.1.length :=
    PolyT.ite' (PolyT.decide_eq hN' (PolyT.const 0)) (PolyT.const 2)
      (PolyT.ite' (PolyT.decide_eq hN' (PolyT.const 1)) (PolyT.const 1) hsearch)
  unfold selectedFactor selectedFactorNat
  with_reducible exact PolyT.bits hfactor

set_option maxHeartbeats 2000000 in
theorem selectedFactor_correct :
    ∃ (x C : ℕ) (run : U × U × PaddedDescent.F.PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ N (pad : List Unit), 2 ^ selectedExponent x C (Computability.encodeNat N) ≤ pad.length →
        selectedFactor x C run (Computability.encodeNat N, pad) = Computability.encodeNat (Factoring.factor N) := by
  obtain ⟨x, C, run, hx, hrun, hfactor⟩ := pure_factor_reference_spec
  refine ⟨x, C, run, hx, hrun, ?_⟩
  intro N pad hpad
  have hP : 2 ≤ K x (selectedLevel x N.size) := K_two x (selectedLevel x N.size) (by omega)
  have hpad' : 2 ^ masterExponent x C N.size ≤ pad.length := by
    simpa [selectedExponent, encodeNat_length] using hpad
  have hq : (fun k y => selectedPrefix x C run N.size pad N k y) =
      prefixFromEvaluator x C run (selectedLevel x N.size) N.size N hP := by
    funext k y
    exact selectedPrefix_eq x C N.size N k y run hP pad hpad'
  unfold selectedFactor selectedFactorNat
  dsimp only
  rw [wordCode_encodeNat, encodeNat_length]
  rw [hq]
  have h := hfactor (selectedLevel x N.size) N hP
  unfold pureFactor at h
  rw [← h]
  exact (encodeNat_eq_bits _).symm

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selectedExponent_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedFactor_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedFactor_correct
