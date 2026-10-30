/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time finite row-state summation using a supplied unary pad. The
factoring-specific state decoder and summand remain to be constructed.
-/
import GroupingSearchPoly
import PaddedDescent.F.FuelCert

open Groundwork

namespace IntegerFactoring.HalfAppendix

def paddedRowSum (S : ℕ) (pad : List Unit) (f : ℕ → ℕ) : ℕ :=
  ((List.range (min S pad.length)).map f).sum

theorem paddedRowSum_polyT {α : Type} [Enc α] {S : α → ℕ} {f : (α × List Unit) → ℕ → ℕ}
    (hS : PolyT S) (hf : PolyT fun p : (α × List Unit) × ℕ => f p.1 p.2) :
    PolyT fun p : α × List Unit => paddedRowSum (S p.1) p.2 (f p) := by
  have hS' : PolyT fun p : α × List Unit => S p.1 := PolyT.comp hS PolyT.id'.fst
  have hpad : PolyT fun p : α × List Unit => p.2.take (S p.1) := PolyT.take hS' PolyT.id'.snd
  have hbound : Unary fun p : α × List Unit => min (S p.1) p.2.length := by
    exact (Unary.of_length hpad).congr (fun p => by simp [List.length_take])
  have hr : PolyT fun p : α × List Unit => List.range (min (S p.1) p.2.length) := PolyT.range hbound
  have hm : PolyT fun p : α × List Unit => (List.range (min (S p.1) p.2.length)).map (f p) :=
    PolyT.map (f := f) hf hr
  exact PolyT.sum hm

theorem paddedRowSum_eq_finSum (S : ℕ) (pad : List Unit) (f : ℕ → ℕ) (hpad : S ≤ pad.length) :
    paddedRowSum S pad f = ∑ i : Fin S, f i.val := by
  unfold paddedRowSum
  rw [Nat.min_eq_left hpad]
  rw [Fin.sum_univ_eq_sum_range]
  rw [← List.toFinset_range, List.sum_toFinset f List.nodup_range]

def paddedRowSumZ (S : ℕ) (pad : List Unit) (f : ℕ → ℤ) : ℤ :=
  ((List.range (min S pad.length)).map f).sum

theorem paddedRowSumZ_polyT {α : Type} [Enc α] {S : α → ℕ} {f : (α × List Unit) → ℕ → ℤ}
    (hS : PolyT S) (hf : PolyT fun p : (α × List Unit) × ℕ => f p.1 p.2) :
    PolyT fun p : α × List Unit => paddedRowSumZ (S p.1) p.2 (f p) := by
  have hS' : PolyT fun p : α × List Unit => S p.1 := PolyT.comp hS PolyT.id'.fst
  have hpad : PolyT fun p : α × List Unit => p.2.take (S p.1) := PolyT.take hS' PolyT.id'.snd
  have hbound : Unary fun p : α × List Unit => min (S p.1) p.2.length := by
    exact (Unary.of_length hpad).congr (fun p => by simp [List.length_take])
  have hr : PolyT fun p : α × List Unit => List.range (min (S p.1) p.2.length) := PolyT.range hbound
  have hm : PolyT fun p : α × List Unit => (List.range (min (S p.1) p.2.length)).map (f p) :=
    PolyT.map (f := f) hf hr
  exact PaddedDescent.F.sumZ_polyT hm

theorem paddedRowSumZ_eq_finSum (S : ℕ) (pad : List Unit) (f : ℕ → ℤ) (hpad : S ≤ pad.length) :
    paddedRowSumZ S pad f = ∑ i : Fin S, f i.val := by
  unfold paddedRowSumZ
  rw [Nat.min_eq_left hpad]
  rw [Fin.sum_univ_eq_sum_range]
  rw [← List.toFinset_range, List.sum_toFinset f List.nodup_range]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.paddedRowSum_polyT
#print axioms IntegerFactoring.HalfAppendix.paddedRowSum_eq_finSum
#print axioms IntegerFactoring.HalfAppendix.paddedRowSumZ_polyT
#print axioms IntegerFactoring.HalfAppendix.paddedRowSumZ_eq_finSum
