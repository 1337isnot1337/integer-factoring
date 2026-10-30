/-
SPDX-License-Identifier: Apache-2.0

A uniform polynomial-time crossing-width scan, conditional on polynomial-time
numeric color and edge-width functions. The factoring edge-width oracle and
its equality to the reference network are not supplied here.
-/
import TopParamsPoly

open Groundwork

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def crossWidth (g : ℕ) (color : ℕ → ℕ) (edge : ℕ → ℕ → ℕ) (u v : ℕ) : ℕ :=
  ((List.range g).flatMap fun i => (List.range g).map fun j =>
    if i < j ∧ ((color i = u ∧ color j = v) ∨ (color i = v ∧ color j = u)) then edge i j else 0).sum

theorem crossWidth_polyT {α : Type} [Enc α] {g u v : α → ℕ}
    {color : α → ℕ → ℕ} {edge : α → ℕ → ℕ → ℕ}
    (hg : Unary g) (hu : PolyT u) (hv : PolyT v)
    (hc : PolyT fun p : α × ℕ => color p.1 p.2)
    (he : PolyT fun p : (α × ℕ) × ℕ => edge p.1.1 p.1.2 p.2) :
    PolyT fun a => crossWidth (g a) (color a) (edge a) (u a) (v a) := by
  have hci : PolyT fun p : (α × ℕ) × ℕ => color p.1.1 p.1.2 := PolyT.comp hc PolyT.id'.fst
  have hcj : PolyT fun p : (α × ℕ) × ℕ => color p.1.1 p.2 :=
    PolyT.comp hc (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hu' : PolyT fun p : (α × ℕ) × ℕ => u p.1.1 := PolyT.comp hu PolyT.id'.fst.fst
  have hv' : PolyT fun p : (α × ℕ) × ℕ => v p.1.1 := PolyT.comp hv PolyT.id'.fst.fst
  have ht : PolyT fun p : (α × ℕ) × ℕ =>
      if p.1.2 < p.2 ∧ ((color p.1.1 p.1.2 = u p.1.1 ∧ color p.1.1 p.2 = v p.1.1) ∨
        (color p.1.1 p.1.2 = v p.1.1 ∧ color p.1.1 p.2 = u p.1.1)) then edge p.1.1 p.1.2 p.2 else 0 := by
    have hlt := PolyT.decide_lt (PolyT.id'.fst.snd : PolyT fun p : (α × ℕ) × ℕ => p.1.2) PolyT.id'.snd
    have hfirst := PolyT.decide_and (PolyT.decide_eq hci hu') (PolyT.decide_eq hcj hv')
    have hsecond := PolyT.decide_and (PolyT.decide_eq hci hv') (PolyT.decide_eq hcj hu')
    exact PolyT.ite' (PolyT.decide_and hlt (PolyT.decide_or hfirst hsecond)) he (PolyT.const 0)
  have hr : PolyT fun p : α × ℕ => List.range (g p.1) := PolyT.range (Unary.comp hg PolyT.id'.fst)
  have hm := PolyT.map (f := fun (p : α × ℕ) j =>
    if p.2 < j ∧ ((color p.1 p.2 = u p.1 ∧ color p.1 j = v p.1) ∨
      (color p.1 p.2 = v p.1 ∧ color p.1 j = u p.1)) then edge p.1 p.2 j else 0) ht hr
  have hf := PolyT.flatMap (f := fun a i => (List.range (g a)).map fun j =>
    if i < j ∧ ((color a i = u a ∧ color a j = v a) ∨
      (color a i = v a ∧ color a j = u a)) then edge a i j else 0) hm (PolyT.range hg)
  exact PolyT.sum hf

def crossingOK (K g L : ℕ) (color : ℕ → ℕ) (edge : ℕ → ℕ → ℕ) : Bool :=
  (List.range K).all fun u => (List.range K).all fun v =>
    decide (u < v → crossWidth g color edge u v ≤ L)

theorem crossingOK_polyT {α : Type} [Enc α] {K g L : α → ℕ}
    {color : α → ℕ → ℕ} {edge : α → ℕ → ℕ → ℕ}
    (hK : Unary K) (hg : Unary g) (hL : PolyT L)
    (hc : PolyT fun p : α × ℕ => color p.1 p.2)
    (he : PolyT fun p : (α × ℕ) × ℕ => edge p.1.1 p.1.2 p.2) :
    PolyT fun a => crossingOK (K a) (g a) (L a) (color a) (edge a) := by
  let β := (α × ℕ) × ℕ
  have hg' : Unary fun p : β => g p.1.1 := Unary.comp hg PolyT.id'.fst.fst
  have hc' : PolyT fun p : β × ℕ => color p.1.1.1 p.2 :=
    PolyT.comp hc (PolyT.prodMk PolyT.id'.fst.fst.fst PolyT.id'.snd)
  have he' : PolyT fun p : (β × ℕ) × ℕ => edge p.1.1.1.1 p.1.2 p.2 :=
    PolyT.comp he (PolyT.prodMk (PolyT.prodMk PolyT.id'.fst.fst.fst.fst PolyT.id'.fst.snd) PolyT.id'.snd)
  have hwidth : PolyT fun p : β => crossWidth (g p.1.1) (color p.1.1) (edge p.1.1) p.1.2 p.2 :=
    crossWidth_polyT hg' PolyT.id'.fst.snd PolyT.id'.snd hc' he'
  have hL' : PolyT fun p : β => L p.1.1 := PolyT.comp hL PolyT.id'.fst.fst
  have hlt := PolyT.decide_lt (PolyT.id'.fst.snd : PolyT fun p : β => p.1.2) PolyT.id'.snd
  have hle := PolyT.decide_le hwidth hL'
  have hpred : PolyT fun p : β => decide (p.1.2 < p.2 →
      crossWidth (g p.1.1) (color p.1.1) (edge p.1.1) p.1.2 p.2 ≤ L p.1.1) := by
    have hbool := PolyT.or' (PolyT.not' hlt) hle
    exact hbool.congr fun p => by simp
  have hr' : PolyT fun p : α × ℕ => List.range (K p.1) := PolyT.range (Unary.comp hK PolyT.id'.fst)
  have hall := PolyT.all (f := fun p v => decide (p.2 < v →
    crossWidth (g p.1) (color p.1) (edge p.1) p.2 v ≤ L p.1)) hpred hr'
  exact PolyT.all (f := fun a u => (List.range (K a)).all fun v => decide (u < v →
    crossWidth (g a) (color a) (edge a) u v ≤ L a)) hall (PolyT.range hK)

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.crossWidth_polyT
#print axioms IntegerFactoring.HalfAppendix.crossingOK_polyT
