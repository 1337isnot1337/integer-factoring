/-
SPDX-License-Identifier: Apache-2.0

Uniform polynomial-time decoding of a grouping code, capped at the unary
number of source groups. This is a component of the grouping predicate.
-/
import TopParamsPoly

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def cappedGroupingDigit (K g code i : ℕ) : ℕ := code / K ^ min i g % K

theorem cappedGroupingDigit_polyT {α : Type} [Enc α] {K g code : α → ℕ}
    (hK : PolyT K) (hg : Unary g) (hc : PolyT code) :
    PolyT fun p : α × ℕ => cappedGroupingDigit (K p.1) (g p.1) (code p.1) p.2 := by
  have hgp : PolyT g := Unary.poly hg
  have hcap : Unary fun p : α × ℕ => min p.2 (g p.1) := by
    have hlim : Unary fun p : α × ℕ => g p.1 := Unary.comp hg (by fun_prop)
    have hp : PolyT fun p : α × ℕ => min p.2 (g p.1) := by fun_prop
    exact Unary.of_le hlim hp (fun p => Nat.min_le_right _ _)
  have hbase : PolyT fun p : α × ℕ => K p.1 := PolyT.comp hK (by fun_prop)
  have hpow := PolyT.pow_unary hbase hcap
  have hcode : PolyT fun p : α × ℕ => code p.1 := PolyT.comp hc (by fun_prop)
  exact PolyT.mod (PolyT.div hcode hpow) hbase

theorem cappedGroupingDigit_eq (K g code i : ℕ) (hi : i ≤ g) :
    cappedGroupingDigit K g code i = code / K ^ i % K := by
  simp [cappedGroupingDigit, Nat.min_eq_left hi]

theorem cappedGroupingDigit_chiOf (K g code : ℕ) (hK : 0 < K) (i : Fin g) :
    cappedGroupingDigit K g code i = (chiOf K g code hK i).val := by
  simp [cappedGroupingDigit, chiOf, Nat.min_eq_left (Nat.le_of_lt i.isLt)]

def groupMembers (K g code group : ℕ) : ℕ :=
  ((List.range g).filter fun i => decide (cappedGroupingDigit K g code i = group)).length

def groupSizesOK (K g code : ℕ) : Bool :=
  (List.range K).all fun group => decide (K * groupMembers K g code group ≤ 2 * g)

theorem groupMembers_polyT {α : Type} [Enc α] {K g code : α → ℕ}
    (hK : PolyT K) (hg : Unary g) (hc : PolyT code) :
    PolyT fun p : α × ℕ => groupMembers (K p.1) (g p.1) (code p.1) p.2 := by
  have hd := cappedGroupingDigit_polyT hK hg hc
  have hd' : PolyT fun p : (α × ℕ) × ℕ => cappedGroupingDigit (K p.1.1) (g p.1.1) (code p.1.1) p.2 :=
    PolyT.comp hd (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hr : PolyT fun p : α × ℕ => List.range (g p.1) := PolyT.range (Unary.comp hg PolyT.id'.fst)
  have hpred := PolyT.decide_eq hd' (PolyT.id'.fst.snd : PolyT fun p : (α × ℕ) × ℕ => p.1.2)
  exact PolyT.length (PolyT.filter hpred hr)

theorem groupSizesOK_polyT {α : Type} [Enc α] {K g code : α → ℕ}
    (hK : Unary K) (hg : Unary g) (hc : PolyT code) :
    PolyT fun a => groupSizesOK (K a) (g a) (code a) := by
  have hKp := Unary.poly hK
  have hgp := Unary.poly hg
  have hm := groupMembers_polyT hKp hg hc
  have hr : PolyT fun a => List.range (K a) := PolyT.range hK
  have hK' : PolyT fun p : α × ℕ => K p.1 := PolyT.comp hKp PolyT.id'.fst
  have hg' : PolyT fun p : α × ℕ => g p.1 := PolyT.comp hgp PolyT.id'.fst
  have hpred := PolyT.decide_le (PolyT.mul hK' hm) (PolyT.mul (PolyT.const 2) hg')
  exact PolyT.all hpred hr

theorem groupMembers_eq_card (K g code : ℕ) (hK : 0 < K) (group : Fin K) :
    groupMembers K g code group.val =
      (Finset.univ.filter fun i : Fin g => chiOf K g code hK i = group).card := by
  let xs := (List.range g).filter fun i => decide (cappedGroupingDigit K g code i = group.val)
  let fs := Finset.univ.filter fun i : Fin g => chiOf K g code hK i = group
  have hnodup : xs.Nodup := (List.nodup_range).filter _
  have heq : xs.toFinset = fs.image Fin.val := by
    ext i
    simp only [xs, List.mem_toFinset, List.mem_filter, List.mem_range, decide_eq_true_eq,
      fs, Finset.mem_image, Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨hi, hd⟩
      refine ⟨⟨i, hi⟩, ?_, rfl⟩
      apply Fin.ext
      exact (cappedGroupingDigit_chiOf K g code hK ⟨i, hi⟩).symm.trans hd
    · rintro ⟨j, hj, rfl⟩
      refine ⟨j.isLt, ?_⟩
      rw [cappedGroupingDigit_chiOf K g code hK j, hj]
  have hc := Finset.card_image_of_injective fs Fin.val_injective
  rw [← hc, ← heq, List.toFinset_card_of_nodup hnodup]
  rfl

theorem groupSizesOK_eq (K g code : ℕ) (hK : 0 < K) :
    groupSizesOK K g code = true ↔
      ∀ group : Fin K, K * (Finset.univ.filter fun i : Fin g => chiOf K g code hK i = group).card ≤ 2 * g := by
  simp only [groupSizesOK, List.all_eq_true, List.mem_range, decide_eq_true_eq]
  constructor
  · intro h group
    simpa [groupMembers_eq_card K g code hK group] using h group.val group.isLt
  · intro h i hi
    simpa [groupMembers_eq_card K g code hK ⟨i, hi⟩] using h ⟨i, hi⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.cappedGroupingDigit_polyT
#print axioms IntegerFactoring.HalfAppendix.cappedGroupingDigit_chiOf
#print axioms IntegerFactoring.HalfAppendix.groupMembers_polyT
#print axioms IntegerFactoring.HalfAppendix.groupSizesOK_polyT
#print axioms IntegerFactoring.HalfAppendix.groupMembers_eq_card
#print axioms IntegerFactoring.HalfAppendix.groupSizesOK_eq
