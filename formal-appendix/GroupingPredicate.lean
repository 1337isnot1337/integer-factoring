/-
SPDX-License-Identifier: Apache-2.0

The executable reference grouping predicate, conditional on a polynomial-time
numeric oracle for the original network's edge widths.
-/
import GroupingCodePoly
import CrossQueryCorrect

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def goodCode (K g L code : ℕ) (edge : ℕ → ℕ → ℕ) : Bool :=
  groupSizesOK K g code && crossingOK K g L (cappedGroupingDigit K g code) edge

theorem goodCode_polyT {α : Type} [Enc α] {K g L code : α → ℕ}
    {edge : α → ℕ → ℕ → ℕ}
    (hK : Unary K) (hg : Unary g) (hL : PolyT L) (hc : PolyT code)
    (he : PolyT fun p : (α × ℕ) × ℕ => edge p.1.1 p.1.2 p.2) :
    PolyT fun a => goodCode (K a) (g a) (L a) (code a) (edge a) := by
  have hs : PolyT fun a => groupSizesOK (K a) (g a) (code a) := groupSizesOK_polyT hK hg hc
  have hcolor : PolyT fun p : α × ℕ => cappedGroupingDigit (K p.1) (g p.1) (code p.1) p.2 :=
    cappedGroupingDigit_polyT hK.poly hg hc
  have hx : PolyT fun a => crossingOK (K a) (g a) (L a) (cappedGroupingDigit (K a) (g a) (code a)) (edge a) :=
    crossingOK_polyT (color := fun a i => cappedGroupingDigit (K a) (g a) (code a) i) hK hg hL hcolor he
  exact PolyT.and' hs hx

theorem goodCode_eq_okG {R : Type} [CommRing R] {g K : ℕ}
    (N : Net R g) (L code : ℕ) (hK : 0 < K) (edge : ℕ → ℕ → ℕ)
    (he : ∀ i j : Fin g, edge i.val j.val = N.c i j) :
    goodCode K g L code edge = true ↔ okG N K L hK code := by
  have hcolor : ∀ i : Fin g, cappedGroupingDigit K g code i.val = (chiOf K g code hK i).val :=
    fun i => cappedGroupingDigit_chiOf K g code hK i
  unfold goodCode okG
  rw [Bool.and_eq_true, groupSizesOK_eq K g code hK,
    crossingOK_eq N (chiOf K g code hK) (cappedGroupingDigit K g code) edge L hcolor he]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.goodCode_polyT
#print axioms IntegerFactoring.HalfAppendix.goodCode_eq_okG
