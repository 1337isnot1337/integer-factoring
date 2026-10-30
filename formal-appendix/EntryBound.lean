/-
SPDX-License-Identifier: Apache-2.0

A general entry-size bound for a reference clique whose network weights and
factors lie between zero and one. Verifying these hypotheses for the grouped
multiplication network is a separate obligation.
-/
import FactoringAux.Clique
import Mathlib

open Finset Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem clique_entry_le_card {K L : ℕ} (N : Net ℤ K)
    (hW : ∀ a z, 0 ≤ N.W a z ∧ N.W a z ≤ 1)
    (hleft : ∀ a b z t, 0 ≤ (N.padF L a b).l z t ∧ (N.padF L a b).l z t ≤ 1)
    (hright : ∀ a b z t, 0 ≤ (N.padF L a b).r z t ∧ (N.padF L a b).r z t ≤ 1)
    (a : Fin K) (θ : Fin K → Bits L) :
    ((N.toClique L) a θ).natAbs ≤ Fintype.card (N.S a) := by
  let f : N.S a → ℤ := fun z => N.W a z *
    (∏ b ∈ univ.filter (a < ·), (N.padF L a b).l z (θ b)) *
    ∏ b ∈ univ.filter (· < a), (N.padF L b a).r z (θ b)
  have hf : ∀ z, 0 ≤ f z ∧ f z ≤ 1 := by
    intro z
    have hL0 : 0 ≤ ∏ b ∈ univ.filter (a < ·), (N.padF L a b).l z (θ b) :=
      Finset.prod_nonneg fun b _ => (hleft a b z (θ b)).1
    have hL1 : (∏ b ∈ univ.filter (a < ·), (N.padF L a b).l z (θ b)) ≤ 1 :=
      Finset.prod_le_one (fun b _ => (hleft a b z (θ b)).1) (fun b _ => (hleft a b z (θ b)).2)
    have hR0 : 0 ≤ ∏ b ∈ univ.filter (· < a), (N.padF L b a).r z (θ b) :=
      Finset.prod_nonneg fun b _ => (hright b a z (θ b)).1
    have hR1 : (∏ b ∈ univ.filter (· < a), (N.padF L b a).r z (θ b)) ≤ 1 :=
      Finset.prod_le_one (fun b _ => (hright b a z (θ b)).1) (fun b _ => (hright b a z (θ b)).2)
    have hw := hW a z
    dsimp [f]
    constructor
    · exact mul_nonneg (mul_nonneg hw.1 hL0) hR0
    · nlinarith [mul_nonneg (sub_nonneg.mpr hw.2) hL0,
        mul_nonneg (sub_nonneg.mpr hL1) hR0]
  have hsum0 : 0 ≤ ∑ z, f z := Finset.sum_nonneg fun z _ => (hf z).1
  have hsum1 : (∑ z, f z) ≤ (Fintype.card (N.S a) : ℤ) := by
    calc
      (∑ z, f z) ≤ ∑ _z : N.S a, (1 : ℤ) := Finset.sum_le_sum fun z _ => (hf z).2
      _ = (Fintype.card (N.S a) : ℤ) := by simp
  change (∑ z, f z).natAbs ≤ Fintype.card (N.S a)
  have habs : ((∑ z, f z).natAbs : ℤ) ≤ (Fintype.card (N.S a) : ℤ) := by
    simpa only [Int.natCast_natAbs, abs_of_nonneg hsum0] using hsum1
  exact_mod_cast habs

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.clique_entry_le_card
