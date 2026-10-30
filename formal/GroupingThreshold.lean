/-
SPDX-License-Identifier: Apache-2.0

Adapted from FactoringAux/Grouping.lean and FactoringAux/Lemma1.lean in
https://github.com/IsaacRudich/SubexponentialDeterministicFactoring
at commit c91cc3f4947caca8977e81a5d14244f73b5ae2b1.
Source authors: Isaac Rudich and Louis-Martin Rousseau; source proof and code
generated with Claude, as described in the upstream README.

The original moments argument proves more than its published hypothesis
`g = 3 * K^8`: the final union bound needs only `K^2 + 4*K^6 < g`.
More generally, an independent accuracy parameter D needs only
`K^2 + 4*D^2*K^4 < g`. This file reuses the original moment lemmas, rescales
the deviation by D, and adapts the union-bound threshold. No stronger
probabilistic tail bound is assumed here.

This proves a finite grouping lemma. It does not by itself change or prove
a running-time bound for a factoring machine.
-/
module

public import FactoringAux.Lemma1

@[expose] public section

namespace Factoring.Improvements

open Finset Factoring.Aux

variable {g K : ℕ}

/-- A grouping with discrepancy controlled by the independent parameter D.
All inequalities are integral and contain no division. -/
def GoodAccuracy (c : Fin g → Fin g → ℕ) (P D : ℕ)
    (χ : Fin g → Fin K) : Prop :=
  (∀ γ, (K : ℤ) * grpSize χ γ ≤ 2 * g) ∧
    ∀ γ γ', γ ≠ γ' → (D : ℤ) * (K : ℤ) ^ 2 * grpY c χ γ γ' ≤
      D * totU c + g * P

/-- Independent accuracy D costs `O(D^2 K^4)` vertices by the original
second-moment estimates. No concentration inequality is assumed. -/
theorem exists_good_accuracy_of_threshold (hK : 2 ≤ K) (c : Fin g → Fin g → ℕ) (hdiag : ∀ a, c a a = 0)
    (hsym : ∀ a b, c a b = c b a) (P D : ℕ) (hrow : ∀ a, ∑ b, c a b ≤ P) (hg : K ^ 2 + 4 * D ^ 2 * K ^ 4 < g) :
    ∃ χ : Fin g → Fin K, GoodAccuracy c P D χ := by
  have hKsq : 4 ≤ K ^ 2 := by
    calc 4 = 2 ^ 2 := by norm_num
      _ ≤ K ^ 2 := Nat.pow_le_pow_left hK 2
  have hg4 : 4 ≤ g := by omega
  have hgpos : (0 : ℤ) < g := by exact_mod_cast (by omega : 0 < g)
  have hKpos : (0 : ℤ) < K := by exact_mod_cast (by omega : 0 < K)
  set B1 : Fin K → Finset (Fin g → Fin K) :=
    fun γ => univ.filter (fun χ => 2 * (g : ℤ) < K * grpSize χ γ) with hB1def
  set B2 : Fin K × Fin K → Finset (Fin g → Fin K) :=
    fun q => univ.filter (fun χ => (D : ℤ) * totU c + g * P < (D : ℤ) * (K : ℤ) ^ 2 * grpY c χ q.1 q.2) with hB2def
  -- the group-size bound
  have hB1 : ∀ γ, ((B1 γ).card : ℤ) * g ≤ (K : ℤ) ^ (g + 1) := by
    intro γ
    have h1 := card_mul_le_sum_sq (fun χ => (K : ℤ) * grpSize χ γ - g) g (B1 γ) (fun χ hχ => by
      have := (Finset.mem_filter.mp hχ).2
      exact pow_le_pow_left₀ hgpos.le (by linarith) 2)
    have h2 := sum_size_dev_sq (K := K) (by omega : 2 ≤ g) γ
    have h3 : ((B1 γ).card : ℤ) * g * g ≤ (K : ℤ) ^ (g + 1) * g := by nlinarith
    exact le_of_mul_le_mul_right h3 hgpos
  -- the shared-bits bound
  have hB2 : ∀ q : Fin K × Fin K, ((B2 q).card : ℤ) * g ≤ 4 * (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) := by
    intro q
    have hU := totU_le c P hrow
    have hU0 := totU_nonneg c
    have hP0 : (0 : ℤ) ≤ P := by positivity
    have h1 := card_mul_le_sum_sq (fun χ => (D : ℤ) * ((K : ℤ) ^ 2 * grpY c χ q.1 q.2 - totU c))
      ((g : ℤ) * P + 1) (B2 q) (fun χ hχ => by
        have := (Finset.mem_filter.mp hχ).2
        have h' : (g : ℤ) * P + 1 ≤ (D : ℤ) * ((K : ℤ) ^ 2 * grpY c χ q.1 q.2 - totU c) := by
          have : (D : ℤ) * ((K : ℤ) ^ 2 * grpY c χ q.1 q.2 - totU c) =
              (D : ℤ) * (K : ℤ) ^ 2 * grpY c χ q.1 q.2 - D * totU c := by ring
          rw [this]; linarith
        exact pow_le_pow_left₀ (by positivity) h' 2)
    have h2 := sum_dev_sq_le c hdiag hsym P hrow hg4 q.1 q.2
    have h3 : ∑ χ : Fin g → Fin K, ((D : ℤ) * ((K : ℤ) ^ 2 * grpY c χ q.1 q.2 - totU c)) ^ 2 =
        (D : ℤ) ^ 2 * ∑ χ : Fin g → Fin K, ((K : ℤ) ^ 2 * grpY c χ q.1 q.2 - totU c) ^ 2 := by
      rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun χ _ => by ring
    have hK4 : (0 : ℤ) ≤ (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) := by positivity
    have h4 : ((B2 q).card : ℤ) * ((g : ℤ) * P + 1) ^ 2 ≤ (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) * (4 * P * (g * P)) := by
      rw [h3] at h1
      calc _ ≤ _ := h1
        _ ≤ (D : ℤ) ^ 2 * ((K : ℤ) ^ (g + 2) * (4 * P * totU c)) :=
            mul_le_mul_of_nonneg_left h2 (by positivity)
        _ = (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) * (4 * P * totU c) := by ring
        _ ≤ _ := mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hU (by positivity)) hK4
    have hc0 : (0 : ℤ) ≤ (B2 q).card := by positivity
    have h5 : ((B2 q).card : ℤ) * g * ((g : ℤ) * P + 1) ^ 2 ≤ 4 * (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) * ((g : ℤ) * P + 1) ^ 2 := by
      have hgP : ((g : ℤ) * P) ^ 2 ≤ ((g : ℤ) * P + 1) ^ 2 := by nlinarith
      calc ((B2 q).card : ℤ) * g * ((g : ℤ) * P + 1) ^ 2 = g * (((B2 q).card : ℤ) * ((g : ℤ) * P + 1) ^ 2) := by
            ring
        _ ≤ g * ((D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) * (4 * P * (g * P))) := mul_le_mul_of_nonneg_left h4 hgpos.le
        _ = 4 * (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) * ((g : ℤ) * P) ^ 2 := by ring
        _ ≤ _ := mul_le_mul_of_nonneg_left hgP (by positivity)
    exact le_of_mul_le_mul_right h5 (by positivity)
  -- every grouping is bad, if none is good
  by_contra hno
  push Not at hno
  have hcover : (univ : Finset (Fin g → Fin K)) ⊆ (univ.biUnion B1) ∪ (univ.biUnion B2) := by
    intro χ _
    have hχ := hno χ
    unfold GoodAccuracy at hχ
    simp only [not_and_or, not_forall, not_le, exists_prop] at hχ
    rcases hχ with ⟨γ, hγ⟩ | ⟨γ, γ', _, hγ⟩
    · exact Finset.mem_union_left _ (Finset.mem_biUnion.mpr ⟨γ, mem_univ _, by simp [hB1def, hγ]⟩)
    · exact Finset.mem_union_right _ (Finset.mem_biUnion.mpr ⟨(γ, γ'), mem_univ _, by simp [hB2def, hγ]⟩)
  have hcard : (K : ℤ) ^ g ≤ ∑ γ, ((B1 γ).card : ℤ) + ∑ q, ((B2 q).card : ℤ) := by
    have := Finset.card_le_card hcover
    have h2 := (Finset.card_union_le _ _).trans (add_le_add (Finset.card_biUnion_le (s := univ) (t := B1))
      (Finset.card_biUnion_le (s := univ) (t := B2)))
    rw [Finset.card_univ, Fintype.card_fun, Fintype.card_fin, Fintype.card_fin] at this
    have := this.trans h2
    exact_mod_cast this
  have hsum : (g : ℤ) * (K : ℤ) ^ g ≤ (K : ℤ) * (K : ℤ) ^ (g + 1) + (K : ℤ) ^ 2 * (4 * (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2)) := by
    calc (g : ℤ) * (K : ℤ) ^ g ≤ (g : ℤ) * (∑ γ, ((B1 γ).card : ℤ) + ∑ q, ((B2 q).card : ℤ)) :=
          mul_le_mul_of_nonneg_left hcard hgpos.le
      _ = ∑ γ, ((B1 γ).card : ℤ) * g + ∑ q, ((B2 q).card : ℤ) * g := by
          rw [mul_add, Finset.mul_sum, Finset.mul_sum]
          congr 1 <;> exact Finset.sum_congr rfl fun _ _ => mul_comm _ _
      _ ≤ ∑ _γ : Fin K, (K : ℤ) ^ (g + 1) + ∑ _q : Fin K × Fin K, 4 * (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2) :=
          add_le_add (Finset.sum_le_sum fun γ _ => hB1 γ) (Finset.sum_le_sum fun q _ => hB2 q)
      _ = _ := by
          rw [Finset.sum_const, Finset.sum_const, card_univ, card_univ, Fintype.card_prod, Fintype.card_fin,
            nsmul_eq_mul, nsmul_eq_mul]
          push_cast; ring
  have hA : (0 : ℤ) < (K : ℤ) ^ g := by positivity
  have e1 : (K : ℤ) * (K : ℤ) ^ (g + 1) = (K : ℤ) ^ 2 * (K : ℤ) ^ g := by ring
  have e2 : (K : ℤ) ^ 2 * (4 * (D : ℤ) ^ 2 * (K : ℤ) ^ (g + 2)) = 4 * (D : ℤ) ^ 2 * (K : ℤ) ^ 4 * (K : ℤ) ^ g := by ring
  rw [e1, e2] at hsum
  have hsum' : (g : ℤ) * (K : ℤ) ^ g ≤
      ((K : ℤ) ^ 2 + 4 * (D : ℤ) ^ 2 * (K : ℤ) ^ 4) * (K : ℤ) ^ g := by
    linarith
  have hbound := le_of_mul_le_mul_right hsum' hA
  have hthreshold : (K : ℤ) ^ 2 + 4 * (D : ℤ) ^ 2 * (K : ℤ) ^ 4 < g := by exact_mod_cast hg
  linarith

/-- A simple explicit size for the independent-accuracy guarantee. -/
theorem exists_good_accuracy_five_fourth (hK : 2 ≤ K)
    (c : Fin g → Fin g → ℕ) (hdiag : ∀ a, c a a = 0)
    (hsym : ∀ a b, c a b = c b a) (P D : ℕ) (hD : 1 ≤ D)
    (hrow : ∀ a, ∑ b, c a b ≤ P) (hg : g = 5 * D ^ 2 * K ^ 4) :
    ∃ χ : Fin g → Fin K, GoodAccuracy c P D χ := by
  apply exists_good_accuracy_of_threshold hK c hdiag hsym P D hrow
  have hK2 : 4 ≤ K ^ 2 :=
    le_trans (by norm_num : 4 ≤ 2 ^ 2) (Nat.pow_le_pow_left hK 2)
  have hK24 : K ^ 2 < K ^ 4 := by
    calc K ^ 2 < K ^ 2 * 2 := by omega
      _ ≤ K ^ 2 * K ^ 2 := Nat.mul_le_mul_left _ (by omega)
      _ = K ^ 4 := by ring
  have hD2 : 1 ≤ D ^ 2 := Nat.one_le_pow _ _ hD
  have hsmall : K ^ 2 < D ^ 2 * K ^ 4 :=
    lt_of_lt_of_le hK24 (Nat.le_mul_of_pos_left _ hD2)
  rw [hg]
  calc K ^ 2 + 4 * D ^ 2 * K ^ 4 < D ^ 2 * K ^ 4 + 4 * D ^ 2 * K ^ 4 :=
      Nat.add_lt_add_right hsmall _
    _ = 5 * D ^ 2 * K ^ 4 := by ring

/-- The original grouping guarantee at the sharper threshold. -/
theorem exists_good_of_threshold (hK : 2 ≤ K)
    (c : Fin g → Fin g → ℕ) (hdiag : ∀ a, c a a = 0)
    (hsym : ∀ a b, c a b = c b a) (P : ℕ)
    (hrow : ∀ a, ∑ b, c a b ≤ P) (hg : K ^ 2 + 4 * K ^ 6 < g) :
    ∃ χ : Fin g → Fin K, Good c P χ := by
  have hsize : K ^ 2 + 4 * K ^ 2 * K ^ 4 < g := by
    calc K ^ 2 + 4 * K ^ 2 * K ^ 4 = K ^ 2 + 4 * K ^ 6 := by ring
      _ < g := hg
  obtain ⟨χ, hχ⟩ := exists_good_accuracy_of_threshold hK c hdiag hsym P K hrow hsize
  refine ⟨χ, hχ.1, fun γ γ' hne => ?_⟩
  have h := hχ.2 γ γ' hne
  have he : (K : ℤ) * (K : ℤ) ^ 2 = (K : ℤ) ^ 3 := by ring
  simpa only [he] using h

/-- A simple sixth-power replacement for the original eighth-power size. -/
theorem exists_good_five_sixth (hK : 2 ≤ K)
    (c : Fin g → Fin g → ℕ) (hdiag : ∀ a, c a a = 0)
    (hsym : ∀ a b, c a b = c b a) (P : ℕ)
    (hrow : ∀ a, ∑ b, c a b ≤ P) (hg : g = 5 * K ^ 6) :
    ∃ χ : Fin g → Fin K, Good c P χ := by
  apply exists_good_of_threshold hK c hdiag hsym P hrow
  have hK4 : 2 ≤ K ^ 4 :=
    le_trans (by norm_num : 2 ≤ 2 ^ 4) (Nat.pow_le_pow_left hK 4)
  have hK2pos : 0 < K ^ 2 := by positivity
  have hK2lt : K ^ 2 < K ^ 6 := by
    calc K ^ 2 < K ^ 2 * 2 := by omega
      _ ≤ K ^ 2 * K ^ 4 := Nat.mul_le_mul_left _ hK4
      _ = K ^ 6 := by ring
  rw [hg]
  omega

/-- The smallest integer node count admitted by this union-bound certificate. -/
theorem exists_good_exact_threshold (hK : 2 ≤ K)
    (c : Fin g → Fin g → ℕ) (hdiag : ∀ a, c a a = 0)
    (hsym : ∀ a b, c a b = c b a) (P : ℕ)
    (hrow : ∀ a, ∑ b, c a b ≤ P)
    (hg : g = K ^ 2 + 4 * K ^ 6 + 1) :
    ∃ χ : Fin g → Fin K, Good c P χ := by
  apply exists_good_of_threshold hK c hdiag hsym P hrow
  omega

/-- Independent accuracy preserves the network coarsening interface, including
the equality of network values. -/
theorem network_grouping_accuracy_of_threshold {R : Type*} [CommRing R]
    (N : Net R g) (P D : ℕ) (hfit : N.Fits P) (hK : 2 ≤ K)
    (hg : K ^ 2 + 4 * D ^ 2 * K ^ 4 < g) :
    ∃ χ : Fin g → Fin K,
      (∀ γ, K * (univ.filter (χ · = γ)).card ≤ 2 * g) ∧
      (∀ γ γ', γ ≠ γ' → D * K ^ 2 * (N.coarsen χ).c γ γ' ≤ D * N.total + g * P) ∧
      ∀ L, (∀ γ γ', γ < γ' → (N.coarsen χ).c γ γ' ≤ L) →
        ((N.coarsen χ).toClique L).value = N.value := by
  obtain ⟨χ, h1, h2⟩ := exists_good_accuracy_of_threshold hK N.cc N.cc_self N.cc_symm P D hfit hg
  refine ⟨χ, fun γ => ?_, fun γ γ' hne => ?_, fun L hL => Net.value_group N χ L hL⟩
  · have := h1 γ
    have hs : grpSize χ γ = ((univ.filter (χ · = γ)).card : ℤ) := by
      simp [grpSize, Finset.sum_boole]
    rw [hs] at this
    exact_mod_cast this
  · have := h2 γ γ' hne
    rw [← Net.coarsen_c_eq N χ hne, Net.totU_cc] at this
    exact_mod_cast this

/-- The improved node threshold preserves the complete interface needed by
network coarsening: group sizes, pair widths, and equality of values. -/
theorem network_grouping_of_threshold {R : Type*} [CommRing R]
    (N : Net R g) (P : ℕ) (hfit : N.Fits P) (hK : 2 ≤ K)
    (hg : K ^ 2 + 4 * K ^ 6 < g) :
    ∃ χ : Fin g → Fin K,
      (∀ γ, K * (univ.filter (χ · = γ)).card ≤ 2 * g) ∧
      (∀ γ γ', γ ≠ γ' → K ^ 3 * (N.coarsen χ).c γ γ' ≤ K * N.total + g * P) ∧
      ∀ L, (∀ γ γ', γ < γ' → (N.coarsen χ).c γ γ' ≤ L) →
        ((N.coarsen χ).toClique L).value = N.value := by
  obtain ⟨χ, h1, h2⟩ := exists_good_of_threshold hK N.cc N.cc_self N.cc_symm P hfit hg
  refine ⟨χ, fun γ => ?_, fun γ γ' hne => ?_, fun L hL => Net.value_group N χ L hL⟩
  · have := h1 γ
    have hs : grpSize χ γ = ((univ.filter (χ · = γ)).card : ℤ) := by
      simp [grpSize, Finset.sum_boole]
    rw [hs] at this
    exact_mod_cast this
  · have := h2 γ γ' hne
    rw [← Net.coarsen_c_eq N χ hne, Net.totU_cc] at this
    exact_mod_cast this

end Factoring.Improvements

#print axioms Factoring.Improvements.exists_good_accuracy_of_threshold
#print axioms Factoring.Improvements.exists_good_accuracy_five_fourth
#print axioms Factoring.Improvements.network_grouping_accuracy_of_threshold
#print axioms Factoring.Improvements.exists_good_of_threshold
#print axioms Factoring.Improvements.exists_good_five_sixth
#print axioms Factoring.Improvements.exists_good_exact_threshold
#print axioms Factoring.Improvements.network_grouping_of_threshold
