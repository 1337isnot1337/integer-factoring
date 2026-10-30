/-
SPDX-License-Identifier: Apache-2.0

A binary clique summand is one exactly when its weight and every factor are one.
-/
import SelectedEntry
import RowSumBridge

open Factoring.Aux Finset

namespace IntegerFactoring.HalfAppendix

private theorem binary_prod_eq_one {ι : Type} [DecidableEq ι] (s : Finset ι) (f : ι → ℤ)
    (hf : ∀ i ∈ s, Is01 (f i)) : (∏ i ∈ s, f i) = 1 ↔ ∀ i ∈ s, f i = 1 := by
  constructor
  · intro h i hi
    rcases hf i hi with h0 | h1
    · have hz : (∏ i ∈ s, f i) = 0 := Finset.prod_eq_zero hi h0
      omega
    · exact h1
  · exact Finset.prod_eq_one

theorem cliqueTerm_01 {K L : ℕ} (N : Net ℤ K) (hN : Net01 N) (γ : Fin K)
    (θ : Fin K → Bits L) (z : N.S γ) : Is01 (cliqueTerm N γ θ z) := by
  unfold cliqueTerm
  apply is01_mul (is01_mul (hN.weights γ z) ?_) ?_
  · apply prod_01
    intro δ _
    exact (padF_01 N hN γ δ).1 z (θ δ)
  · apply prod_01
    intro δ _
    exact (padF_01 N hN δ γ).2 z (θ δ)

theorem cliqueTerm_one_iff {K L : ℕ} (N : Net ℤ K) (hN : Net01 N) (γ : Fin K)
    (θ : Fin K → Bits L) (z : N.S γ) : cliqueTerm N γ θ z = 1 ↔
      N.W γ z = 1 ∧
      (∀ δ : Fin K, γ < δ → (N.padF L γ δ).l z (θ δ) = 1) ∧
      (∀ δ : Fin K, δ < γ → (N.padF L δ γ).r z (θ δ) = 1) := by
  let U := (Finset.univ : Finset (Fin K)).filter (γ < ·)
  let D := (Finset.univ : Finset (Fin K)).filter (· < γ)
  let f := fun δ => (N.padF L γ δ).l z (θ δ)
  let h := fun δ => (N.padF L δ γ).r z (θ δ)
  have hU : ∀ δ ∈ U, Is01 (f δ) := by intro δ _; exact (padF_01 N hN γ δ).1 z (θ δ)
  have hD : ∀ δ ∈ D, Is01 (h δ) := by intro δ _; exact (padF_01 N hN δ γ).2 z (θ δ)
  have hW := hN.weights γ z
  have hPU := prod_01 U f hU
  have hPD := prod_01 D h hD
  have hmul : N.W γ z * (∏ δ ∈ U, f δ) * (∏ δ ∈ D, h δ) = 1 ↔
      N.W γ z = 1 ∧ (∏ δ ∈ U, f δ) = 1 ∧ (∏ δ ∈ D, h δ) = 1 := by
    rcases hW with h0 | h1 <;> rcases hPU with h2 | h3 <;> rcases hPD with h4 | h5 <;> simp_all
  change (N.W γ z * (∏ δ ∈ U, f δ) * (∏ δ ∈ D, h δ) = 1) ↔ _
  rw [hmul, binary_prod_eq_one U f hU, binary_prod_eq_one D h hD]
  simp [U, D, f, h]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.cliqueTerm_01
#print axioms IntegerFactoring.HalfAppendix.cliqueTerm_one_iff
