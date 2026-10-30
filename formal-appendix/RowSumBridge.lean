/-
SPDX-License-Identifier: Apache-2.0

Value bridge from padded numeric state enumeration to a reference clique row.
Computability of the supplied state codec and summand remains an obligation.
-/
import PaddedRowSum
import FactoringAux.Clique

open Factoring.Aux Finset

namespace IntegerFactoring.HalfAppendix

noncomputable def cliqueTerm {K L : ℕ} (N : Net ℤ K) (γ : Fin K)
    (θ : Fin K → Bits L) (z : N.S γ) : ℤ :=
  N.W γ z * (∏ γ' ∈ univ.filter (γ < ·), (N.padF L γ γ').l z (θ γ')) *
    ∏ γ' ∈ univ.filter (· < γ), (N.padF L γ' γ).r z (θ γ')

noncomputable def cliqueSummand {K L S : ℕ} (N : Net ℤ K) (γ : Fin K)
    (θ : Fin K → Bits L) (codec : N.S γ ≃ Fin S) (i : ℕ) : ℤ :=
  if h : i < S then cliqueTerm N γ θ (codec.symm ⟨i, h⟩) else 0

theorem paddedRowSumZ_eq_clique {K L S : ℕ} (N : Net ℤ K) (γ : Fin K)
    (θ : Fin K → Bits L) (codec : N.S γ ≃ Fin S) (pad : List Unit) (hpad : S ≤ pad.length) :
    paddedRowSumZ S pad (cliqueSummand N γ θ codec) = N.toClique L γ θ := by
  classical
  rw [paddedRowSumZ_eq_finSum S pad _ hpad]
  unfold Net.toClique
  apply Fintype.sum_equiv codec.symm
  intro i
  simp [cliqueSummand, cliqueTerm, i.isLt]

/-- A padded range may overestimate the number of states: invalid codes
    contribute zero. The summand is still not certified executable here. -/
theorem paddedRowSumZ_bound_eq_clique {K L S T : ℕ} (N : Net ℤ K) (γ : Fin K)
    (θ : Fin K → Bits L) (codec : N.S γ ≃ Fin S) (pad : List Unit)
    (hST : S ≤ T) (hpad : T ≤ pad.length) :
    paddedRowSumZ T pad (cliqueSummand N γ θ codec) = N.toClique L γ θ := by
  classical
  have hsmall : S ≤ pad.length := hST.trans hpad
  have heq : paddedRowSumZ T pad (cliqueSummand N γ θ codec) =
      paddedRowSumZ S pad (cliqueSummand N γ θ codec) := by
    rw [paddedRowSumZ_eq_finSum T pad _ hpad, paddedRowSumZ_eq_finSum S pad _ hsmall]
    rw [Fin.sum_univ_eq_sum_range, Fin.sum_univ_eq_sum_range]
    symm
    apply Finset.sum_subset (Finset.range_mono hST)
    intro i _ hi
    have hn : ¬i < S := by simpa using hi
    simp [cliqueSummand, hn]
  rw [heq]
  exact paddedRowSumZ_eq_clique N γ θ codec pad hsmall

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.paddedRowSumZ_eq_clique
#print axioms IntegerFactoring.HalfAppendix.paddedRowSumZ_bound_eq_clique
