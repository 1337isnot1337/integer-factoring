/-
SPDX-License-Identifier: Apache-2.0

The generic network-to-row conversion collapses to a single reference
clique entry. This avoids treating its formal finite sum as an algorithm.
-/
import CliqueRows

open Finset Groundwork Groundwork.Grid

namespace IntegerFactoring.HalfAppendix

theorem cliqueRows_eq {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L)
    (a : Fin K) (θ : Fin K → List Bool) :
    cliqueRows T a θ = T a (fun b => if b = a then Factoring.Aux.zeroBits L else bitsOf L (θ b)) := by
  classical
  let z0 : Fin K → Factoring.Aux.Bits L :=
    fun b => if b = a then Factoring.Aux.zeroBits L else bitsOf L (θ b)
  unfold cliqueRows Groundwork.Network.toInst
  simp only [cliqueNetwork, mem_univ, ite_true]
  rw [Finset.sum_eq_single z0]
  · have hp : (∏ b ∈ univ.erase a, if z0 b = bitsOf L (θ b) then (1 : ℤ) else 0) = 1 := by
      apply Finset.prod_eq_one
      intro b hb
      have hne : b ≠ a := (mem_erase.mp hb).1
      simp [z0, hne]
    simp only [show z0 a = Factoring.Aux.zeroBits L by simp [z0], ite_true, hp, mul_one]
    rfl
  · intro z _ hne
    by_cases hz : z a = Factoring.Aux.zeroBits L
    · have hdiff : ∃ b, b ≠ a ∧ z b ≠ bitsOf L (θ b) := by
        by_contra h
        push Not at h
        apply hne
        funext b
        by_cases hb : b = a
        · subst b; simpa [z0] using hz
        · simpa [z0, hb] using h b hb
      obtain ⟨b, hb, hzb⟩ := hdiff
      have hzero : (∏ b ∈ univ.erase a, if z b = bitsOf L (θ b) then (1 : ℤ) else 0) = 0 := by
        apply Finset.prod_eq_zero (Finset.mem_erase.mpr ⟨hb, mem_univ b⟩)
        simp [hzb]
      simp [hz, hzero]
    · simp [hz]
  · simp

/-- A bound on the reference entries transfers unchanged to encoded rows. -/
theorem cliqueQuery_bound {K L B : ℕ} (T : Factoring.Aux.Clique ℤ K L)
    (hT : ∀ a z, (T a z).natAbs ≤ 2 ^ B)
    (a : Fin K) (θ : Fin K → List Bool) :
    (cliqueQuery T (a.1, List.ofFn θ)).natAbs ≤ 2 ^ B := by
  rw [cliqueQuery_ofFn, cliqueRows_eq]
  exact hT _ _

/-- Direct row lookup, without enumerating the auxiliary clique-network states. -/
def cliqueQueryFast {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) (q : Groundwork.Prog.IQ) : ℤ :=
  if h : q.1 < K then
    T ⟨q.1, h⟩ (fun b => if b.1 = q.1 then Factoring.Aux.zeroBits L else bitsOf L (q.2.getD b.1 []))
  else 0

theorem cliqueQueryFast_eq {K L : ℕ} (T : Factoring.Aux.Clique ℤ K L) :
    cliqueQueryFast T = cliqueQuery T := by
  funext q
  unfold cliqueQueryFast cliqueQuery
  split
  · rw [cliqueRows_eq]
    congr 1
    funext b
    congr 1
    simp only [Fin.ext_iff]
  · rfl

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.cliqueRows_eq
#print axioms IntegerFactoring.HalfAppendix.cliqueQuery_bound
#print axioms IntegerFactoring.HalfAppendix.cliqueQueryFast_eq
