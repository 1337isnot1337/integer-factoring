/-
SPDX-License-Identifier: Apache-2.0

Extensional factor correctness for the explicit top-evaluator composition.
This does not prove that the composition is a machine or has any time bound.
-/
import PrefixEvaluator

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

def prefixFromEvaluator (x C : ℕ) (run : U × U × PInst → ℤ)
    (l n N : ℕ) (hK : 2 ≤ K x l) (k y : ℕ) : ℕ :=
  let s := sOf n
  let P := K x l
  let L := topWidth s P + 1
  let B := topEntryBits s P
  let χ := selectedTopGrouping s N k y P hK
  let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
  let tabPad := List.replicate (2 ^ ((P - 1) * L)) ()
  (run (List.replicate (2 ^ padExponent P L B C (2 ^ (l / 2))) (),
    List.replicate l (), writtenInstance tabPad P L (cliqueQueryFast T)) -
      (if 1 / 2 ^ k = y then 1 else 0)).toNat

/-- The same evaluator, composed extensionally with the reference search,
returns `Nat.minFac`. The choice of a level and all computational costs are
outside this statement. -/
theorem pure_factor_correct :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l n N : ℕ) (hK : 2 ≤ K x l), 2 ≤ N → N < 2 ^ n →
        search (prefixFromEvaluator x C run l n N hK) n = N.minFac := by
  obtain ⟨x, C, run, hx, hpoly, heval⟩ := evaluator_counts_reference_prefix
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l n N hK hN hNbits
  have hcount : prefixFromEvaluator x C run l n N hK = cntPre N := by
    funext k y
    unfold prefixFromEvaluator
    dsimp only
    rw [heval l n N k y hK hN hNbits]
    simp
  rw [hcount]
  exact search_eq_minFac hN hNbits

def pureFactor (x C : ℕ) (run : U × U × PInst → ℤ)
    (l : ℕ) (hK : 2 ≤ K x l) (N : ℕ) : ℕ :=
  if N = 0 then 2 else if N = 1 then 1 else
    search (prefixFromEvaluator x C run l N.size N hK) N.size

/-- The same value-level composition also handles the two special inputs.
There is still no machine or time certificate for this function. -/
theorem pure_factor_all :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l N : ℕ) (hK : 2 ≤ K x l), pureFactor x C run l hK N = N.minFac := by
  obtain ⟨x, C, run, hx, hpoly, hfactor⟩ := pure_factor_correct
  refine ⟨x, C, run, hx, hpoly, ?_⟩
  intro l N hK
  by_cases h0 : N = 0
  · subst N; simp [pureFactor, Nat.minFac_zero]
  by_cases h1 : N = 1
  · subst N; simp [pureFactor, Nat.minFac_one]
  simp only [pureFactor, h0, h1, ite_false]
  exact hfactor l N.size N hK (by omega) (Nat.lt_size_self N)

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.pure_factor_correct
#print axioms IntegerFactoring.HalfAppendix.pure_factor_all
