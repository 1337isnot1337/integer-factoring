/-
SPDX-License-Identifier: Apache-2.0

Bridge the explicit pad arithmetic to the pinned Appendix B evaluator.
This only produces a padded clique evaluator, not a factoring machine.
-/
import GeneralPieces.Main
import HalfPad

open Groundwork PaddedDescent.F Levels.Gen ClauseBound.Alg
open IntegerFactoring.HalfPad

set_option maxHeartbeats 2000000

namespace IntegerFactoring.HalfAppendix

theorem K_pos (x l : ℕ) : 0 < K x l := by
  induction l with
  | zero => simp [K]
  | succ l ih => simp only [K]; positivity

theorem K_two (x l : ℕ) (hx : 1 ≤ x) : 2 ≤ K x l := by
  induction l with
  | zero =>
      simpa [K] using (Nat.pow_le_pow_right (by norm_num : 0 < 2) hx)
  | succ l ih =>
      simp only [K]
      have hb : 1 ≤ 2 ^ (3 * l + x) := Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
      have hh : 2 ≤ 2 ^ (3 * l + x) * K x l := by nlinarith
      nlinarith

/-- Replace the existential pad bound in the source interface. -/
theorem explicit_of_parts (x C : ℕ)
    (ha : AppRunningTime.PartA (1 / Real.sqrt 2) x)
    (hb : AppRunningTime.PartB (1 / Real.sqrt 2) x C) :
    ∃ (run : U × U × PInst → ℤ), PolyT run ∧
      ∀ (l lam B : ℕ) (I : PInst), 1 ≤ lam → 1 ≤ B → WFe (K x l) lam B I →
        run (List.replicate (2 ^ padExponent (K x l) lam B C (2 ^ (l / 2))) (),
          List.replicate l (), I) = zval (K x l) I := by
  obtain ⟨run, need, hpoly, hcorrect, hneed⟩ := hb
  refine ⟨run, hpoly, ?_⟩
  intro l lam B I hlam hB hI
  have harg : 0 < K x l * (lam + B) := by have := K_pos x l; positivity
  have hpad := half_need_le_padExponent (K x l) lam B C l (need l lam B)
    harg (e (1 / Real.sqrt 2) x l) (ha.e_bounds l).2 (hneed l lam B hlam hB)
  apply hcorrect _ l lam B I hlam hB hI
  simpa using Nat.pow_le_pow_right (by norm_num : 0 < 2) hpad

/-- A single Appendix B evaluator is correct with a computable pad exponent.
The theorem does not bound the cost of constructing the input clique or
translate this fragment into the factoring reference's machine. -/
theorem exists_explicit_padded_evaluator_ge :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l lam B : ℕ) (I : PInst), 1 ≤ lam → 1 ≤ B → WFe (K x l) lam B I →
        run (List.replicate (2 ^ padExponent (K x l) lam B C (2 ^ (l / 2))) (),
          List.replicate l (), I) = zval (K x l) I := by
  obtain ⟨c, P, A, val, hp, _⟩ := AppRunningTime.thmHalfPieces
  obtain ⟨C, hC⟩ := AppRunningTime.thmRunningTime.b c
  let x := c + 23
  have hx : 18 ≤ x := by dsimp [x]; omega
  have hpow : c + 5 ≤ 2 ^ x := by
    have h := x.lt_two_pow_self
    dsimp [x] at *
    omega
  have hex : AppRunningTime.Exist (1 / Real.sqrt 2) c := ⟨P, A, val, hp⟩
  have ha := AppRunningTime.thmRunningTime.a _ x hp.rho_ge hp.rho_lt hx
  obtain ⟨run, hpoly, hcorrect⟩ := explicit_of_parts x C ha (hC _ x hex hx hpow)
  exact ⟨x, C, run, hx, hpoly, hcorrect⟩

theorem exists_explicit_padded_evaluator :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), PolyT run ∧
      ∀ (l lam B : ℕ) (I : PInst), 1 ≤ lam → 1 ≤ B → WFe (K x l) lam B I →
        run (List.replicate (2 ^ padExponent (K x l) lam B C (2 ^ (l / 2))) (),
          List.replicate l (), I) = zval (K x l) I := by
  obtain ⟨x, C, run, _, hpoly, hcorrect⟩ := exists_explicit_padded_evaluator_ge
  exact ⟨x, C, run, hpoly, hcorrect⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.exists_explicit_padded_evaluator
