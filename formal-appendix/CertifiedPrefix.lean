/-
SPDX-License-Identifier: Apache-2.0

Correctness and square-root-logarithmic exponent for the evaluator fragment
on the already-written top table. The table constructor and complete
factoring-machine runtime remain outside this theorem.
-/
import PrefixSteps
import EvaluatorRate

open Groundwork Groundwork.Prog Groundwork.TM PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem prefix_evaluator_sqrt_steps :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ) (Pgm : Fr) (D : ℕ), 18 ≤ x ∧
      ∀ (l n N k y : ℕ) (hK : 2 ≤ K x l),
        let s := sOf n
        let P := K x l
        let L := topWidth s P + 1
        let B := topEntryBits s P
        let χ := selectedTopGrouping s N k y P hK
        let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
        let tabPad := List.replicate (2 ^ ((P - 1) * L)) ()
        let I := writtenInstance tabPad P L (cliqueQueryFast T)
        let r := 2 ^ (l / 2)
        let input := (List.replicate (2 ^ padExponent P L B C r) (), List.replicate l (), I)
        1 < n → 2 ≤ N → N < 2 ^ n → P ^ 40 ≤ n → n ≤ (K x (l + 1)) ^ 40 →
        (run input - (if 1 / 2 ^ k = y then 1 else 0) = (cntPre N k y : ℤ)) ∧
        (((D * (n / r + 1) : ℕ) : ℝ) ≤
          (D : ℝ) * Real.sqrt (160 * (3 * (x : ℝ) + 6) * Real.log 2) *
            (n : ℝ) / Real.sqrt (Real.log n) + D) ∧
        ∀ (M S T' N' : List Sym), ∃ t ≤ 2 ^ (D * (n / r + 1)),
          Pgm.Runs none (st (wd input ++ M) S T' N') t (st (wd (run input) ++ M) S T' N') := by
  obtain ⟨x, C, run, Pgm, D, hx, hpref⟩ := prefix_evaluator_correct_steps
  obtain ⟨c, P', A', val, hp, _⟩ := AppRunningTime.thmHalfPieces
  have ha := AppRunningTime.thmRunningTime.a _ x hp.rho_ge hp.rho_lt hx
  refine ⟨x, C, run, Pgm, D, hx, ?_⟩
  intro l n N k y hK
  dsimp only
  intro hn hN hNbits hnlo hnhi
  have h := hpref l n N k y hK hN hNbits hnlo hnhi
  exact ⟨h.1, evaluator_step_exponent_sqrt (1 / Real.sqrt 2) x l n D ha hn hnhi, h.2⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.prefix_evaluator_sqrt_steps
