/-
SPDX-License-Identifier: Apache-2.0

Correctness and step bound for the Appendix evaluator on the selected,
already-written factoring table. Constructing the table and the complete
factoring machine remain outside this theorem.
-/
import PrefixEvaluator
import EvaluatorSteps
import BoundedLength

open Groundwork Groundwork.Prog Groundwork.TM PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux
open IntegerFactoring.HalfPad

namespace IntegerFactoring.HalfAppendix

theorem prefix_evaluator_correct_steps :
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
        2 ≤ N → N < 2 ^ n → P ^ 40 ≤ n → n ≤ (K x (l + 1)) ^ 40 →
        (run input - (if 1 / 2 ^ k = y then 1 else 0) = (cntPre N k y : ℤ)) ∧
        ∀ (M S T' N' : List Sym), ∃ t ≤ 2 ^ (D * (n / r + 1)),
          Pgm.Runs none (st (wd input ++ M) S T' N') t (st (wd (run input) ++ M) S T' N') := by
  obtain ⟨x, C, run, hx, hpoly, hcorrect⟩ := evaluator_counts_reference_prefix
  obtain ⟨Pgm, D, hsteps⟩ := evaluator_steps_at_schedule (1 / Real.sqrt 2) x C run hpoly
  obtain ⟨c, P', A', val, hp, _⟩ := AppRunningTime.thmHalfPieces
  have ha := AppRunningTime.thmRunningTime.a _ x hp.rho_ge hp.rho_lt hx
  refine ⟨x, C, run, Pgm, D, hx, ?_⟩
  intro l n N k y hK
  dsimp only
  intro hN hNbits hnlo hnhi
  constructor
  · exact hcorrect l n N k y hK hN hNbits
  · apply hsteps l n _ (by omega) ha hnlo hnhi
    exact selected_instance_len (sOf n) N k y (K x l) hK

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.prefix_evaluator_correct_steps
