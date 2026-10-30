/-
SPDX-License-Identifier: Apache-2.0

The closed theorem for the degree-twelve grouping variant. This module must
be compiled against the patched reference program, never against the old
schedule. Its statement retains the reference's exact machine, alphabet,
input/output encoding, and smallest-prime-factor definitions.
-/
module

public import FactoringK6.Quantitative
public import FactoringK6.MachineBridge
public import FactoringK6.RateInterpretation

@[expose] public section

namespace Factoring.K6

open Factoring.Aux
open IntegerFactoring.K6RateArithmetic
open IntegerFactoring.K6Quantitative
open IntegerFactoring.K6Machine
open IntegerFactoring.K6RateInterpretation

/-- One finite-alphabet machine, correct on every input, with an exact
integer bound exhibiting a thirtieth-root logarithmic saving. -/
theorem factoring_root_rate :
    ∃ M : FactoringMachine, ∃ C : ℕ,
      FactorsWithin M
        (fun n => C * 2 ^ (C * (n / (rateRoot n + 1) + 1))) := by
  have hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12 := by
    intro d
    exact (Kb_succ K0u d).le
  apply factoring_rate_of_Zu_bound (fun n => rateRoot n + 1) hstep
  exact Zu_rate_bound (by norm_num [K0u]) hstep (by norm_num [Kb, K0u])

/-- The same machine bound in the usual real-valued asymptotic form.
Both the algorithm and the conversion to exponent 1/30 are kernel checked. -/
theorem factoring_real_rate :
    ∃ M : FactoringMachine, ∃ C : ℝ,
      M.FiniteAlphabets ∧ 0 < C ∧ ∀ n : ℕ, 2 ≤ n →
        (M.time n : ℝ) ≤ C * Real.exp
          (C * (n : ℝ) / (Real.log (n : ℝ)) ^ (1 / 30 : ℝ)) := by
  obtain ⟨M, c, hfin, ht⟩ := factoring_root_rate
  obtain ⟨C, hC, htime⟩ := exists_real_rate_bound M.time c ht
  exact ⟨M, C, hfin, hC, htime⟩

/-- A strict improvement over the explicit exponent extracted from the
unchanged reference, together with the complete new machine bound. -/
theorem certified_improvement :
    Real.log (12 / 11 : ℝ) / Real.log 16 < (1 / 30 : ℝ) ∧
    ∃ M : FactoringMachine, ∃ C : ℝ,
      M.FiniteAlphabets ∧ 0 < C ∧ ∀ n : ℕ, 2 ≤ n →
        (M.time n : ℝ) ≤ C * Real.exp
          (C * (n : ℝ) / (Real.log (n : ℝ)) ^ (1 / 30 : ℝ)) :=
  ⟨baseline_log_rate_lt_one_thirtieth, factoring_real_rate⟩

end Factoring.K6

#print axioms Factoring.K6.factoring_root_rate
#print axioms Factoring.K6.factoring_real_rate
#print axioms Factoring.K6.certified_improvement
