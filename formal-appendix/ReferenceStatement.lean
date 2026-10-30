/-
SPDX-License-Identifier: Apache-2.0

Connect value-level correctness to the unchanged factoring challenge's
audited output specification. The source challenge also contains an
unrelated placeholder theorem; this theorem does not depend on it.
-/
import PureFactor
import Challenge

open Groundwork PaddedDescent.F Levels.Gen

namespace IntegerFactoring.HalfAppendix

theorem pure_factor_reference_spec :
    ∃ (x C : ℕ) (run : U × U × PInst → ℤ), 18 ≤ x ∧ PolyT run ∧
      ∀ (l N : ℕ) (hK : 2 ≤ K x l), pureFactor x C run l hK N = Factoring.factor N := by
  simpa only [Factoring.factor] using pure_factor_all

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.pure_factor_reference_spec
