/-
SPDX-License-Identifier: Apache-2.0

The written-table interface used by Appendix B. The rows encode only
off-diagonal labels; the source's tabulation theorem supplies the exact
lookup and width facts. This does not implement the factoring-network table.
-/
import PaddedDescent.F.Table
import ClauseBound.Alg.Interfaces

open Groundwork Groundwork.Prog PaddedDescent.F ClauseBound.Alg

namespace IntegerFactoring.HalfAppendix

def writtenInstance (pad : U) (K lam : ℕ) (F : IQ → ℤ) : PInst :=
  (List.replicate lam (), tabulate pad K lam F)

theorem writtenInstance_value (pad : U) (K lam : ℕ) (F : IQ → ℤ)
    (hpad : 2 ^ ((K - 1) * lam) ≤ pad.length) :
    zval K (writtenInstance pad K lam F) =
      Inst.ZW lam [] (fun (v : Fin K) θ => F (v.1, List.ofFn θ)) := by
  unfold zval writtenInstance
  simp only [List.length_replicate]
  apply Table.stT1 K lam
  intro v θ hθ
  exact Table.stT2 pad K lam F hpad v θ hθ

theorem writtenInstance_wfe (pad : U) (K lam B : ℕ) (F : IQ → ℤ)
    (hpad : 2 ^ ((K - 1) * lam) ≤ pad.length) (hlam : 1 ≤ lam)
    (hbound : ∀ (v : Fin K) (θ : Fin K → List Bool), WFQ lam v θ →
      (F (v.1, List.ofFn θ)).natAbs ≤ 2 ^ B) :
    WFe K lam B (writtenInstance pad K lam F) := by
  refine ⟨Table.stT3a pad K lam lam B F hpad hlam le_rfl hbound, ?_⟩
  simp [writtenInstance]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.writtenInstance_value
#print axioms IntegerFactoring.HalfAppendix.writtenInstance_wfe
