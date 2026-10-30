/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time interface for writing the explicit table, conditional on a
polynomial-time row query and polynomial-time parameter and pad preparation.
-/
import Tabulate
import PaddedDescent.F.CertNet

open Groundwork Groundwork.Prog PaddedDescent.F

namespace IntegerFactoring.HalfAppendix

theorem writtenInstance_polyT {α : Type} [Enc α]
    {pad : α → U} {K L : α → ℕ} {F : α → IQ → ℤ}
    (hp : PolyT pad) (hK : Unary K) (hL : Unary L)
    (hF : PolyT fun p : α × IQ => F p.1 p.2) :
    PolyT fun a => writtenInstance (pad a) (K a) (L a) (F a) := by
  have ht := PaddedDescent.F.CertNet.tabulate_polyT hp hK hL hF
  exact PolyT.prodMk hL ht

theorem writtenInstance_len (pad : U) (K L B : ℕ) (F : IQ → ℤ)
    (hpad : 2 ^ ((K - 1) * L) ≤ pad.length) (hL : 1 ≤ L)
    (hbound : ∀ (a : Fin K) (θ : Fin K → List Bool), WFQ L a θ →
      (F (a.1, List.ofFn θ)).natAbs ≤ 2 ^ B) :
    len (writtenInstance pad K L F) + 2 ≤ instLen K L B := by
  exact PaddedDescent.F.Table.stT3b K L B _ (writtenInstance_wfe pad K L B F hpad hL hbound).1

/-- Once the row query and both pads are supplied by polynomial-time
fragments, the written table and evaluator compose as one fragment. The
polynomial-time hypotheses are not proved here for the factoring constructor. -/
theorem written_evaluator_polyT {α : Type} [Enc α]
    {tabPad evalPad depth : α → U} {K L : α → ℕ} {F : α → IQ → ℤ}
    {run : U × U × PInst → ℤ}
    (htab : PolyT tabPad) (heval : PolyT evalPad) (hdepth : PolyT depth)
    (hK : Unary K) (hL : Unary L) (hF : PolyT fun p : α × IQ => F p.1 p.2)
    (hrun : PolyT run) :
    PolyT fun a => run (evalPad a, depth a, writtenInstance (tabPad a) (K a) (L a) (F a)) := by
  have hw := writtenInstance_polyT htab hK hL hF
  exact PolyT.comp hrun (PolyT.prodMk heval (PolyT.prodMk hdepth hw))

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.writtenInstance_polyT
#print axioms IntegerFactoring.HalfAppendix.writtenInstance_len
#print axioms IntegerFactoring.HalfAppendix.written_evaluator_polyT
