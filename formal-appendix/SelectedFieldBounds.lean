/-
SPDX-License-Identifier: Apache-2.0

The selected top-state bound suffices to decode the reference input and
residue-vector gate fields. Full gate checks and row summands remain open.
-/
import GateFields
import CodingBounds

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem selected_faP_field (s N k y P : ℕ) (hP : 2 ≤ P) (γ : Fin P)
    (cluster : Fin (3 * P ^ 8)) (r : Fin (pA s))
    (hc : selectedTopGrouping s N k y P hP cluster = γ)
    (hg : clMG s (3 * P ^ 8) (by positivity) (.faP r) = cluster)
    (w : Fin (ph s) → Fin (2 ^ pL s)) (i : Fin (ph s)) :
    fixedRadixDigit (topEntryBits s P) (ph s) (2 ^ pL s) (dCode s w).val i.val = (w i).val := by
  have hs : pD s ≤ 2 ^ topEntryBits s P := by
    simpa [gsz] using selected_gate_kind_size_bound s N k y P hP γ cluster (.faP r) hc hg
  exact fixed_input_digit_of_size s (topEntryBits s P) w i hs

theorem selected_faQ_field (s N k y P : ℕ) (hP : 2 ≤ P) (γ : Fin P)
    (cluster : Fin (3 * P ^ 8)) (r : Fin (pA s))
    (hc : selectedTopGrouping s N k y P hP cluster = γ)
    (hg : clMG s (3 * P ^ 8) (by positivity) (.faQ r) = cluster)
    (w : Fin (ph s) → Fin (2 ^ pL s)) (i : Fin (ph s)) :
    fixedRadixDigit (topEntryBits s P) (ph s) (2 ^ pL s) (dCode s w).val i.val = (w i).val := by
  have hs : pD s ≤ 2 ^ topEntryBits s P := by
    simpa [gsz] using selected_gate_kind_size_bound s N k y P hP γ cluster (.faQ r) hc hg
  exact fixed_input_digit_of_size s (topEntryBits s P) w i hs

theorem selected_ib_field (s N k y P : ℕ) (hP : 2 ≤ P) (γ : Fin P)
    (cluster : Fin (3 * P ^ 8)) (r : Fin (pA s))
    (hc : selectedTopGrouping s N k y P hP cluster = γ)
    (hg : clMG s (3 * P ^ 8) (by positivity) (.ib r) = cluster)
    (z : Fin (pA s) → ZMod (pM s)) (i : Fin (pA s)) :
    fixedRadixDigit (topEntryBits s P) (pA s) (pM s) (zCode s z).val i.val = (z i).val := by
  have hs : pZ s ≤ 2 ^ topEntryBits s P := by
    simpa [gsz] using selected_gate_kind_size_bound s N k y P hP γ cluster (.ib r) hc hg
  exact fixed_residue_digit_of_size s (topEntryBits s P) z i hs

theorem selected_mid_field (s N k y P : ℕ) (hP : 2 ≤ P) (γ : Fin P)
    (cluster : Fin (3 * P ^ 8)) (r : Fin (pA s))
    (hc : selectedTopGrouping s N k y P hP cluster = γ)
    (hg : clMG s (3 * P ^ 8) (by positivity) (.mid r) = cluster)
    (z : Fin (pA s) → ZMod (pM s)) (i : Fin (pA s)) :
    fixedRadixDigit (topEntryBits s P) (pA s) (pM s) (zCode s z).val i.val = (z i).val := by
  have hs : pZ s * pZ s ≤ 2 ^ topEntryBits s P := by
    simpa [gsz] using selected_gate_kind_size_bound s N k y P hP γ cluster (.mid r) hc hg
  have hz : pZ s ≤ 2 ^ topEntryBits s P := (Nat.le_mul_of_pos_right (pZ s) (pZ_pos s)).trans hs
  exact fixed_residue_digit_of_size s (topEntryBits s P) z i hz

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_faP_field
#print axioms IntegerFactoring.HalfAppendix.selected_faQ_field
#print axioms IntegerFactoring.HalfAppendix.selected_ib_field
#print axioms IntegerFactoring.HalfAppendix.selected_mid_field
