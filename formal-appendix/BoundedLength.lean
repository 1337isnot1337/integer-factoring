/-
SPDX-License-Identifier: Apache-2.0

The selected table satisfies the Appendix instance-length envelope. This
does not establish the time required to construct the table.
-/
import BoundedTop
import WrittenPoly

open Groundwork Groundwork.Prog PaddedDescent.F Levels.Gen ClauseBound.Alg Factoring.Aux

namespace IntegerFactoring.HalfAppendix

theorem selected_instance_len (s N k y P : ℕ) (hP : 2 ≤ P) :
    let L := topWidth s P + 1
    let B := topEntryBits s P
    let χ := selectedTopGrouping s N k y P hP
    let T := ((clusterNet s ℤ N k y (3 * P ^ 8) (by positivity)).coarsen χ).toClique L
    let pad := List.replicate (2 ^ ((P - 1) * L)) ()
    len (writtenInstance pad P L (cliqueQueryFast T)) + 2 ≤ PaddedDescent.F.instLen P L B := by
  dsimp only
  apply writtenInstance_len
  · simp
  · omega
  · intro a θ hθ
    rw [cliqueQueryFast_eq]
    apply cliqueQuery_bound
    intro a z
    exact (selectedTop_entry_bound s N k y P (topWidth s P + 1) hP a z).trans
      (Nat.pow_le_pow_right (by norm_num) (by simp [topEntryBits]))

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.selected_instance_len
