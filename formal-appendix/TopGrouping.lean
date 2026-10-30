/-
SPDX-License-Identifier: Apache-2.0

The factoring reference's exhaustive search selects a good top grouping.
This file establishes its combinatorial specification, not its cost.
-/
import FactoringAux.EngineGroup
import FactoringAux.Cluster

open Factoring.Aux

namespace IntegerFactoring.HalfAppendix

def topWidth (s K : ℕ) : ℕ :=
  (K * (10 * (pA s * pA s) * fw s) + (3 * K ^ 8) * clP s (3 * K ^ 8)) / K ^ 3

def topGroupingCode (s N k y K : ℕ) (hK : 2 ≤ K) : ℕ :=
  let g := 3 * K ^ 8
  let net := clusterNet s ℤ N k y g (by positivity)
  firstGood (okG net K (topWidth s K) (by omega)) (K ^ g) (K ^ g)

theorem topGroupingCode_good (s N k y K : ℕ) (hK : 2 ≤ K) :
    let g := 3 * K ^ 8
    let net := clusterNet s ℤ N k y g (by positivity)
    topGroupingCode s N k y K hK < K ^ g ∧
      okG net K (topWidth s K) (by omega) (topGroupingCode s N k y K hK) := by
  dsimp only
  let g := 3 * K ^ 8
  have hg : 0 < g := by dsimp [g]; positivity
  let net := clusterNet s ℤ N k y g hg
  obtain ⟨q, hq, hok⟩ := good_exists net (clP s g) (10 * (pA s * pA s) * fw s)
    (clusterNet_fits s ℤ N k y g hg) (clusterNet_total_le s ℤ N k y g hg) hK rfl
  have he : ∃ q < K ^ g, okG net K (topWidth s K) (by omega) q := by
    refine ⟨q, hq, ?_⟩
    simpa [topWidth, g] using hok
  simpa [topGroupingCode, g, net] using
    (firstGood_spec (okG net K (topWidth s K) (by omega)) (le_refl (K ^ g)) he)

def selectedTopGrouping (s N k y K : ℕ) (hK : 2 ≤ K) : Fin (3 * K ^ 8) → Fin K :=
  chiOf K (3 * K ^ 8) (topGroupingCode s N k y K hK) (by omega)

theorem selectedTopGrouping_fit (s N k y K : ℕ) (hK : 2 ≤ K) :
    ∀ a b, a < b →
      ((clusterNet s ℤ N k y (3 * K ^ 8) (by positivity)).coarsen
        (selectedTopGrouping s N k y K hK)).c a b ≤ topWidth s K + 1 := by
  have h := (topGroupingCode_good s N k y K hK).2.2
  intro a b hab
  exact (h a b hab).trans (Nat.le_succ _)

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.topGroupingCode_good
#print axioms IntegerFactoring.HalfAppendix.selectedTopGrouping_fit
