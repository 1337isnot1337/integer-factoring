/-
SPDX-License-Identifier: Apache-2.0

A quantitative bound on the actual uniform engine exponent. The two
schedule hypotheses are discharged by the modified implementation in the
final theorem. The integer denominator is positive for every input length.
-/
module

public import FactoringK6.BudgetBounds
public import FactoringK6.ChainEnvelope
public import FactoringK6.BottomException
public import FactoringK6.RateArithmetic

@[expose] public section

namespace IntegerFactoring.K6Quantitative

open Factoring.Aux Factoring.K6
open IntegerFactoring.K6Budget IntegerFactoring.K6RateArithmetic

def totalCoefficient (c : ℕ) : ℕ :=
  topCoefficient + 40000 + (25 * c) * volumeCoefficient + 9623

theorem weighted_zC_le {ell n R c : ℕ} (hK : 48 ≤ K0u)
    (h : UGood ell n) (hR : R ≤ Kb K0u ell)
    (hdecay : R * 23 ^ ell ≤ c * 25 ^ ell) :
    R * zC ell n ≤ ((25 * c) * volumeCoefficient + 9623) * n := by
  have hchain := Nat.mul_le_mul_left R (zC_le_envelope48 hK ell n)
  have hbottom := weighted_bottomBudget_of_full_depth48 K0u ell n R c h.1 hdecay
  have htop := Nat.mul_le_mul_left (25 * c) (top_volume_le h)
  have hdepth := weighted_depth_padding_le h hR
  have hpad := weighted_chain_padding_le h hR
  have hG : Gs K0u ell ^ 4 = Gm4 K0u ell := by unfold Gm4; ring
  unfold chainEnvelope at hchain
  rw [hG] at hchain
  nlinarith only [hchain, hbottom, htop, hdepth, hpad]

theorem weighted_Zof_le {ell n R c : ℕ} (hK : 48 ≤ K0u)
    (h : UGood ell n) (hR : R ≤ Kb K0u ell)
    (hdecay : R * 23 ^ ell ≤ c * 25 ^ ell) :
    R * Zof ell n ≤ totalCoefficient c * n := by
  have ht := weighted_zT_le h hR
  have hc := weighted_zC_le hK h hR hdecay
  have hp := weighted_zP_le h hR
  unfold Zof totalCoefficient
  nlinarith only [ht, hc, hp]

theorem denominator_Zof_le {ell n : ℕ} (hK : 48 ≤ K0u)
    (hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12)
    (hbase : Kb K0u 0 ≤ 1024) (h : UGood ell n) :
    (rateRoot n + 1) * Zof ell n ≤ totalCoefficient 3 * n := by
  have hR : rateRoot n + 1 ≤ Kb K0u ell :=
    rateRoot_add_one_le_parameter (Kb48_ge48 ell) (good_facts h).2.2.2.1
  have hw := root_weighted_bound_small_start (Kb K0u) n ell hstep hbase h.2.2
  have hd := denominator_weighted_bound n ell 2 hw
  exact weighted_Zof_le hK h hR hd

/-- Division introduces only a fixed additive unit, uniformly in the divisor. -/
theorem quotient_bound {Z n D c : ℕ} (hD : 0 < D) (hZ : D * Z ≤ c * n) :
    Z ≤ c * (n / D + 1) := by
  have hn : n ≤ D * (n / D + 1) := by
    have hmod := Nat.mod_lt n hD
    have hdiv := Nat.div_add_mod n D
    nlinarith only [hmod, hdiv]
  have hprod : D * Z ≤ D * (c * (n / D + 1)) := by
    calc D * Z ≤ c * n := hZ
      _ ≤ c * (D * (n / D + 1)) := Nat.mul_le_mul_left _ hn
      _ = _ := by ring
  exact Nat.le_of_mul_le_mul_left hprod hD

/-- Global quantitative bound, including the finite range at depth zero. -/
theorem Zu_rate_bound (hK : 48 ≤ K0u)
    (hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12)
    (hbase : Kb K0u 0 ≤ 1024) :
    ∃ A : ℕ, ∀ n, Zu n ≤ A * (n / (rateRoot n + 1) + 1) := by
  let B : ℕ := ∑ m ∈ Finset.range (Kb K0u 1 ^ 40), Zu m
  refine ⟨B + totalCoefficient 3, fun n => ?_⟩
  by_cases hell : 1 ≤ ellOf n
  · have h := denominator_Zof_le hK hstep hbase (good_ellOf hell)
    have hz : Zu n ≤ totalCoefficient 3 * (n / (rateRoot n + 1) + 1) :=
      quotient_bound (by omega) h
    exact le_trans hz (Nat.mul_le_mul_right _ (by omega))
  · have he : ellOf n = 0 := by omega
    have hn : n < Kb K0u 1 ^ 40 := by
      have h := lt_Kb_ellOf n
      simpa [he] using h
    have hb : Zu n ≤ B := Finset.single_le_sum
      (f := fun m => Zu m) (fun _ _ => Nat.zero_le _) (Finset.mem_range.mpr hn)
    have hA : B ≤ B + totalCoefficient 3 := Nat.le_add_right _ _
    have hmul : B + totalCoefficient 3 ≤
        (B + totalCoefficient 3) * (n / (rateRoot n + 1) + 1) :=
      Nat.le_mul_of_pos_right _ (Nat.succ_pos _)
    exact le_trans hb (le_trans hA hmul)

end IntegerFactoring.K6Quantitative

#print axioms IntegerFactoring.K6Quantitative.weighted_Zof_le
#print axioms IntegerFactoring.K6Quantitative.denominator_Zof_le
#print axioms IntegerFactoring.K6Quantitative.Zu_rate_bound
