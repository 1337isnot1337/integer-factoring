/-
SPDX-License-Identifier: Apache-2.0

An all-level exponent envelope with contraction 23/25, valid for any base
group count at least 600. The interfaces are generic in Kb/Gs and therefore
also apply to the k6 variant once its source libraries have been rebuilt.

Based on the checked chain argument in the pinned upstream repository.
No asymptotic conversion or final machine theorem is asserted in this file.
-/
module

public import FactoringAux.UniExp

@[expose] public section

namespace Factoring.K6

open Factoring.Aux

/-- The generic group-count recurrence grows by at least a factor of two. -/
theorem kb_ge_double (K0 : ℕ) (hK0 : 2 ≤ K0) : ∀ d,
    K0 * 2 ^ d ≤ Kb K0 d
  | 0 => by simp [Kb]
  | d + 1 => by
      have ih := kb_ge_double K0 hK0 d
      have hk := Kb_ge K0 hK0 d
      have hg : Kb K0 d ≤ Gs K0 d := K_le_G _
      have hg2 : 2 * Kb K0 d ≤ Gs K0 d * Gs K0 d :=
        le_trans (Nat.mul_le_mul_right _ hk) (Nat.mul_le_mul hg hg)
      change K0 * 2 ^ (d + 1) ≤ Gs K0 d * Gs K0 d
      rw [pow_succ]
      nlinarith

/-- The source's exact width formula contracts by 23/25 once K >= 600. -/
theorem N_step_23_25 (K G lam : ℕ) (hK : 600 ≤ K) (hKG : K ≤ G) :
    25 * (K ^ 2 * Lnext K G lam) ≤
      23 * (G * G * (G * G) * lam) + 50 * (G * G * (G * G)) := by
  have h := N_next K G lam hKG
  have e : G ^ 4 = G * G * (G * G) := by ring
  rw [e] at h
  set N := G * G * (G * G) * lam
  set N' := K ^ 2 * Lnext K G lam
  set G4 := G * G * (G * G)
  have hkN : 600 * N ≤ K * N := Nat.mul_le_mul_right _ hK
  have hs := Nat.mul_le_mul_left 25 h
  have hc : (12 * K) * (25 * N') ≤ (12 * K) * (23 * N + 50 * G4) := by
    nlinarith
  exact Nat.le_of_mul_le_mul_left hc (by positivity)

/-- The source chain estimate with the sharper 23/25 contraction. -/
theorem chainMax_le_23_25 (K0 Y mb W W0 Gm : ℕ) (hK : 600 ≤ K0)
    (hLev : ∀ d N, N + 2 * (Gm * Gm * (Gm * Gm)) ≤ W →
      23 ^ (d + 1) * N ≤ 25 ^ (d + 1) * W0 →
      14 * (N + 2 * (Gm * Gm * (Gm * Gm))) / Kb K0 d + 7 +
        (Gm * Gm * (Gm * Gm) + 9) *
          (Nat.log 2 (W + mb + Gm + Gm + 64) + 1) ≤ Y)
    (hBr : W0 + 2 * (Nat.log 2 (mb + K0 + 64) + 1) ≤ Y) :
    ∀ d L, Kb K0 d ≤ Gm →
      Kb K0 d ^ 2 * L + 2 * d * (Gm * Gm * (Gm * Gm)) ≤ W →
      23 ^ d * (Kb K0 d ^ 2 * L) +
        25 ^ d * (2 * d * (Gm * Gm * (Gm * Gm))) ≤ 25 ^ d * W0 →
      chainMax mb K0 d L ≤ Y
  | 0, L, _, _, h2 => by
      show L * (K0 * K0) + 2 * (Nat.log 2 (mb + K0 + 64) + 1) ≤ Y
      have he : Kb K0 0 ^ 2 * L = L * (K0 * K0) := by show K0 ^ 2 * L = _; ring
      simp only [pow_zero, one_mul, mul_zero, zero_mul, add_zero, he] at h2
      omega
  | d + 1, L, hGm, h1, h2 => by
      show max (eLev (Kb K0 d) (Gs K0 d) L mb)
        (chainMax mb K0 d (Lnext (Kb K0 d) (Gs K0 d) L)) ≤ Y
      set K := Kb K0 d
      set G := Gs K0 d
      have hKd : K0 ≤ K := Kb_mono K0 (by omega) (Nat.zero_le d)
      have hKG : K ≤ G := K_le_G K
      have hGGm : G ≤ Gm := le_trans (show G ≤ Kb K0 (d + 1) by
        show G ≤ G * G
        exact Nat.le_mul_of_pos_left _ (by omega)) hGm
      have hG4 : G * G * (G * G) ≤ Gm * Gm * (Gm * Gm) :=
        Nat.mul_le_mul (Nat.mul_le_mul hGGm hGGm) (Nat.mul_le_mul hGGm hGGm)
      have eN : Kb K0 (d + 1) ^ 2 * L = G * G * (G * G) * L := by
        show (G * G) ^ 2 * L = _; ring
      rw [eN] at h1 h2
      set N := G * G * (G * G) * L
      set N' := K ^ 2 * Lnext K G L
      set Gm4 := Gm * Gm * (Gm * Gm)
      have hs := N_step_23_25 K G L (by omega) hKG
      have hN'N : N' ≤ N + 2 * Gm4 := by omega
      have hN2 : N + 2 * Gm4 ≤ W := by
        have he : 2 * (d + 1) * Gm4 = 2 * d * Gm4 + 2 * Gm4 := by ring
        omega
      have h1' : N' + 2 * d * Gm4 ≤ W := by
        have he : 2 * (d + 1) * Gm4 = 2 * d * Gm4 + 2 * Gm4 := by ring
        omega
      have h23 : 23 ^ (d + 1) * N ≤ 25 ^ (d + 1) * W0 := by omega
      have hG0 : 0 < G := by omega
      have hl : L ≤ W := by
        have hlN : L ≤ N := Nat.le_mul_of_pos_left _
          (Nat.mul_pos (Nat.mul_pos hG0 hG0) (Nat.mul_pos hG0 hG0))
        omega
      refine max_le ?_ ?_
      · refine le_trans (eLev_le K K G L mb (N + 2 * Gm4)
          (by omega) le_rfl hKG (by omega) hN'N) ?_
        have hp : ePolyL K G L mb ≤ (Gm4 + 9) *
            (Nat.log 2 (W + mb + Gm + Gm + 64) + 1) :=
          ePolyL_mono (mb := mb) (le_trans hKG hGGm) hGGm hl
        have hb := hLev d N hN2 h23
        change 14 * (N + 2 * Gm4) / K + 7 +
          (Gm4 + 9) * (Nat.log 2 (W + mb + Gm + Gm + 64) + 1) ≤ Y at hb
        omega
      · refine chainMax_le_23_25 K0 Y mb W W0 Gm hK hLev hBr d _
          (le_trans (le_trans hKG (by
            show G ≤ G * G; exact Nat.le_mul_of_pos_left _ (by omega))) hGm) h1' ?_
        have a1 : 25 * (23 ^ d * N') ≤
            23 ^ (d + 1) * N + 25 ^ (d + 1) * (2 * Gm4) := by
          have b1 : 23 ^ d * (25 * N') ≤
              23 ^ d * (23 * N + 50 * (G * G * (G * G))) := Nat.mul_le_mul_left _ hs
          have b2 : 23 ^ d * (50 * (G * G * (G * G))) ≤ 25 ^ d * (50 * Gm4) :=
            Nat.mul_le_mul (Nat.pow_le_pow_left (by norm_num) d) (Nat.mul_le_mul_left _ hG4)
          have e1 : 25 * (23 ^ d * N') = 23 ^ d * (25 * N') := by ring
          have e2 : 23 ^ d * (23 * N + 50 * (G * G * (G * G))) =
              23 ^ (d + 1) * N + 23 ^ d * (50 * (G * G * (G * G))) := by ring
          have e3 : 25 ^ d * (50 * Gm4) = 25 ^ (d + 1) * (2 * Gm4) := by ring
          omega
        have e4 : 25 ^ (d + 1) * (2 * (d + 1) * Gm4) =
            25 ^ (d + 1) * (2 * Gm4) + 25 * (25 ^ d * (2 * d * Gm4)) := by ring
        have e5 : 25 ^ (d + 1) * W0 = 25 * (25 ^ d * W0) := by ring
        have hf : 25 * (23 ^ d * N' + 25 ^ d * (2 * d * Gm4)) ≤ 25 * (25 ^ d * W0) := by
          rw [mul_add]; omega
        exact Nat.le_of_mul_le_mul_left hf (by norm_num)

theorem group_count_dominates_25 (K0 d : ℕ) (hK : 600 ≤ K0) :
    14 * 25 ^ (d + 1) ≤ Kb K0 d * 23 ^ (d + 1) := by
  have hpow : (46 : ℕ) ^ d = 2 ^ d * 23 ^ d := by
    rw [show (46 : ℕ) = 2 * 23 by norm_num, mul_pow]
  have hlow : 600 * 2 ^ d ≤ Kb K0 d :=
    le_trans (Nat.mul_le_mul_right _ hK) (kb_ge_double K0 (by omega) d)
  calc 14 * 25 ^ (d + 1) = 350 * 25 ^ d := by ring
    _ ≤ 13800 * 46 ^ d :=
      Nat.mul_le_mul (by norm_num) (Nat.pow_le_pow_left (by norm_num : 25 ≤ 46) d)
    _ = (600 * 2 ^ d) * 23 ^ (d + 1) := by rw [hpow]; ring
    _ ≤ Kb K0 d * 23 ^ (d + 1) := Nat.mul_le_mul_right _ hlow

theorem level_volume_le_bottom_25 {K0 d N W0 : ℕ} (hK : 600 ≤ K0)
    (hN : 23 ^ (d + 1) * N ≤ 25 ^ (d + 1) * W0) :
    14 * N ≤ Kb K0 d * W0 := by
  have h := group_count_dominates_25 K0 d hK
  have hp : 0 < (23 : ℕ) ^ (d + 1) := by positivity
  have hscaled : 23 ^ (d + 1) * (14 * N) ≤ 23 ^ (d + 1) * (Kb K0 d * W0) := by
    calc 23 ^ (d + 1) * (14 * N) = 14 * (23 ^ (d + 1) * N) := by ring
      _ ≤ 14 * (25 ^ (d + 1) * W0) := Nat.mul_le_mul_left _ hN
      _ = (14 * 25 ^ (d + 1)) * W0 := by ring
      _ ≤ (Kb K0 d * 23 ^ (d + 1)) * W0 := Nat.mul_le_mul_right _ h
      _ = 23 ^ (d + 1) * (Kb K0 d * W0) := by ring
  exact Nat.le_of_mul_le_mul_left hscaled hp

theorem level_divided_volume_le_25 {K0 d N W0 G4 : ℕ} (hK : 600 ≤ K0)
    (hN : 23 ^ (d + 1) * N ≤ 25 ^ (d + 1) * W0) :
    14 * (N + 2 * G4) / Kb K0 d ≤ W0 + G4 := by
  have hmain := level_volume_le_bottom_25 hK hN
  have hbase : K0 ≤ Kb K0 d := Kb_mono K0 (by omega) (Nat.zero_le d)
  have hpad : 28 * G4 ≤ Kb K0 d * G4 := Nat.mul_le_mul_right _ (by omega)
  apply Nat.div_le_of_le_mul
  calc 14 * (N + 2 * G4) = 14 * N + 28 * G4 := by ring
    _ ≤ Kb K0 d * W0 + Kb K0 d * G4 := Nat.add_le_add hmain hpad
    _ = Kb K0 d * (W0 + G4) := by ring

def chainEnvelope (K0 mb W W0 G : ℕ) : ℕ :=
  W0 + G ^ 4 + 7 + (G ^ 4 + 9) * (Nat.log 2 (W + mb + G + G + 64) + 1) +
    2 * (Nat.log 2 (mb + K0 + 64) + 1)

theorem chainMax_le_envelope (K0 mb W W0 G : ℕ) (hK : 600 ≤ K0) :
    ∀ d L, Kb K0 d ≤ G →
      Kb K0 d ^ 2 * L + 2 * d * (G * G * (G * G)) ≤ W →
      23 ^ d * (Kb K0 d ^ 2 * L) + 25 ^ d * (2 * d * (G * G * (G * G))) ≤
        25 ^ d * W0 →
      chainMax mb K0 d L ≤ chainEnvelope K0 mb W W0 G := by
  apply chainMax_le_23_25 K0 (chainEnvelope K0 mb W W0 G) mb W W0 G hK
  · intro d N _ hN
    have hmain := level_divided_volume_le_25 (G4 := G * G * (G * G)) hK hN
    have hG : G * G * (G * G) = G ^ 4 := by ring
    rw [hG] at hmain ⊢
    unfold chainEnvelope
    omega
  · unfold chainEnvelope
    omega

/-- The explicit bottom budget, with sharper decay and unchanged padding. -/
def bottomBudget (K0 ell n : ℕ) : ℕ :=
  23 ^ ell * Ntn K0 ell n / 25 ^ ell + 1 + 2 * ell * Gm4 K0 ell

theorem chain_hyp2_25 (K0 ell n : ℕ) (hK : 600 ≤ K0) :
    23 ^ ell * (Kb K0 ell ^ 2 * Ltop (Kb K0 ell) (Gs K0 ell) (sOf n)) +
      25 ^ ell * (2 * ell * Gm4 K0 ell) ≤ 25 ^ ell * bottomBudget K0 ell n := by
  have h := Ntop_S (Kb K0 ell) (Gs K0 ell) (sOf n)
    (le_trans (by omega) (Kb_ge K0 (by omega) ell))
  have h1 : 23 ^ ell * (Kb K0 ell ^ 2 * Ltop (Kb K0 ell) (Gs K0 ell) (sOf n)) ≤
      23 ^ ell * Ntn K0 ell n := Nat.mul_le_mul_left _ h
  have h2 : 23 ^ ell * Ntn K0 ell n <
      25 ^ ell * (23 ^ ell * Ntn K0 ell n / 25 ^ ell + 1) := by
    have hd := Nat.div_add_mod (23 ^ ell * Ntn K0 ell n) (25 ^ ell)
    have hm := Nat.mod_lt (23 ^ ell * Ntn K0 ell n)
      (Nat.pow_pos (n := ell) (show 0 < 25 by norm_num))
    rw [mul_add]; omega
  unfold bottomBudget
  have he : 25 ^ ell * (23 ^ ell * Ntn K0 ell n / 25 ^ ell + 1 + 2 * ell * Gm4 K0 ell) =
      25 ^ ell * (23 ^ ell * Ntn K0 ell n / 25 ^ ell + 1) +
        25 ^ ell * (2 * ell * Gm4 K0 ell) := by ring
  omega

theorem chainMax_bound_on_input (K0 ell n : ℕ) (hK : 600 ≤ K0) :
    chainMax (mbn K0 ell n) K0 ell (Ltop (Kb K0 ell) (Gs K0 ell) (sOf n)) ≤
      chainEnvelope K0 (mbn K0 ell n) (Wn K0 ell n) (bottomBudget K0 ell n) (Gs K0 ell) := by
  exact chainMax_le_envelope _ _ _ _ _ hK ell _
    (Kb_le_Gs ell (by omega)) (chain_hyp1 ell (by omega) n) (chain_hyp2_25 K0 ell n hK)

theorem zC_le_envelope (hK : 600 ≤ K0u) (ell n : ℕ) :
    zC ell n ≤ chainEnvelope K0u (mbn K0u ell n) (Wn K0u ell n)
      (bottomBudget K0u ell n) (Gs K0u ell) := chainMax_bound_on_input _ _ _ hK

theorem weighted_bottomBudget_le (K0 ell n R C : ℕ)
    (hR : R * 23 ^ ell ≤ C * 25 ^ ell) :
    R * bottomBudget K0 ell n ≤ C * Ntn K0 ell n + R + 2 * R * ell * Gm4 K0 ell := by
  have hmain : R * (23 ^ ell * Ntn K0 ell n / 25 ^ ell) ≤ C * Ntn K0 ell n := by
    refine le_trans (Nat.mul_div_le_mul_div_assoc _ _ _) (Nat.div_le_of_le_mul ?_)
    calc R * (23 ^ ell * Ntn K0 ell n) = (R * 23 ^ ell) * Ntn K0 ell n := by ring
      _ ≤ (C * 25 ^ ell) * Ntn K0 ell n := Nat.mul_le_mul_right _ hR
      _ = 25 ^ ell * (C * Ntn K0 ell n) := by ring
  unfold bottomBudget
  nlinarith

end Factoring.K6

#print axioms Factoring.K6.N_step_23_25
#print axioms Factoring.K6.chainMax_le_23_25
#print axioms Factoring.K6.chainMax_le_envelope
#print axioms Factoring.K6.chainMax_bound_on_input
#print axioms Factoring.K6.zC_le_envelope
#print axioms Factoring.K6.weighted_bottomBudget_le
