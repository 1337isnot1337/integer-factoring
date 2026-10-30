/-
SPDX-License-Identifier: Apache-2.0

Quantitative bounds on the top construction and bookkeeping exponents.
Uses the pinned reference's proved cost inequalities. These statements are
independent of the stronger chain contraction and the degree-twelve schedule.
-/
module

public import FactoringAux.UniExp

@[expose] public section

namespace IntegerFactoring.K6Budget

open Factoring.Aux

def topCoefficient : ℕ := 17280 + 8 + 64 * 1662 ^ 4 + 1602
def volumeCoefficient : ℕ := 5760 + 8 + 64 * 306 ^ 4

theorem parameter_power_le_input {ell n a : ℕ} (h : UGood ell n) (ha : a ≤ 40) :
    Kb K0u ell ^ a ≤ n := by
  exact le_trans (Nat.pow_le_pow_right (by have := Kb48_ge48 ell; omega) ha) h.2.1

theorem parameter_le_input {ell n : ℕ} (h : UGood ell n) : Kb K0u ell ≤ n := by
  simpa using parameter_power_le_input h (by norm_num : 1 ≤ 40)

theorem grid_fourth_le_input {ell n : ℕ} (h : UGood ell n) :
    (2 ^ sOf n) ^ 4 ≤ 8 * n := by
  have hS := S_le_n n
  have hn := (good_facts h).2.2.1
  omega

/-- A convenient polynomial bound with spare room for either schedule. -/
theorem zP_le_poly {ell n : ℕ} (h : UGood ell n) :
    zP ell n ≤ 40000 * Kb K0u ell ^ 38 := by
  have hP := Kb48_ge48 ell
  obtain ⟨hG, _, h37, _⟩ := good_bounds h
  obtain ⟨a1, a2, a3, _, _, _, a7⟩ := good_logs h
  have a4 : Nat.log 2 (Kb K0u ell) + 1 ≤ 801 * Kb K0u ell :=
    le_trans (Nat.add_le_add_right (Nat.log_mono_right
      (by omega : Kb K0u ell ≤ Kb K0u ell + 2)) 1) a7
  set P := Kb K0u ell
  have hpw : ∀ i j, i ≤ j → P ^ i ≤ P ^ j :=
    fun i j hij => Nat.pow_le_pow_right (by omega) hij
  have p1 : P ≤ P ^ 38 := by simpa using hpw 1 38 (by norm_num)
  have p9 : 50 ≤ P ^ 9 := by
    have h2 : 48 ^ 2 ≤ P ^ 2 := Nat.pow_le_pow_left hP 2
    exact le_trans (by norm_num : 50 ≤ 48 ^ 2) (le_trans h2 (hpw 2 9 (by norm_num)))
  have t1 : ell * Gm4 K0u ell * (Nat.log 2 (Wn K0u ell n + 1) + 1)
      ≤ 801 * P ^ 38 := by
    calc _ ≤ P ^ 37 * (801 * P) := Nat.mul_le_mul h37 a2
      _ = _ := by ring
  have t5 : (Gs K0u ell + 50) * (Nat.log 2 (P + 2) + 1)
      ≤ 1602 * P ^ 38 := by
    calc _ ≤ (2 * P ^ 9) * (801 * P) := Nat.mul_le_mul (by omega) a7
      _ = 1602 * P ^ 10 := by ring
      _ ≤ _ := Nat.mul_le_mul_left _ (hpw 10 38 (by norm_num))
  change ell * Gm4 K0u ell * (Nat.log 2 (Wn K0u ell n + 1) + 1) +
    (Nat.log 2 (mbn K0u ell n) + 1) + (Nat.log 2 (n + ell + 64) + 1) +
    40 * (Nat.log 2 P + 1) + (Gs K0u ell + 50) * (Nat.log 2 (P + 2) + 1)
      ≤ 40000 * P ^ 38
  nlinarith only [t1, t5, a1, a3, a4, p1]

theorem weighted_zP_le {ell n R : ℕ} (h : UGood ell n) (hR : R ≤ Kb K0u ell) :
    R * zP ell n ≤ 40000 * n := by
  calc R * zP ell n ≤ Kb K0u ell * (40000 * Kb K0u ell ^ 38) :=
      Nat.mul_le_mul hR (zP_le_poly h)
    _ = 40000 * Kb K0u ell ^ 39 := by ring
    _ ≤ _ := Nat.mul_le_mul_left _ (parameter_power_le_input h (by norm_num))

/-- Uniform weighted bounds for the three terms in the top cost. -/
theorem weighted_zT_le {ell n R : ℕ} (h : UGood ell n) (hR : R ≤ Kb K0u ell) :
    R * zT ell n ≤ topCoefficient * n := by
  set P := Kb K0u ell
  set G := Gs K0u ell
  set S := 2 ^ sOf n
  have hP := Kb48_ge48 ell
  have hG : G ≤ P ^ 9 := (good_bounds h).1
  have hS : S ^ 4 ≤ 8 * n := grid_fourth_le_input h
  have hF1 : R * (2160 * S ^ 4 / P) ≤ 17280 * n := by
    calc _ ≤ P * (2160 * S ^ 4 / P) := Nat.mul_le_mul_right _ hR
      _ ≤ 2160 * S ^ 4 := Nat.mul_div_le _ _
      _ ≤ 2160 * (8 * n) := Nat.mul_le_mul_left _ hS
      _ = _ := by ring
  have hRG : R * G ≤ P ^ 10 := by
    calc _ ≤ P * P ^ 9 := Nat.mul_le_mul hR hG
      _ = _ := by ring
  have hbeta : 1662 * R * G ≤ 1662 * P ^ 10 := by
    simpa [mul_assoc] using Nat.mul_le_mul_left 1662 hRG
  have hbeta4 : (1662 * R * G) ^ 4 ≤ 1662 ^ 4 * n := by
    calc _ ≤ (1662 * P ^ 10) ^ 4 := Nat.pow_le_pow_left hbeta 4
      _ = 1662 ^ 4 * P ^ 40 := by ring
      _ ≤ _ := Nat.mul_le_mul_left _ h.2.1
  have hF2 : R * (1662 * (G * S ^ 3)) ≤ (8 + 64 * 1662 ^ 4) * n := by
    have hab := absorb_cube (1662 * R * G) S
    have hd := Nat.div_le_self (S ^ 4) 4
    have he : (1662 * R * G) * (4 * (1662 * R * G)) ^ 3 =
        64 * (1662 * R * G) ^ 4 := by ring
    have he' : R * (1662 * (G * S ^ 3)) = (1662 * R * G) * S ^ 3 := by ring
    rw [he] at hab
    rw [he']
    nlinarith only [hab, hd, hS, hbeta4]
  have hF3 : R * ePolyT P G (sOf n) n (mbn K0u ell n) ≤ 1602 * n := by
    have hlog := (good_logs h).2.2.2.1
    have hp9 : 4 ≤ P ^ 9 := le_trans (by omega : 4 ≤ P)
      (Nat.le_self_pow (by norm_num) _)
    have hG4 : G + 4 ≤ 2 * P ^ 9 := by omega
    unfold ePolyT
    calc _ ≤ P * ((2 * P ^ 9) * (801 * P)) :=
        Nat.mul_le_mul hR (Nat.mul_le_mul hG4 hlog)
      _ = 1602 * P ^ 11 := by ring
      _ ≤ _ := Nat.mul_le_mul_left _ (parameter_power_le_input h (by norm_num))
  have ht := Nat.mul_le_mul_left R (eTop_split ell n (mbn K0u ell n))
  change R * eTop K0u ell n (mbn K0u ell n) ≤ topCoefficient * n
  unfold topCoefficient
  dsimp only [P, G, S] at hF1 hF2 hF3
  nlinarith only [ht, hF1, hF2, hF3]

/-- The top network volume is bounded linearly for every good input pair. -/
theorem top_volume_le {ell n : ℕ} (h : UGood ell n) :
    Ntn K0u ell n ≤ volumeCoefficient * n := by
  set P := Kb K0u ell
  set G := Gs K0u ell
  set S := 2 ^ sOf n
  have hG : G ≤ P ^ 9 := (good_bounds h).1
  have hS : S ^ 4 ≤ 8 * n := grid_fourth_le_input h
  have hbeta : (306 * G) ^ 4 ≤ 306 ^ 4 * n := by
    calc _ ≤ (306 * P ^ 9) ^ 4 :=
        Nat.pow_le_pow_left (Nat.mul_le_mul_left _ hG) 4
      _ = 306 ^ 4 * P ^ 36 := by ring
      _ ≤ _ := Nat.mul_le_mul_left _ (parameter_power_le_input h (by norm_num))
  have hab := absorb_cube (306 * G) S
  have hd := Nat.div_le_self (S ^ 4) 4
  have he : (306 * G) * (4 * (306 * G)) ^ 3 = 64 * (306 * G) ^ 4 := by ring
  rw [he] at hab
  change 720 * S ^ 4 + 306 * (G * S ^ 3) ≤ volumeCoefficient * n
  unfold volumeCoefficient
  nlinarith only [hab, hd, hS, hbeta]

/-- All polynomial parts of the chain envelope fit after multiplication by R. -/
theorem weighted_chain_padding_le {ell n R : ℕ} (h : UGood ell n)
    (hR : R ≤ Kb K0u ell) :
    R * (Gm4 K0u ell + 7 + (Gm4 K0u ell + 9) *
      (Nat.log 2 (Wn K0u ell n + mbn K0u ell n + Gs K0u ell + Gs K0u ell + 64) + 1) +
      2 * (Nat.log 2 (mbn K0u ell n + K0u + 64) + 1)) ≤ 9620 * n := by
  set P := Kb K0u ell
  have hP := Kb48_ge48 ell
  have hG4 : Gm4 K0u ell ≤ P ^ 36 := (good_bounds h).2.1
  obtain ⟨_, _, _, _, hlog, hbase, _⟩ := good_logs h
  have hp36 : 1 ≤ P ^ 36 := Nat.one_le_pow _ _ (by omega)
  have a : R * Gm4 K0u ell ≤ n := by
    calc _ ≤ P * P ^ 36 := Nat.mul_le_mul hR hG4
      _ = P ^ 37 := by ring
      _ ≤ _ := parameter_power_le_input h (by norm_num)
  have b : R ≤ n := le_trans hR (parameter_le_input h)
  have c : R * ((Gm4 K0u ell + 9) *
      (Nat.log 2 (Wn K0u ell n + mbn K0u ell n + Gs K0u ell + Gs K0u ell + 64) + 1))
      ≤ 8010 * n := by
    calc _ ≤ P * ((10 * P ^ 36) * (801 * P)) :=
        Nat.mul_le_mul hR (Nat.mul_le_mul (by omega) hlog)
      _ = 8010 * P ^ 38 := by ring
      _ ≤ _ := Nat.mul_le_mul_left _ (parameter_power_le_input h (by norm_num))
  have d : R * (2 * (Nat.log 2 (mbn K0u ell n + K0u + 64) + 1)) ≤ 1602 * n := by
    calc _ ≤ P * (2 * (801 * P)) :=
        Nat.mul_le_mul hR (Nat.mul_le_mul_left _ hbase)
      _ = 1602 * P ^ 2 := by ring
      _ ≤ _ := Nat.mul_le_mul_left _ (parameter_power_le_input h (by norm_num))
  nlinarith only [a, b, c, d]

theorem weighted_depth_padding_le {ell n R : ℕ} (h : UGood ell n)
    (hR : R ≤ Kb K0u ell) : R * (1 + 2 * ell * Gm4 K0u ell) ≤ 3 * n := by
  have he := (good_bounds h).2.2.1
  have hr : R ≤ n := le_trans hR (parameter_le_input h)
  have hm : R * (ell * Gm4 K0u ell) ≤ n := by
    calc _ ≤ Kb K0u ell * Kb K0u ell ^ 37 := Nat.mul_le_mul hR he
      _ = Kb K0u ell ^ 38 := by ring
      _ ≤ _ := parameter_power_le_input h (by norm_num)
  nlinarith only [hr, hm]

end IntegerFactoring.K6Budget

#print axioms IntegerFactoring.K6Budget.weighted_zT_le
#print axioms IntegerFactoring.K6Budget.weighted_zP_le
#print axioms IntegerFactoring.K6Budget.top_volume_le
#print axioms IntegerFactoring.K6Budget.weighted_chain_padding_le
