/-
SPDX-License-Identifier: Apache-2.0

An explicit upper bound on every exponent along the upstream counting chain.
Uses the development in IsaacRudich/SubexponentialDeterministicFactoring,
commit c91cc3f4947caca8977e81a5d14244f73b5ae2b1, under Apache-2.0.

Unlike the source's small_zC proof, this estimate uses no split at a depth
depending on the requested exponential rate. It keeps the bottom-volume
bound W0 and polynomial terms visible. The final theorem carries this
explicit envelope through the unchanged counting program, compiler, and
finite-alphabet machine simulator. The conversion of that envelope to a
logarithmic-power asymptotic bound is a separate obligation.
-/
module

public import FactoringAux.UniMain

@[expose] public section

namespace Factoring.Improvements

open Factoring.Aux

/-- The growth in group counts dominates the backwards volume ratio. -/
theorem group_count_dominates_backwards_ratio (d : ℕ) :
    14 * 24 ^ (d + 1) ≤ Kb K0u d * 23 ^ (d + 1) := by
  have hpow : (46 : ℕ) ^ d = 2 ^ d * 23 ^ d := by
    rw [show (46 : ℕ) = 2 * 23 by norm_num, mul_pow]
  calc 14 * 24 ^ (d + 1) = 336 * 24 ^ d := by ring
    _ ≤ 1104 * 46 ^ d :=
      Nat.mul_le_mul (by norm_num) (Nat.pow_le_pow_left (by norm_num : 24 ≤ 46) d)
    _ = (48 * 2 ^ d) * 23 ^ (d + 1) := by rw [hpow]; ring
    _ ≤ Kb K0u d * 23 ^ (d + 1) := Nat.mul_le_mul_right _ (Kb48_ge d)

/-- Every level's main volume cost is bounded by the bottom-volume budget. -/
theorem level_volume_le_bottom {d N W0 : ℕ}
    (hN : 23 ^ (d + 1) * N ≤ 24 ^ (d + 1) * W0) :
    14 * N ≤ Kb K0u d * W0 := by
  have h := group_count_dominates_backwards_ratio d
  have hp : 0 < (23 : ℕ) ^ (d + 1) := by positivity
  have hscaled : 23 ^ (d + 1) * (14 * N) ≤ 23 ^ (d + 1) * (Kb K0u d * W0) := by
    calc 23 ^ (d + 1) * (14 * N) = 14 * (23 ^ (d + 1) * N) := by ring
      _ ≤ 14 * (24 ^ (d + 1) * W0) := Nat.mul_le_mul_left _ hN
      _ = (14 * 24 ^ (d + 1)) * W0 := by ring
      _ ≤ (Kb K0u d * 23 ^ (d + 1)) * W0 := Nat.mul_le_mul_right _ h
      _ = 23 ^ (d + 1) * (Kb K0u d * W0) := by ring
  exact Nat.le_of_mul_le_mul_left hscaled hp

/-- The additive padding term also fits a fixed explicit budget. -/
theorem level_divided_volume_le {d N W0 G4 : ℕ}
    (hN : 23 ^ (d + 1) * N ≤ 24 ^ (d + 1) * W0) :
    14 * (N + 2 * G4) / Kb K0u d ≤ W0 + G4 := by
  have hmain := level_volume_le_bottom hN
  have hpad : 28 * G4 ≤ Kb K0u d * G4 :=
    Nat.mul_le_mul_right _ (le_trans (by norm_num : 28 ≤ 48) (Kb48_ge48 d))
  apply Nat.div_le_of_le_mul
  calc 14 * (N + 2 * G4) = 14 * N + 28 * G4 := by ring
    _ ≤ Kb K0u d * W0 + Kb K0u d * G4 := Nat.add_le_add hmain hpad
    _ = Kb K0u d * (W0 + G4) := by ring

/-- An explicit bound on the maximum exponent in the chain. -/
def chainEnvelope (mb W W0 G : ℕ) : ℕ :=
  W0 + G ^ 4 + 7 + (G ^ 4 + 9) * (Nat.log 2 (W + mb + G + G + 64) + 1) +
    2 * (Nat.log 2 (mb + K0u + 64) + 1)

/-- Every level and the final brute-force evaluation fit the same explicit
envelope; the depth does not appear in an additional cutoff. -/
theorem chainMax_le_explicit (mb W W0 G : ℕ) :
    ∀ d L, Kb K0u d ≤ G →
      Kb K0u d ^ 2 * L + 2 * d * (G * G * (G * G)) ≤ W →
      23 ^ d * (Kb K0u d ^ 2 * L) + 24 ^ d * (2 * d * (G * G * (G * G))) ≤
        24 ^ d * W0 →
      chainMax mb K0u d L ≤ chainEnvelope mb W W0 G := by
  apply chainMax_le K0u (chainEnvelope mb W W0 G) mb W W0 G (by norm_num)
  · intro d N _ hN
    have hmain := level_divided_volume_le (G4 := G * G * (G * G)) hN
    have hG : G * G * (G * G) = G ^ 4 := by ring
    rw [hG] at hmain ⊢
    unfold chainEnvelope
    omega
  · unfold chainEnvelope
    omega

/-- The explicit estimate applies to the source's actual chain exponent at
every selected pair of input length and depth. -/
theorem zC_le_explicit (ell n : ℕ) :
    zC ell n ≤ chainEnvelope (mbn K0u ell n) (Wn K0u ell n)
      (W0n K0u ell n) (Gs K0u ell) := by
  exact chainMax_le_explicit _ _ _ _ ell _
    (Kb_le_Gs ell (by norm_num)) (chain_hyp1 ell (by norm_num) n)
    (chain_hyp2 ell (by norm_num) n)

/-- In particular the bound holds at the runtime-selected depth. -/
theorem uniform_chain_explicit (n : ℕ) :
    zC (ellOf n) n ≤ chainEnvelope (mbn K0u (ellOf n) n) (Wn K0u (ellOf n) n)
      (W0n K0u (ellOf n) n) (Gs K0u (ellOf n)) := zC_le_explicit _ _

/-- Preserve the top and bookkeeping exponents while replacing the hidden
maximum over levels by its explicit envelope. -/
def explicitZu (n : ℕ) : ℕ :=
  zT (ellOf n) n +
    chainEnvelope (mbn K0u (ellOf n) n) (Wn K0u (ellOf n) n)
      (W0n K0u (ellOf n) n) (Gs K0u (ellOf n)) +
    zP (ellOf n) n

theorem Zu_le_explicit (n : ℕ) : Zu n ≤ explicitZu n := by
  have h := uniform_chain_explicit n
  unfold Zu Zof explicitZu
  omega

/-! ## The unchanged factoring program has a machine with the explicit time bound. -/

/-- Reconstruct the source's value requirements with a named exponent instead
of using the existential witness in `uniReq_exp`. -/
theorem uniReq_of_Zu_le (n Z : ℕ) (hZ : Zu n ≤ Z) :
    UniReq (2 ^ (16 * Z + 127)) n := by
  have hT : zT (ellOf n) n ≤ Z := by unfold Zu Zof at hZ; omega
  have hC : zC (ellOf n) n ≤ Z := by unfold Zu Zof at hZ; omega
  have hP : zP (ellOf n) n ≤ Z := by unfold Zu Zof at hZ; omega
  have hmain : MainReq (2 ^ Z) K0u (ellOf n) n :=
    mainReq_of (by norm_num) hT (deepE_chainMax _ _ _ _ _ hC) le_rfl
  have hNP : NPof K0u (ellOf n)
      (Ltop (Kb K0u (ellOf n)) (Gs K0u (ellOf n)) (sOf n)) ≤ 2 ^ Z := by
    unfold zP at hP
    refine le_trans (NPof_le K0u (Wn K0u (ellOf n) n) (Gs K0u (ellOf n))
      (le_refl _) (ellOf n) _ (Kb_le_Gs (K0 := K0u) (ellOf n) (by norm_num))
      (chain_hyp1 (K0 := K0u) (ellOf n) (by norm_num) n)) ?_
    refine le_trans (pow_le_two_pow_log _ _) (two_pow_le ?_)
    show ellOf n * Gm4 K0u (ellOf n) *
      (Nat.log 2 (Wn K0u (ellOf n) n + 1) + 1) ≤ Z
    omega
  have hS : n + ellOf n + 64 ≤ 2 ^ Z := by
    unfold zP at hP
    exact le_trans (le_two_pow_log _).le (two_pow_le (by omega))
  have hK : Kb K0u (ellOf n) ^ 40 ≤ 2 ^ Z := by
    unfold zP at hP
    refine le_trans (pow_le_two_pow_log _ _) (two_pow_le ?_)
    omega
  have hL : lvlV (Kb K0u (ellOf n)) ≤ 2 ^ Z := by
    unfold zP at hP
    refine le_trans (pow_le_two_pow_log _ _) (two_pow_le ?_)
    omega
  have hmono : 2 ^ Z ≤ 2 ^ (16 * Z + 127) := two_pow_le (by omega)
  refine ⟨hmain.mono hmono, le_trans hNP hmono, le_trans hL hmono, le_trans hS hmono, ?_⟩
  rw [Kb_succ_pow]
  calc 9 ^ 40 * (Kb K0u (ellOf n) ^ 40) ^ 16 ≤ 2 ^ 127 * (2 ^ Z) ^ 16 :=
      Nat.mul_le_mul (by norm_num) (Nat.pow_le_pow_left hK 16)
    _ = 2 ^ (16 * Z + 127) := by rw [← pow_mul, ← pow_add]; ring_nf

def explicitValueBound (n : ℕ) : ℕ := 2 ^ (16 * explicitZu n + 127)

theorem explicit_value_requirements (n : ℕ) : UniReq (explicitValueBound n) n :=
  uniReq_of_Zu_le n (explicitZu n) (Zu_le_explicit n)

/-- A fixed polynomial counting cost stays a fixed polynomial after the
source's quadratic compilation and cubic machine simulation. -/
theorem simulation_cost_polynomial (n B V c E Q k0 cc : ℕ)
    (hn : n + 1 ≤ V) (hB : B ≤ c * V ^ E) :
    cc * (Q * (n + k0 + 200 * (n + 1) ^ 3 * (B + 1) + 1) ^ 2 + n + 1) ^ 3 ≤
      cc * (Q * (k0 + 200 * (c + 1) + 1) ^ 2 + 1) ^ 3 * V ^ (6 * (E + 3)) := by
  have hV : 1 ≤ V := by omega
  have hVE : 1 ≤ V ^ E := Nat.one_le_pow _ _ hV
  have hB1 : B + 1 ≤ (c + 1) * V ^ E := by nlinarith
  have hN3 : (n + 1) ^ 3 ≤ V ^ 3 := Nat.pow_le_pow_left hn 3
  have hS : 200 * (n + 1) ^ 3 * (B + 1) ≤ 200 * (c + 1) * V ^ (E + 3) := by
    calc 200 * (n + 1) ^ 3 * (B + 1) ≤ 200 * V ^ 3 * ((c + 1) * V ^ E) := by gcongr
      _ = 200 * (c + 1) * V ^ (E + 3) := by rw [pow_add]; ring
  have hVF : V ≤ V ^ (E + 3) := by
    calc V = V ^ 1 := by simp
      _ ≤ V ^ (E + 3) := Nat.pow_le_pow_right hV (by omega)
  have hF : 1 ≤ V ^ (E + 3) := Nat.one_le_pow _ _ hV
  have hN : n + 1 ≤ V ^ (E + 3) := le_trans hn hVF
  have hk0 : k0 ≤ k0 * V ^ (E + 3) := Nat.le_mul_of_pos_right _ hF
  have hbase : n + k0 + 200 * (n + 1) ^ 3 * (B + 1) + 1 ≤
      (k0 + 200 * (c + 1) + 1) * V ^ (E + 3) := by nlinarith
  have hsq : (n + k0 + 200 * (n + 1) ^ 3 * (B + 1) + 1) ^ 2 ≤
      (k0 + 200 * (c + 1) + 1) ^ 2 * V ^ (2 * (E + 3)) := by
    calc _ ≤ ((k0 + 200 * (c + 1) + 1) * V ^ (E + 3)) ^ 2 := Nat.pow_le_pow_left hbase 2
      _ = _ := by rw [mul_pow, ← pow_mul]; congr 2; omega
  have hN2 : n + 1 ≤ V ^ (2 * (E + 3)) := le_trans hN
    (Nat.pow_le_pow_right hV (by omega))
  have houter : Q * (n + k0 + 200 * (n + 1) ^ 3 * (B + 1) + 1) ^ 2 + n + 1 ≤
      (Q * (k0 + 200 * (c + 1) + 1) ^ 2 + 1) * V ^ (2 * (E + 3)) := by
    have := Nat.mul_le_mul_left Q hsq
    nlinarith
  calc cc * (Q * (n + k0 + 200 * (n + 1) ^ 3 * (B + 1) + 1) ^ 2 + n + 1) ^ 3 ≤
      cc * ((Q * (k0 + 200 * (c + 1) + 1) ^ 2 + 1) * V ^ (2 * (E + 3))) ^ 3 :=
        Nat.mul_le_mul_left _ (Nat.pow_le_pow_left houter 3)
    _ = _ := by rw [mul_pow, ← pow_mul, ← mul_assoc]; congr 2; ring

open Factoring.Aux.Copied Factoring.Aux.Copied.PLib Computability

/-- One finite-alphabet machine factors in a fixed polynomial of the explicit
value bound. The compiler and simulator are the source's checked ones. -/
theorem factoring_explicit_polynomial :
    ∃ C E : ℕ, FactoringInTime (fun n => C * (explicitValueBound n + 1) ^ E) := by
  obtain ⟨c, E, hcost⟩ := uniC_cost
  let B : ℕ → ℕ := fun n => c * (explicitValueBound n + 1) ^ E
  have hc : CSpec (uniC UP) B := by
    refine ⟨fun st => ?_, uniC_le_one, fun st hkn hN hD _ => uniC_val st hkn hN hD,
      fun st hkn hA => hcost _ st hkn hA (explicit_value_requirements _)⟩
    obtain ⟨h1, h2⟩ := uniC_keeps st
    exact ⟨h1 0 (by norm_num) (by simp [Uni.Wu]), h1 2 (by norm_num) (by simp [Uni.Wu]),
      h1 4 (by norm_num) (by simp [Uni.Wu]), fun a ha => h2 a (by omega)⟩
  obtain ⟨hout, hdrv⟩ := drvC_spec (uniC UP) B hc
  obtain ⟨Q, k0, hT⟩ := timeComp_wrapFn (drvC (uniC UP)) 2
    (fun xs => outOf ((srcOut (drvC (uniC UP)) xs).ar 2))
    (fun xs => 200 * (xs.length + 1) ^ 3 * (B xs.length + 1)) (fun _ => rfl) hdrv
  obtain ⟨cc, M, hfin, htime⟩ := Sim.factoringInTime_of_timeComp (wrapFn (drvC (uniC UP)) 2) _ _ hT hout
    (fun n => Q * (n + k0 + 200 * (n + 1) ^ 3 * (B n + 1) + 1) ^ 2) (fun N => le_rfl)
  refine ⟨cc * (Q * (k0 + 200 * (c + 1) + 1) ^ 2 + 1) ^ 3, 6 * (E + 3), M, hfin, fun n => ?_⟩
  have hn : n + 1 ≤ explicitValueBound n + 1 := by
    have := (explicit_value_requirements n).small
    omega
  exact le_trans (htime n)
    (simulation_cost_polynomial n (B n) (explicitValueBound n + 1) c E Q k0 cc hn le_rfl)

/-- An end-to-end bound obtained from the unchanged factoring program, with
every exponent quantity explicitly defined above. No logarithmic-power
asymptotic conversion is asserted in this theorem. -/
theorem factoring_explicit_exponential :
    ∃ M : FactoringMachine, ∃ C : ℕ,
      FactorsWithin M (fun n => C * 2 ^ (C * (explicitZu n + 1))) := by
  obtain ⟨c, E, M, hfin, htime⟩ := factoring_explicit_polynomial
  refine ⟨M, c + 128 * E + 1, hfin, fun n => ?_⟩
  have hpow : explicitValueBound n + 1 ≤ 2 ^ (16 * explicitZu n + 128) := by
    unfold explicitValueBound
    have hp := Nat.one_le_two_pow (n := 16 * explicitZu n + 127)
    calc 2 ^ (16 * explicitZu n + 127) + 1 ≤ 2 ^ (16 * explicitZu n + 127) * 2 := by omega
      _ = 2 ^ (16 * explicitZu n + 128) := (pow_succ (2 : ℕ) (16 * explicitZu n + 127)).symm
  calc M.time n ≤ c * (explicitValueBound n + 1) ^ E := htime n
    _ ≤ c * (2 ^ (16 * explicitZu n + 128)) ^ E :=
      Nat.mul_le_mul_left _ (Nat.pow_le_pow_left hpow E)
    _ = c * 2 ^ (E * (16 * explicitZu n + 128)) := by rw [← pow_mul]; congr 2; ring
    _ ≤ (c + 128 * E + 1) * 2 ^ ((c + 128 * E + 1) * (explicitZu n + 1)) :=
      Nat.mul_le_mul (by omega) (two_pow_le (by nlinarith))

end Factoring.Improvements

#print axioms Factoring.Improvements.group_count_dominates_backwards_ratio
#print axioms Factoring.Improvements.level_divided_volume_le
#print axioms Factoring.Improvements.chainMax_le_explicit
#print axioms Factoring.Improvements.zC_le_explicit
#print axioms Factoring.Improvements.uniform_chain_explicit
#print axioms Factoring.Improvements.Zu_le_explicit
#print axioms Factoring.Improvements.explicit_value_requirements
#print axioms Factoring.Improvements.factoring_explicit_polynomial
#print axioms Factoring.Improvements.factoring_explicit_exponential
