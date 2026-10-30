/-
SPDX-License-Identifier: Apache-2.0

A quantitative bridge from the uniform engine exponent to one factoring
machine. The growth assumption is explicit, so this module can be checked
against both the pinned reference interfaces and the modified k6 engine.

The counting program, compiler, and finite-alphabet simulator are the checked
ones from IsaacRudich/SubexponentialDeterministicFactoring, commit
c91cc3f4947caca8977e81a5d14244f73b5ae2b1, under Apache-2.0.
-/
module

public import FactoringAux.UniMain

@[expose] public section

namespace IntegerFactoring.K6Machine

open Factoring Factoring.Aux

/-- Reconstruct the source's value requirements with a named exponent instead
of using the existential witness in `uniReq_exp`. -/
theorem uniReq_of_Zu_le
    (hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12) (n Z : ℕ) (hZ : Zu n ≤ Z) :
    UniReq (2 ^ (12 * Z + 186)) n := by
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
  have hmono : 2 ^ Z ≤ 2 ^ (12 * Z + 186) := two_pow_le (by omega)
  refine ⟨hmain.mono hmono, le_trans hNP hmono, le_trans hL hmono, le_trans hS hmono, ?_⟩
  calc Kb K0u (ellOf n + 1) ^ 40 ≤ (25 * Kb K0u (ellOf n) ^ 12) ^ 40 :=
      Nat.pow_le_pow_left (hstep _) 40
    _ = 25 ^ 40 * (Kb K0u (ellOf n) ^ 40) ^ 12 := by ring
    _ ≤ 2 ^ 186 * (2 ^ Z) ^ 12 :=
      Nat.mul_le_mul (by norm_num) (Nat.pow_le_pow_left hK 12)
    _ = 2 ^ (12 * Z + 186) := by rw [← pow_mul, ← pow_add]; ring_nf

/-- The uniform value bound for a schedule with growth degree twelve. -/
def valueBound (n : ℕ) : ℕ := 2 ^ (12 * Zu n + 186)

theorem value_requirements
    (hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12) (n : ℕ) :
    UniReq (valueBound n) n := uniReq_of_Zu_le hstep n (Zu n) le_rfl

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
theorem factoring_Zu_polynomial
    (hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12) :
    ∃ C E : ℕ, FactoringInTime (fun n => C * (valueBound n + 1) ^ E) := by
  obtain ⟨c, E, hcost⟩ := uniC_cost
  let B : ℕ → ℕ := fun n => c * (valueBound n + 1) ^ E
  have hc : CSpec (uniC UP) B := by
    refine ⟨fun st => ?_, uniC_le_one, fun st hkn hN hD _ => uniC_val st hkn hN hD,
      fun st hkn hA => hcost _ st hkn hA (value_requirements hstep _)⟩
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
  have hn : n + 1 ≤ valueBound n + 1 := by
    have := (value_requirements hstep n).small
    omega
  exact le_trans (htime n)
    (simulation_cost_polynomial n (B n) (valueBound n + 1) c E Q k0 cc hn le_rfl)

/-- A single finite-alphabet machine factors within an exponential in the
actual uniform cost exponent, including compilation and simulation. -/
theorem factoring_Zu_exponential
    (hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12) :
    ∃ M : FactoringMachine, ∃ C : ℕ,
      FactorsWithin M (fun n => C * 2 ^ (C * (Zu n + 1))) := by
  obtain ⟨c, E, M, hfin, htime⟩ := factoring_Zu_polynomial hstep
  refine ⟨M, c + 187 * E + 1, hfin, fun n => ?_⟩
  have hpow : valueBound n + 1 ≤ 2 ^ (12 * Zu n + 187) := by
    unfold valueBound
    have hp := Nat.one_le_two_pow (n := 12 * Zu n + 186)
    calc 2 ^ (12 * Zu n + 186) + 1 ≤ 2 ^ (12 * Zu n + 186) * 2 := by omega
      _ = 2 ^ (12 * Zu n + 187) := (pow_succ (2 : ℕ) (12 * Zu n + 186)).symm
  calc M.time n ≤ c * (valueBound n + 1) ^ E := htime n
    _ ≤ c * (2 ^ (12 * Zu n + 187)) ^ E :=
      Nat.mul_le_mul_left _ (Nat.pow_le_pow_left hpow E)
    _ = c * 2 ^ (E * (12 * Zu n + 187)) := by rw [← pow_mul]; congr 2; ring
    _ ≤ (c + 187 * E + 1) * 2 ^ ((c + 187 * E + 1) * (Zu n + 1)) :=
      Nat.mul_le_mul (by omega) (two_pow_le (by nlinarith))

/-- Any uniform bound on the engine exponent transfers to the same form of
factoring time. The denominator may be any natural-valued function. -/
theorem factoring_rate_of_Zu_bound (D : ℕ → ℕ)
    (hstep : ∀ d, Kb K0u (d + 1) ≤ 25 * Kb K0u d ^ 12)
    (hZ : ∃ A : ℕ, ∀ n, Zu n ≤ A * (n / D n + 1)) :
    ∃ M : FactoringMachine, ∃ C : ℕ,
      FactorsWithin M (fun n => C * 2 ^ (C * (n / D n + 1))) := by
  obtain ⟨M, c, hfin, htime⟩ := factoring_Zu_exponential hstep
  obtain ⟨A, hA⟩ := hZ
  refine ⟨M, c * (A + 1) + 1, hfin, fun n => ?_⟩
  have hc : c ≤ c * (A + 1) + 1 := by nlinarith
  have hz : Zu n + 1 ≤ (A + 1) * (n / D n + 1) := by
    calc Zu n + 1 ≤ A * (n / D n + 1) + 1 := Nat.add_le_add_right (hA n) 1
      _ ≤ A * (n / D n + 1) + (n / D n + 1) :=
        Nat.add_le_add_left (Nat.succ_le_succ (Nat.zero_le (n / D n))) _
      _ = (A + 1) * (n / D n + 1) := by ring
  have he : c * (Zu n + 1) ≤ (c * (A + 1) + 1) * (n / D n + 1) := by
    calc c * (Zu n + 1) ≤ c * ((A + 1) * (n / D n + 1)) :=
        Nat.mul_le_mul_left c hz
      _ = (c * (A + 1)) * (n / D n + 1) := by ring
      _ ≤ (c * (A + 1) + 1) * (n / D n + 1) :=
        Nat.mul_le_mul_right _ (Nat.le_succ _)
  exact le_trans (htime n) (Nat.mul_le_mul hc (two_pow_le he))

end IntegerFactoring.K6Machine

#print axioms IntegerFactoring.K6Machine.value_requirements
#print axioms IntegerFactoring.K6Machine.factoring_Zu_exponential
#print axioms IntegerFactoring.K6Machine.factoring_rate_of_Zu_bound
