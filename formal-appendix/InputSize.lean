/-
SPDX-License-Identifier: Apache-2.0

Arithmetic bounds used when charging the written evaluator input. These
lemmas do not yet give a complete input-length or machine-step bound.
-/
import SelectedPadCost
import WrittenPoly

open Levels.Gen Factoring.Aux IntegerFactoring.HalfPad
open Groundwork PaddedDescent.F

namespace IntegerFactoring.HalfAppendix

set_option maxRecDepth 10000

theorem le_div_envelope (a n r M : ℕ) (hr : 0 < r) (h : r * a ≤ M * n) :
    a ≤ M * (n / r + 1) := by
  by_cases hM : M = 0
  · subst M
    have : a = 0 := by nlinarith
    simp [this]
  · have hMpos : 0 < M := by omega
    have hmod := Nat.mod_lt n hr
    have hdiv := Nat.mod_add_div n r
    have hn : n < r * (n / r + 1) := by nlinarith
    have hmul := Nat.mul_lt_mul_of_pos_left hn hMpos
    have hlt : r * a < r * (M * (n / r + 1)) := by nlinarith
    exact (Nat.lt_of_mul_lt_mul_left hlt).le

theorem selected_table_exponents (P n r : ℕ) (hP : 2 ≤ P)
    (hn : P ^ 40 ≤ n) (hr : 0 < r) (hrP : r ≤ P) :
    (P - 1) * (topWidth (sOf n) P + 1) ≤ 6000 * (n / r + 1) ∧
      topEntryBits (sOf n) P ≤ 8000 * (n / r + 1) := by
  have h := selected_written_exponents P n r hP hn hrP
  exact ⟨le_div_envelope _ n r 6000 hr h.1, le_div_envelope _ n r 8000 hr h.2⟩

theorem n_le_exp_of_parameters (P n r : ℕ) (hP : 2 ≤ P) (hr : 0 < r)
    (hrP : r ≤ P) (hnlo : P ^ 40 ≤ n) (hnhi : n ≤ P ^ 400) :
    n ≤ 2 ^ (400 * (n / r + 1)) := by
  have hP2 : P ^ 2 ≤ n := (Nat.pow_le_pow_right (by omega : 0 < P) (by omega : 2 ≤ 40)).trans hnlo
  have hPr : r * P ≤ n := by
    calc
      r * P ≤ P * P := Nat.mul_le_mul_right P hrP
      _ = P ^ 2 := by ring
      _ ≤ n := hP2
  have hPq : P ≤ n / r + 1 := by
    have h := le_div_envelope P n r 1 hr (by simpa using hPr)
    simpa using h
  calc
    n ≤ P ^ 400 := hnhi
    _ ≤ (2 ^ P) ^ 400 := Nat.pow_le_pow_left (Nat.lt_two_pow_self.le) 400
    _ = 2 ^ (400 * P) := by rw [← pow_mul]; congr 1; ring
    _ ≤ 2 ^ (400 * (n / r + 1)) := Nat.pow_le_pow_right (by norm_num) (by nlinarith)

/-- The encoded instance-length envelope is exponential in O(n/r), provided
the selected width and entry bounds hold. This bounds input size only. -/
theorem instLen_le_exp (P L B n r : ℕ) (hP : 2 ≤ P) (hr : 0 < r) (hrP : r ≤ P)
    (hnlo : P ^ 40 ≤ n) (hnhi : n ≤ P ^ 400)
    (hwidth : P ^ 2 * L ≤ 6000 * n) (hbits : P * B ≤ 8000 * n) :
    PaddedDescent.F.instLen P L B ≤ 2 ^ (7000 * (n / r + 1)) := by
  let q := n / r + 1
  let Q := (P - 1) * L
  let u := 2 ^ Q
  have hq : 1 ≤ q := Nat.le_add_left 1 (n / r)
  have hn1 : 1 ≤ n := by
    have hp : 0 < P ^ 40 := by positivity
    omega
  have hP2 : P ^ 2 ≤ n := (Nat.pow_le_pow_right (by omega : 0 < P) (by omega : 2 ≤ 40)).trans hnlo
  have hP_n : P ≤ n := (Nat.le_self_pow (by norm_num : 40 ≠ 0) P).trans hnlo
  have hPr : r * P ≤ n := by
    calc
      r * P ≤ P * P := Nat.mul_le_mul_right P hrP
      _ = P ^ 2 := by ring
      _ ≤ n := hP2
  have hPq : P ≤ q := by
    have h := le_div_envelope P n r 1 hr (by simpa using hPr)
    simpa [q] using h
  have hQmul : r * Q ≤ 6000 * n := by
    calc
      r * Q ≤ P * (P * L) := Nat.mul_le_mul hrP (Nat.mul_le_mul_right L (Nat.sub_le P 1))
      _ = P ^ 2 * L := by ring
      _ ≤ 6000 * n := hwidth
  have hQ : Q ≤ 6000 * q := by simpa [q] using le_div_envelope Q n r 6000 hr hQmul
  have hL : L ≤ 6000 * n := by
    have h := Nat.mul_le_mul_right L (show 1 ≤ P ^ 2 from Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by omega)))
    simpa using h.trans hwidth
  have hB : B ≤ 8000 * n := by
    have h := Nat.mul_le_mul_right B (show 1 ≤ P by omega)
    simpa using h.trans hbits
  have hu : 1 ≤ u := by dsimp [u]; exact Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
  have hpoly : PaddedDescent.F.instLen P L B ≤ (50000 * (n + 1) ^ 2) * u := by
    have hc : 4 * L + P * (3 * B + 14) + 10 ≤ 50000 * (n + 1) ^ 2 := by nlinarith
    have hm := Nat.mul_le_mul_right u hc
    unfold PaddedDescent.F.instLen
    dsimp [Q, u] at *
    nlinarith
  have hPexp : P ^ 400 ≤ 2 ^ (400 * q) := by
    calc
      P ^ 400 ≤ (2 ^ P) ^ 400 := Nat.pow_le_pow_left (Nat.lt_two_pow_self.le) 400
      _ = 2 ^ (400 * P) := by rw [← pow_mul]; congr 1; ring
      _ ≤ 2 ^ (400 * q) := Nat.pow_le_pow_right (by norm_num) (by omega)
  have hnexp : n ≤ 2 ^ (400 * q) := hnhi.trans hPexp
  have hnplus : n + 1 ≤ 2 ^ (400 * q + 1) := by
    have hp : 1 ≤ 2 ^ (400 * q) := Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
    rw [pow_add]
    nlinarith
  have hsq : (n + 1) ^ 2 ≤ 2 ^ (800 * q + 2) := by
    calc
      (n + 1) ^ 2 ≤ (2 ^ (400 * q + 1)) ^ 2 := Nat.pow_le_pow_left hnplus 2
      _ = 2 ^ (800 * q + 2) := by rw [← pow_mul]; congr 1; omega
  have hpolyexp : 50000 * (n + 1) ^ 2 ≤ 2 ^ (800 * q + 18) := by
    calc
      50000 * (n + 1) ^ 2 ≤ 2 ^ 16 * 2 ^ (800 * q + 2) := Nat.mul_le_mul (by norm_num) hsq
      _ = 2 ^ (800 * q + 18) := by rw [pow_add]; ring
  calc
    PaddedDescent.F.instLen P L B ≤ (50000 * (n + 1) ^ 2) * u := hpoly
    _ ≤ 2 ^ (800 * q + 18) * 2 ^ (6000 * q) :=
      Nat.mul_le_mul hpolyexp (Nat.pow_le_pow_right (by norm_num) hQ)
    _ = 2 ^ ((800 * q + 18) + 6000 * q) := (pow_add 2 _ _).symm
    _ = 2 ^ (6800 * q + 18) := congrArg (fun t : ℕ => 2 ^ t) (by omega)
    _ ≤ 2 ^ (7000 * q) := Nat.pow_le_pow_right (by norm_num) (by omega)

/-- Charge the evaluator's entire encoded input, including both unary pads
and the instance. The hypotheses can be supplied by the schedule estimates;
this remains an input-size theorem, not a constructor runtime theorem. -/
theorem evaluator_input_len_le_exp (q A p l : ℕ) (I : PInst)
    (hq : 1 ≤ q) (hp : p ≤ A * q)
    (hl : l ≤ 2 ^ (7000 * q)) (hI : len I ≤ 2 ^ (7000 * q)) :
    len (List.replicate (2 ^ p) (), List.replicate l (), I) ≤
      2 ^ ((A + 7005) * q) := by
  let X := 2 ^ ((A + 7000) * q)
  have hX : 1 ≤ X := by dsimp [X]; exact Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by norm_num))
  have hpX : 2 ^ p ≤ X := by
    dsimp [X]
    exact Nat.pow_le_pow_right (by norm_num) (by nlinarith)
  have hiX : 2 ^ (7000 * q) ≤ X := by
    dsimp [X]
    exact Nat.pow_le_pow_right (by norm_num) (by nlinarith)
  have hlen : len (List.replicate (2 ^ p) (), List.replicate l (), I) ≤ 32 * X := by
    simp only [len_pair, len_units, List.length_replicate]
    omega
  calc
    _ ≤ 32 * X := hlen
    _ = 2 ^ (5 + (A + 7000) * q) := by dsimp [X]; rw [pow_add]; norm_num
    _ ≤ 2 ^ ((A + 7005) * q) := Nat.pow_le_pow_right (by norm_num) (by nlinarith)

/-- Full evaluator-input length at a level satisfying the selection
inequalities, for any instance obeying the source length envelope. -/
theorem scheduled_evaluator_input_len (ρ : ℝ) (x l n C : ℕ) (I : PInst)
    (hx : 1 ≤ x) (ha : AppRunningTime.PartA ρ x)
    (hnlo : (K x l) ^ 40 ≤ n) (hnhi : n ≤ (K x (l + 1)) ^ 40)
    (hI : len I + 2 ≤ PaddedDescent.F.instLen (K x l)
      (topWidth (sOf n) (K x l) + 1) (topEntryBits (sOf n) (K x l))) :
    let r := 2 ^ (l / 2)
    let q := n / r + 1
    let p := padExponent (K x l) (topWidth (sOf n) (K x l) + 1)
      (topEntryBits (sOf n) (K x l)) C r
    len (List.replicate (2 ^ p) (), List.replicate l (), I) ≤
      2 ^ ((13006 + C * 500) * q) := by
  let P := K x l
  let L := topWidth (sOf n) P + 1
  let B := topEntryBits (sOf n) P
  let r := 2 ^ (l / 2)
  let q := n / r + 1
  let p := padExponent P L B C r
  have hP : 2 ≤ P := K_two x l hx
  have hr : 0 < r := by dsimp [r]; positivity
  have hrP : r ≤ P := half_level_le_K x l hx
  have hhi : n ≤ P ^ 400 := by
    calc
      n ≤ (K x (l + 1)) ^ 40 := hnhi
      _ ≤ (P ^ 10) ^ 40 := Nat.pow_le_pow_left (ha.K_lt_le l).2 40
      _ = P ^ 400 := by rw [← pow_mul]
  have hw : P ^ 2 * L ≤ 6000 * n := selected_width_le P n hP hnlo
  have hb : P * B ≤ 8000 * n := selected_state_exponent_le P n hP hnlo
  have hinst : PaddedDescent.F.instLen P L B ≤ 2 ^ (7000 * q) :=
    instLen_le_exp P L B n r hP hr hrP hnlo hhi hw hb
  change len I + 2 ≤ PaddedDescent.F.instLen P L B at hI
  have hi : len I ≤ 2 ^ (7000 * q) := by omega
  have hP_n : P ≤ n := (Nat.le_self_pow (by norm_num : 40 ≠ 0) P).trans hnlo
  have hl_n : l ≤ n := by
    have hpow := two_pow_le_K x l hx
    have hlt := Nat.lt_two_pow_self (n := l)
    omega
  have hnexp := n_le_exp_of_parameters P n r hP hr hrP hnlo hhi
  have hl : l ≤ 2 ^ (7000 * q) := by
    apply hl_n.trans
    apply hnexp.trans
    apply Nat.pow_le_pow_right (by norm_num)
    omega
  have hpmul := scheduled_pad_cost ρ x l n C hx ha hnlo hnhi
  have hp : p ≤ (6001 + C * 500) * q :=
    le_div_envelope p n r (6001 + C * 500) hr hpmul
  have hq : 1 ≤ q := Nat.le_add_left 1 (n / r)
  have h := evaluator_input_len_le_exp q (6001 + C * 500) p l I hq hp hl hi
  simpa [P, L, B, r, q, p, Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using h

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.le_div_envelope
#print axioms IntegerFactoring.HalfAppendix.selected_table_exponents
#print axioms IntegerFactoring.HalfAppendix.n_le_exp_of_parameters
#print axioms IntegerFactoring.HalfAppendix.instLen_le_exp
#print axioms IntegerFactoring.HalfAppendix.evaluator_input_len_le_exp
#print axioms IntegerFactoring.HalfAppendix.scheduled_evaluator_input_len
