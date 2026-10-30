/-
Quantitative arithmetic for a group-count recurrence of degree twelve.

These lemmas are independent of the modified factoring implementation.
All roots and logarithms are natural-number operations; no real analysis,
new axiom, or executable shortcut is used.
-/
module

public import Mathlib.Analysis.SpecialFunctions.Pow.NthRootLemmas
public import Mathlib.Data.Nat.Log
public import Mathlib.Tactic.NormNum
public import Mathlib.Tactic.Ring
public import Mathlib.Tactic.Linarith

@[expose] public section

namespace IntegerFactoring.K6RateArithmetic

/-- The integer denominator scale for a saving of exponent `1/30`. -/
def rateRoot (n : ℕ) : ℕ := Nat.nthRoot 30 (Nat.log 2 n)

theorem rateRoot_pow_le (n : ℕ) : rateRoot n ^ 30 ≤ Nat.log 2 n := by
  exact Nat.pow_nthRoot_le (Or.inl (by norm_num))

theorem rateRoot_le_log (n : ℕ) : rateRoot n ≤ Nat.log 2 n :=
  le_trans (Nat.le_self_pow (by norm_num) _) (rateRoot_pow_le n)

/-- The two finite inequalities separate `1/30` from both reference rates. -/
theorem baseline_rate_comparison : (12 : ℕ) ^ 30 < 16 * 11 ^ 30 := by norm_num

theorem new_rate_comparison : 12 * (23 : ℕ) ^ 30 < 25 ^ 30 := by norm_num

/-- Multiplying a group count by two absorbs the coefficient 25. -/
theorem scaled_groups_bound (K : ℕ → ℕ)
    (hstep : ∀ d, K (d + 1) ≤ 25 * K d ^ 12) :
    ∀ d, 2 * K d ≤ (2 * K 0) ^ (12 ^ d) := by
  intro d
  induction d with
  | zero => simp
  | succ d ih =>
    calc
      2 * K (d + 1) ≤ 2 * (25 * K d ^ 12) := Nat.mul_le_mul_left _ (hstep d)
      _ ≤ (2 * K d) ^ 12 := by
        rw [mul_pow]
        have hc : (50 : ℕ) ≤ 2 ^ 12 := by norm_num
        have hh := Nat.mul_le_mul_right (K d ^ 12) hc
        nlinarith only [hh]
      _ ≤ ((2 * K 0) ^ (12 ^ d)) ^ 12 := Nat.pow_le_pow_left ih _
      _ = (2 * K 0) ^ (12 ^ (d + 1)) := by rw [← pow_mul, pow_succ]

/-- A pure-power upper bound, convenient for the natural logarithm. -/
theorem groups_le_two_pow (K : ℕ → ℕ) (B : ℕ)
    (hstep : ∀ d, K (d + 1) ≤ 25 * K d ^ 12)
    (hbase : 2 * K 0 ≤ 2 ^ B) (d : ℕ) :
    K d ≤ 2 ^ (B * 12 ^ d) := by
  calc
    K d ≤ 2 * K d := by omega
    _ ≤ (2 * K 0) ^ (12 ^ d) := scaled_groups_bound K hstep d
    _ ≤ (2 ^ B) ^ (12 ^ d) := Nat.pow_le_pow_left hbase _
    _ = 2 ^ (B * 12 ^ d) := by rw [pow_mul]

/-- Maximal admissible depth bounds `log₂ n` by a constant times `12^ell`. -/
theorem log_le_depth (K : ℕ → ℕ) (B n ell : ℕ)
    (hstep : ∀ d, K (d + 1) ≤ 25 * K d ^ 12)
    (hbase : 2 * K 0 ≤ 2 ^ B)
    (hupper : n < K (ell + 1) ^ 40) :
    Nat.log 2 n ≤ (480 * B) * 12 ^ ell := by
  have hg := groups_le_two_pow K B hstep hbase (ell + 1)
  have hn : n ≤ 2 ^ ((480 * B) * 12 ^ ell) := by
    calc
      n ≤ K (ell + 1) ^ 40 := hupper.le
      _ ≤ (2 ^ (B * 12 ^ (ell + 1))) ^ 40 := Nat.pow_le_pow_left hg _
      _ = 2 ^ ((480 * B) * 12 ^ ell) := by
        rw [← pow_mul, pow_succ]
        congr 1
        ring
  calc
    Nat.log 2 n ≤ Nat.log 2 (2 ^ ((480 * B) * 12 ^ ell)) := Nat.log_mono_right hn
    _ = (480 * B) * 12 ^ ell := Nat.log_pow (by norm_num) _

/-- The weighted contraction estimate from a direct logarithm bound. -/
theorem root_weighted_of_log_bound (n ell A C : ℕ)
    (hlog : Nat.log 2 n ≤ A * 12 ^ ell) (hA : A ≤ C ^ 30) :
    rateRoot n * 23 ^ ell ≤ C * 25 ^ ell := by
  apply (Nat.pow_le_pow_iff_left (by norm_num : 30 ≠ 0)).mp
  calc
    (rateRoot n * 23 ^ ell) ^ 30
        = rateRoot n ^ 30 * ((23 : ℕ) ^ 30) ^ ell := by
          rw [mul_pow, ← pow_mul, ← pow_mul]
          congr 2
          omega
    _ ≤ (A * 12 ^ ell) * ((23 : ℕ) ^ 30) ^ ell :=
      Nat.mul_le_mul_right _ (le_trans (rateRoot_pow_le n) hlog)
    _ = A * (12 * (23 : ℕ) ^ 30) ^ ell := by rw [mul_pow]; ring
    _ ≤ A * ((25 : ℕ) ^ 30) ^ ell :=
      Nat.mul_le_mul_left _ (Nat.pow_le_pow_left new_rate_comparison.le _)
    _ ≤ C ^ 30 * ((25 : ℕ) ^ 30) ^ ell := Nat.mul_le_mul_right _ hA
    _ = (C * 25 ^ ell) ^ 30 := by
      rw [mul_pow, ← pow_mul, ← pow_mul]
      congr 2
      omega

/-- Generic recurrence-to-rate bridge. -/
theorem root_weighted_bound (K : ℕ → ℕ) (B n ell C : ℕ)
    (hstep : ∀ d, K (d + 1) ≤ 25 * K d ^ 12)
    (hbase : 2 * K 0 ≤ 2 ^ B)
    (hupper : n < K (ell + 1) ^ 40)
    (hC : 480 * B ≤ C ^ 30) :
    rateRoot n * 23 ^ ell ≤ C * 25 ^ ell :=
  root_weighted_of_log_bound n ell (480 * B) C
    (log_le_depth K B n ell hstep hbase hupper) hC

/-- For a starting group count at most 1024, the universal constant is two. -/
theorem root_weighted_bound_small_start (K : ℕ → ℕ) (n ell : ℕ)
    (hstep : ∀ d, K (d + 1) ≤ 25 * K d ^ 12)
    (hbase : K 0 ≤ 1024) (hupper : n < K (ell + 1) ^ 40) :
    rateRoot n * 23 ^ ell ≤ 2 * 25 ^ ell := by
  apply root_weighted_bound K 11 n ell 2 hstep
  · norm_num
    omega
  · exact hupper
  · norm_num

/-- `R+1` is also controlled, without a separate small-input exception. -/
theorem denominator_weighted_bound (n ell C : ℕ)
    (h : rateRoot n * 23 ^ ell ≤ C * 25 ^ ell) :
    (rateRoot n + 1) * 23 ^ ell ≤ (C + 1) * 25 ^ ell := by
  have hp : (23 : ℕ) ^ ell ≤ 25 ^ ell := Nat.pow_le_pow_left (by norm_num) _
  nlinarith only [h, hp]

/-- The weighted estimate survives the floor in the bottom volume. -/
theorem weighted_floor_bound_nat (R ell C N : ℕ)
    (h : R * 23 ^ ell ≤ C * 25 ^ ell) :
    R * (23 ^ ell * N / 25 ^ ell) ≤ C * N := by
  calc
    R * (23 ^ ell * N / 25 ^ ell)
        ≤ R * (23 ^ ell * N) / 25 ^ ell :=
          Nat.mul_div_le_mul_div_assoc _ _ _
    _ ≤ C * N := by
      apply Nat.div_le_of_le_mul
      have hh := Nat.mul_le_mul_right N h
      nlinarith only [hh]

/-- The same bound specialized to the thirtieth-root scale. -/
theorem weighted_floor_bound (n ell C N : ℕ)
    (h : rateRoot n * 23 ^ ell ≤ C * 25 ^ ell) :
    rateRoot n * (23 ^ ell * N / 25 ^ ell) ≤ C * N :=
  weighted_floor_bound_nat (rateRoot n) ell C N h

set_option exponentiation.threshold 1000 in
/-- A simple polynomial-size bound implies a small integer denominator. -/
theorem rateRoot_le_parameter {n P : ℕ} (hP : 48 ≤ P) (hn : n ≤ P ^ 700) :
    rateRoot n ≤ P := by
  have hp : P ≤ 2 ^ P := (Nat.lt_two_pow_self).le
  have hn2 : n ≤ 2 ^ (700 * P) := by
    calc
      n ≤ P ^ 700 := hn
      _ ≤ (2 ^ P) ^ 700 := Nat.pow_le_pow_left hp _
      _ = 2 ^ (700 * P) := by rw [← pow_mul]; congr 1; omega
  have hlog : Nat.log 2 n ≤ 700 * P := by
    calc
      Nat.log 2 n ≤ Nat.log 2 (2 ^ (700 * P)) := Nat.log_mono_right hn2
      _ = 700 * P := Nat.log_pow (by norm_num) _
  have hP2 : 700 ≤ P ^ 2 := by nlinarith
  have hP29 : P ^ 2 ≤ P ^ 29 := Nat.pow_le_pow_right (by omega) (by norm_num)
  have hlarge : 700 * P ≤ P ^ 30 := by
    calc
      700 * P ≤ P ^ 29 * P := Nat.mul_le_mul_right _ (le_trans hP2 hP29)
      _ = P ^ 30 := by ring
  apply (Nat.pow_le_pow_iff_left (by norm_num : 30 ≠ 0)).mp
  exact le_trans (rateRoot_pow_le n) (le_trans hlog hlarge)

set_option exponentiation.threshold 1000 in
/-- The positive final denominator is still at most the parameter. -/
theorem rateRoot_add_one_le_parameter {n P : ℕ}
    (hP : 48 ≤ P) (hn : n ≤ P ^ 700) : rateRoot n + 1 ≤ P := by
  have hp : P ≤ 2 ^ P := (Nat.lt_two_pow_self).le
  have hn2 : n ≤ 2 ^ (700 * P) := by
    calc
      n ≤ P ^ 700 := hn
      _ ≤ (2 ^ P) ^ 700 := Nat.pow_le_pow_left hp _
      _ = 2 ^ (700 * P) := by rw [← pow_mul]; congr 1; omega
  have hlog : Nat.log 2 n ≤ 700 * P := by
    calc
      Nat.log 2 n ≤ Nat.log 2 (2 ^ (700 * P)) := Nat.log_mono_right hn2
      _ = 700 * P := Nat.log_pow (by norm_num) _
  have hP2 : 700 < P ^ 2 := by nlinarith
  have hP29 : P ^ 2 ≤ P ^ 29 := Nat.pow_le_pow_right (by omega) (by norm_num)
  have hlarge : 700 * P < P ^ 30 := by
    calc
      700 * P < P ^ 29 * P :=
        Nat.mul_lt_mul_of_pos_right (lt_of_lt_of_le hP2 hP29) (by omega)
      _ = P ^ 30 := by ring
  have hr : rateRoot n < P := by
    unfold rateRoot
    exact (Nat.nthRoot_lt_iff (by norm_num : 30 ≠ 0)).mpr (lt_of_le_of_lt hlog hlarge)
  omega

end IntegerFactoring.K6RateArithmetic

#print axioms IntegerFactoring.K6RateArithmetic.root_weighted_bound
#print axioms IntegerFactoring.K6RateArithmetic.root_weighted_bound_small_start
#print axioms IntegerFactoring.K6RateArithmetic.rateRoot_le_parameter
#print axioms IntegerFactoring.K6RateArithmetic.rateRoot_add_one_le_parameter
#print axioms IntegerFactoring.K6RateArithmetic.weighted_floor_bound_nat
#print axioms IntegerFactoring.K6RateArithmetic.baseline_rate_comparison
#print axioms IntegerFactoring.K6RateArithmetic.new_rate_comparison
