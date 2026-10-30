/-
SPDX-License-Identifier: Apache-2.0

Exact arithmetic for equations (8)--(9) of research/integration.md.
These lemmas certify widths conditional on the stated network/grouping
bounds. They do not construct a network or an improved factoring machine.
-/
module

public import Mathlib.Tactic.Linarith
public import Mathlib.Tactic.Positivity
public import Mathlib.Tactic.Ring

@[expose] public section

namespace Factoring.Improvements

def orderedTotalBound (grid width : ℕ) : ℕ :=
  281 * grid ^ 4 * width / 312 + grid ^ 3 * width + 3 * grid ^ 4

def commonWidth (groups accuracy grid width : ℕ) : ℕ :=
  (accuracy * orderedTotalBound grid width + grid ^ 4 * width) /
    (accuracy * groups ^ 2)

/-- The floor gives exactly the largest integral width admitted by the bound. -/
theorem common_width_iff (groups accuracy grid width pair : ℕ)
    (hgroups : 0 < groups) (haccuracy : 0 < accuracy) :
    pair ≤ commonWidth groups accuracy grid width ↔
      accuracy * groups ^ 2 * pair ≤
        accuracy * orderedTotalBound grid width + grid ^ 4 * width := by
  unfold commonWidth
  rw [Nat.le_div_iff_mul_le (by positivity : 0 < accuracy * groups ^ 2)]
  simp [Nat.mul_comm]

/-- Any pair obeying the grouping inequality fits the common padded width.
The upper bound depends only on the parameters of the level. -/
theorem grouped_width_le_common (groups accuracy grid width pair total degree : ℕ)
    (hgroups : 0 < groups) (haccuracy : 0 < accuracy)
    (hpair : accuracy * groups ^ 2 * pair ≤ accuracy * total + grid * degree)
    (htotal : total ≤ orderedTotalBound grid width)
    (hdegree : degree ≤ grid ^ 3 * width) :
    pair ≤ commonWidth groups accuracy grid width := by
  have hnum : accuracy * groups ^ 2 * pair ≤
      accuracy * orderedTotalBound grid width + grid ^ 4 * width := by
    calc
      _ ≤ accuracy * total + grid * degree := hpair
      _ ≤ accuracy * orderedTotalBound grid width + grid * (grid ^ 3 * width) :=
        Nat.add_le_add (Nat.mul_le_mul_left _ htotal) (Nat.mul_le_mul_left _ hdegree)
      _ = _ := by ring
  exact (common_width_iff groups accuracy grid width pair hgroups haccuracy).2 hnum

/-- Integral form of the one-level volume recurrence (9), for positive grid:
N_out ≤ (281/312 + 1/grid + 1/accuracy)*N_in + 3*grid^4.
The polynomial inequality also holds when a parameter is zero. -/
theorem common_width_volume (groups accuracy grid width : ℕ) :
    312 * accuracy * (groups ^ 2 * commonWidth groups accuracy grid width) ≤
      281 * accuracy * grid ^ 4 * width + 312 * accuracy * grid ^ 3 * width +
        312 * grid ^ 4 * width + 936 * accuracy * grid ^ 4 := by
  have hpad := Nat.div_mul_le_self
    (accuracy * orderedTotalBound grid width + grid ^ 4 * width) (accuracy * groups ^ 2)
  have hround := Nat.div_mul_le_self (281 * grid ^ 4 * width) 312
  have hpad' := Nat.mul_le_mul_left 312 hpad
  have hround' := Nat.mul_le_mul_left accuracy hround
  unfold commonWidth orderedTotalBound
  unfold orderedTotalBound at hpad
  unfold orderedTotalBound at hpad'
  nlinarith

end Factoring.Improvements

#print axioms Factoring.Improvements.common_width_iff
#print axioms Factoring.Improvements.grouped_width_le_common
#print axioms Factoring.Improvements.common_width_volume
