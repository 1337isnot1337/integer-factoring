/-
SPDX-License-Identifier: Apache-2.0

Uses the coefficient and crossing definitions of Isaac Rudich and
Louis-Martin Rousseau, pinned at c91cc3f4947caca8977e81a5d14244f73b5ae2b1.
This proves a value identity only, not compressed factor evaluation or a
factoring-machine runtime.
-/
module

public import FactoringAux.Pieces
public import Mathlib.Tactic.Linarith

@[expose] public section

namespace Factoring.Improvements

open Factoring.Aux Finset

/-- Scaling both coefficients in every linear factor scales each coefficient
by the same t-fold power. -/
theorem coefV_scale : ∀ (t : ℕ) (c e : Fin t → ℤ) (a : ℤ) (j : ℕ),
    coefV t (fun s => a * c s) (fun s => a * e s) j = a ^ t * coefV t c e j
  | 0, c, e, a, j => by simp [coefV]
  | t + 1, c, e, a, 0 => by
      simp only [coefV]
      change coefV t (fun s => a * Fin.init c s) (fun s => a * Fin.init e s) 0 *
        (a * c (Fin.last t)) = _
      rw [coefV_scale t (Fin.init c) (Fin.init e) a 0]
      rw [pow_succ]
      ring
  | t + 1, c, e, a, j + 1 => by
      simp only [coefV]
      change coefV t (fun s => a * Fin.init c s) (fun s => a * Fin.init e s) (j + 1) *
          (a * c (Fin.last t)) +
        coefV t (fun s => a * Fin.init c s) (fun s => a * Fin.init e s) j *
          (a * e (Fin.last t)) = _
      rw [coefV_scale t (Fin.init c) (Fin.init e) a (j + 1),
          coefV_scale t (Fin.init c) (Fin.init e) a j]
      rw [pow_succ]
      ring

/-- Coefficient built from `(4-3T)R + T(8E-R)` blockwise. -/
def compressedPiece (lam i : ℕ) (x y u v : Bits lam) : ℤ :=
  coefV (tB lam)
    (fun s => 4 * Rx (blk x s) (blk y s) (blk u s) (blk v s))
    (fun s => (8 * (if blk x s = blk v s ∧ blk y s = blk u s then 1 else 0)
      - Rx (blk x s) (blk y s) (blk u s) (blk v s))
      - 3 * Rx (blk x s) (blk y s) (blk u s) (blk v s)) i

/-- The polynomial using the complement has the same coefficients as `4^t A`. -/
theorem compressedPiece_eq_scaled (lam i : ℕ) (x y u v : Bits lam) :
    compressedPiece lam i x y u v = 4 ^ tB lam * pieceA lam i x y u v := by
  unfold compressedPiece pieceA
  have he : (fun s : Fin (tB lam) =>
      (8 * (if blk x s = blk v s ∧ blk y s = blk u s then (1 : ℤ) else 0)
        - Rx (blk x s) (blk y s) (blk u s) (blk v s))
        - 3 * Rx (blk x s) (blk y s) (blk u s) (blk v s)) =
      (fun s => 4 * Cx (blk x s) (blk y s) (blk u s) (blk v s)) := by
    funext s
    have h := crossing (blk x s) (blk y s) (blk u s) (blk v s)
    split_ifs at h ⊢ <;> linarith
  rw [he]
  exact coefV_scale _ _ _ _ _

/-- Summing all scaled pieces restores crossing equality with multiplier `8^t`. -/
theorem sum_compressedPiece (lam : ℕ) (x y u v : Bits lam) :
    ∑ i ∈ range (tB lam + 1), compressedPiece lam i x y u v =
      8 ^ tB lam * if x = v ∧ y = u then 1 else 0 := by
  simp_rw [compressedPiece_eq_scaled, ← Finset.mul_sum, sum_pieceA]
  rw [← mul_assoc, ← mul_pow]
  norm_num

/-- The scaled pieces retain the same absolute-value envelope. -/
theorem abs_compressedPiece_le (lam i : ℕ) (x y u v : Bits lam) :
    |compressedPiece lam i x y u v| ≤ 8 ^ tB lam := by
  rw [compressedPiece_eq_scaled, abs_mul, abs_of_nonneg (by positivity)]
  calc
    (4 : ℤ) ^ tB lam * |pieceA lam i x y u v| ≤ 4 ^ tB lam * 2 ^ tB lam :=
      mul_le_mul_of_nonneg_left (abs_pieceA_le lam i x y u v) (by positivity)
    _ = 8 ^ tB lam := by rw [← mul_pow]; norm_num

#print axioms coefV_scale
#print axioms compressedPiece_eq_scaled
#print axioms sum_compressedPiece
#print axioms abs_compressedPiece_le

end Factoring.Improvements
