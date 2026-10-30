/-
SPDX-License-Identifier: Apache-2.0

Executable extraction of a multiplication-gate code from a cluster-state code.
The gate's internal numeric fields and row summands are not yet certified.
-/
import ClusterSize

open Groundwork Factoring.Aux

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def clusterStateDigit (s g B γ code i : ℕ) : ℕ :=
  cappedDigit (2 ^ B) (gateCount (pA s)) code i (clusterRadix s g γ)

theorem clusterStateDigit_polyT {α : Type} [Enc α] {n g B γ code : α → ℕ}
    (hn : Unary n) (hg : PolyT g) (hB : Unary B) (hγ : PolyT γ) (hc : PolyT code) :
    PolyT fun p : α × ℕ => clusterStateDigit (sOf (n p.1)) (g p.1) (B p.1) (γ p.1) (code p.1) p.2 := by
  have hA := selected_pA_unary hn
  have hcount : Unary fun a => gateCount (pA (sOf (n a))) :=
    Unary.add hA (Unary.add hA (Unary.add hA (Unary.add hA (Unary.mul hA hA))))
  have hr : PolyT fun p : α × ℕ => clusterRadix (sOf (n p.1)) (g p.1) (γ p.1) p.2 :=
    clusterRadix_polyT hn hg hγ
  exact cappedDigit_polyT (B := fun a => 2 ^ B a) (g := fun a => gateCount (pA (sOf (n a))))
    (code := code) (radix := fun a i => clusterRadix (sOf (n a)) (g a) (γ a) i)
    (PolyT.two_pow hB) hcount hc hr

theorem clusterStateDigit_eq (s g B N k y : ℕ) (hg : 0 < g) (hB : Wc s g ≤ B)
    (γ : Fin g) :
    let χ := fun a => clMG s g hg ((MG.enum (pA s)).symm a)
    let κ := mulCodingZ s N k y
    ∀ (c : Fin ((κ.coarsen χ).s γ)) (a : {a : Fin (mulNodes s) // χ a = γ}),
      clusterStateDigit s g B γ.val c.val a.val.val =
        ((κ.e a.val) (((κ.coarsen χ).e γ).symm c a)).val := by
  dsimp only
  let χ := fun a => clMG s g hg ((MG.enum (pA s)).symm a)
  let κ := mulCodingZ s N k y
  intro c a
  have hrad : clusterRadix s g γ.val = Coding.sg κ χ γ := by
    funext i
    exact clusterRadix_eq s g N k y hg γ i
  have hpos : ∀ i, 0 < Coding.sg κ χ γ i := Coding.sg_pos κ χ γ
  have hmono : Monotone (pp (Coding.sg κ χ γ)) := by
    apply monotone_nat_of_le_succ
    intro i
    exact Nat.le_mul_of_pos_right _ (hpos i)
  have hfull : pp (Coding.sg κ χ γ) (gateCount (pA s)) = (κ.coarsen χ).s γ := by
    change pp (Coding.sg κ χ γ) (gateCount (pA s)) =
      if h : γ.val < g then pp (Coding.sg κ χ ⟨γ.val, h⟩) (mulNodes s) else 1
    rw [dif_pos γ.isLt]
    rfl
  have hsize : (κ.coarsen χ).s γ ≤ 2 ^ Wc s g := cluster_s_le s 2 N k y g hg γ
  have ha : a.val.val ≤ gateCount (pA s) := by simpa [gateCount, mulNodes] using Nat.le_of_lt a.val.isLt
  have hprefix : pp (Coding.sg κ χ γ) a.val.val ≤ 2 ^ B :=
    (hmono ha).trans
      (by rw [hfull]; exact hsize.trans (Nat.pow_le_pow_right (by norm_num) hB))
  unfold clusterStateDigit
  rw [cappedDigit_eq _ _ _ _ _ (by simpa [gateCount, mulNodes] using Nat.le_of_lt a.val.isLt) (by
    intro i _
    rw [hrad]
    exact hpos i) (by rw [hrad]; exact hprefix), hrad]
  exact (Coding.coarsen_digit κ χ γ c a).symm

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.clusterStateDigit_polyT
#print axioms IntegerFactoring.HalfAppendix.clusterStateDigit_eq
