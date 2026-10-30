/-
SPDX-License-Identifier: Apache-2.0

Numeric, capped mixed-radix sizes for one multiplication-network cluster.
-/
import GateSize
import CappedProducts
import ExplicitCoding
import ClusterEdge

open Groundwork Factoring.Aux Levels.Gen

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def clusterRadix (s g γ i : ℕ) : ℕ :=
  let A := pA s
  if i < gateCount A ∧ gateCluster A g i = γ then
    numericGateSize A (pD s) (pZ s) (pM s) (fw s) i
  else 1

def cappedClusterSize (s g B γ : ℕ) : ℕ :=
  cappedProduct (2 ^ B) (clusterRadix s g γ) (gateCount (pA s))

theorem clusterRadix_polyT {α : Type} [Enc α] {n g γ : α → ℕ}
    (hn : Unary n) (hg : PolyT g) (hγ : PolyT γ) :
    PolyT fun p : α × ℕ => clusterRadix (sOf (n p.1)) (g p.1) (γ p.1) p.2 := by
  have hA : PolyT fun p : α × ℕ => pA (sOf (n p.1)) := PolyT.comp (f := fun a => pA (sOf (n a))) (selected_pA_unary hn).poly PolyT.id'.fst
  have hD : PolyT fun p : α × ℕ => pD (sOf (n p.1)) := PolyT.comp (f := fun a => pD (sOf (n a))) (selected_pD_polyT hn) PolyT.id'.fst
  have hZ : PolyT fun p : α × ℕ => pZ (sOf (n p.1)) := PolyT.comp (f := fun a => pZ (sOf (n a))) (selected_pZ_polyT hn) PolyT.id'.fst
  have hM : PolyT fun p : α × ℕ => pM (sOf (n p.1)) := PolyT.comp (f := fun a => pM (sOf (n a))) (selected_pM_polyT hn) PolyT.id'.fst
  have hF : Unary fun p : α × ℕ => fw (sOf (n p.1)) := Unary.comp (selected_fw_unary hn) PolyT.id'.fst
  have hg' : PolyT fun p : α × ℕ => g p.1 := PolyT.comp hg PolyT.id'.fst
  have hγ' : PolyT fun p : α × ℕ => γ p.1 := PolyT.comp hγ PolyT.id'.fst
  have hcl := gateCluster_polyT hA hg' (PolyT.id'.snd : PolyT fun p : α × ℕ => p.2)
  have hsz := numericGateSize_polyT hA hD hZ hM hF (PolyT.id'.snd : PolyT fun p : α × ℕ => p.2)
  have hN : PolyT fun p : α × ℕ => gateCount (pA (sOf (n p.1))) :=
    PolyT.add hA (PolyT.add hA (PolyT.add hA (PolyT.add hA (PolyT.mul hA hA))))
  have hlt := PolyT.decide_lt (PolyT.id'.snd : PolyT fun p : α × ℕ => p.2) hN
  have heq := PolyT.decide_eq hcl hγ'
  exact PolyT.ite' (PolyT.decide_and hlt heq) hsz (PolyT.const 1)

theorem cappedClusterSize_polyT {α : Type} [Enc α] {n g B γ : α → ℕ}
    (hn : Unary n) (hg : PolyT g) (hB : Unary B) (hγ : PolyT γ) :
    PolyT fun a => cappedClusterSize (sOf (n a)) (g a) (B a) (γ a) := by
  have hA := selected_pA_unary hn
  have hN : Unary fun a => gateCount (pA (sOf (n a))) :=
    Unary.add hA (Unary.add hA (Unary.add hA (Unary.add hA (Unary.mul hA hA))))
  have hr : PolyT fun p : α × ℕ => clusterRadix (sOf (n p.1)) (g p.1) (γ p.1) p.2 :=
    clusterRadix_polyT hn hg hγ
  exact cappedProduct_polyT (B := fun a => 2 ^ B a) (n := fun a => gateCount (pA (sOf (n a))))
    (f := fun a i => clusterRadix (sOf (n a)) (g a) (γ a) i) (PolyT.two_pow hB) hN hr

theorem clusterRadix_eq (s g N k y : ℕ) (hg : 0 < g) (γ : Fin g) (i : ℕ) :
    clusterRadix s g γ.val i = Coding.sg (mulCodingZ s N k y)
      (fun a => clMG s g hg ((MG.enum (pA s)).symm a)) γ i := by
  by_cases hi : i < gateCount (pA s)
  · have hi' : i < pA s + (pA s + (pA s + (pA s + pA s * pA s))) := by simpa [gateCount] using hi
    let a : Fin (pA s + (pA s + (pA s + (pA s + pA s * pA s)))) := ⟨i, hi'⟩
    let j := (MG.enum (pA s)).symm a
    have hcode : (MG.enum (pA s) j).val = i := by
      exact congrArg Fin.val (Equiv.apply_symm_apply (MG.enum (pA s)) a)
    have hcl : gateCluster (pA s) g i = (clMG s g hg j).val := by
      rw [← hcode]
      exact gateCluster_eq s g hg j
    have hsz : numericGateSize (pA s) (pD s) (pZ s) (pM s) (fw s) i = mulSz s i := by
      exact numericGateSize_eq s a
    simp [clusterRadix, Coding.sg, hi', gateCount, hcl, hsz, Fin.ext_iff, mulCodingZ, a, j]
  · have hi' : ¬ i < pA s + (pA s + (pA s + (pA s + pA s * pA s))) := by simpa [gateCount] using hi
    simp [clusterRadix, Coding.sg, hi', gateCount]

theorem cappedClusterSize_eq (s g B N k y : ℕ) (hg : 0 < g) (γ : Fin g) (hB : Wc s g ≤ B) :
    cappedClusterSize s g B γ.val = (clusterCodingZ s g N k y hg).s γ := by
  have hrad : clusterRadix s g γ.val = Coding.sg (mulCodingZ s N k y)
      (fun a => clMG s g hg ((MG.enum (pA s)).symm a)) γ := by
    funext i
    exact clusterRadix_eq s g N k y hg γ i
  have hpp : pp (clusterRadix s g γ.val) (gateCount (pA s)) = (clusterCodingZ s g N k y hg).s γ := by
    rw [hrad]
    simp [clusterCodingZ, Coding.coarsen, γ.isLt, gateCount, mulNodes]
  have hsize : (clusterCodingZ s g N k y hg).s γ ≤ 2 ^ Wc s g := by
    have h := cluster_s_le s 2 N k y g hg γ
    exact h
  have hbound : ((List.range (gateCount (pA s))).map (clusterRadix s g γ.val)).prod ≤ 2 ^ B := by
    rw [range_prod_eq_pp, hpp]
    exact hsize.trans (Nat.pow_le_pow_right (by norm_num) hB)
  unfold cappedClusterSize
  rw [cappedProduct_eq _ _ _ (by
    intro i _
    rw [hrad]
    exact Coding.sg_pos _ _ _ _) hbound, range_prod_eq_pp, hpp]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.cappedClusterSize_polyT
#print axioms IntegerFactoring.HalfAppendix.clusterRadix_eq
#print axioms IntegerFactoring.HalfAppendix.cappedClusterSize_eq
