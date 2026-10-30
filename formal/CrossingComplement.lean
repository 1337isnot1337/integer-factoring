/-
SPDX-License-Identifier: Apache-2.0

An integer rank-seven factorization of a signed 8 by 8 complement block.
Eight disjoint XOR blocks give the 56-channel construction in
research/crossing-compression.md. This file proves the block factorization
and the scalar crossing identity; it does not formalize the full recursion.
-/
module

public import FactoringAux.Cross
public import Mathlib.Tactic.Ring
public import Mathlib.Tactic.Linarith

@[expose] public section

namespace Factoring.Improvements

open Finset Factoring.Aux

def complementLeft (s : Fin 8 → ℤ) (x : Fin 8) (a : Fin 7) : ℤ :=
  (if x = a.succ then 1 else 0) - (if x = 0 then s a.succ else 0)

def complementRight (s : Fin 8 → ℤ) (a : Fin 7) (v : Fin 8) : ℤ :=
  8 * (if a.succ = v then 1 else 0) - s a.succ * s v

/-- All coefficients are integers. The seven channels are indexed by `Fin 7`. -/
theorem signed_complement_factorization (s : Fin 8 → ℤ) (h0 : s 0 = 1)
    (hs : ∀ a, s a * s a = 1) (x v : Fin 8) :
    ∑ a : Fin 7, complementLeft s x a * complementRight s a v =
      8 * (if x = v then 1 else 0) - s x * s v := by
  have hzero (a : Fin 7) : (0 : Fin 8) ≠ a.succ := (Fin.succ_ne_zero a).symm
  refine Fin.cases ?_ (fun x => ?_) x
  · refine Fin.cases ?_ (fun v => ?_) v
    · simp [complementLeft, complementRight, hzero, Fin.succ_ne_zero, h0, neg_mul, hs]
    · have hterm (a : Fin 7) :
          (-s a.succ) * (8 * (if a = v then 1 else 0) - s a.succ * s v.succ) =
            -8 * (if a = v then s a.succ else 0) + s v.succ := by
        have h : s a.succ * s a.succ * s v.succ = s v.succ := by
          rw [hs, one_mul]
        split_ifs <;> nlinarith
      simp only [complementLeft, complementRight, hzero,
        Fin.succ_inj, ite_false, ite_true, zero_sub]
      simp_rw [hterm]
      simp [Finset.sum_add_distrib, h0]
      ring
  · simp [complementLeft, complementRight, Fin.succ_ne_zero, Fin.succ_inj,
      ite_mul, eq_comm]

/-- The new residual is an integer linear combination of the original two
crossing tensors. No inverse is introduced in any coefficient ring. -/
theorem crossing_complement_identity (x y u v : Blk) :
    8 * (if x = v ∧ y = u then (1 : ℤ) else 0) - Rx x y u v =
      3 * Rx x y u v + 4 * Cx x y u v := by
  have h := crossing x y u v
  split_ifs at h ⊢ <;> linarith

/-- Exact rational envelope used by the recurrence certificate. -/
theorem seven_pow_26_lt_two_pow_73 : (7 : ℕ) ^ 26 < 2 ^ 73 := by
  norm_num

/-- Explicit little-endian encoding of the eight three-bit blocks. -/
def blockOfFin (a : Fin 8) : Blk := fun j => a.val.testBit j.val

theorem blockOfFin_bijective : Function.Bijective blockOfFin := by
  decide +kernel

/-- Explicit bijection, used to choose the seven nonzero positions. -/
noncomputable def blockEquiv : Fin 8 ≃ Blk := Equiv.ofBijective blockOfFin blockOfFin_bijective

theorem blockOfFin_zero : blockOfFin 0 = zeroBits 3 := by
  decide +kernel

/-- The reference sign on the all-zero left block is one. -/
theorem sign_zero_left (k : Blk) : sgn (zeroBits 3) k = 1 := by
  simp [sgn, phiB, zeroBits]

/-- Solving an XOR class fixes the second endpoint. -/
theorem bxor_solve (x y k : Blk) (h : bxor x y = k) : y = bxor x k := by
  rw [blk_eq_iff] at h ⊢
  simp only [bxor] at h ⊢
  obtain ⟨h0, h1, h2⟩ := h
  constructor
  · cases hx : x 0 <;> cases hy : y 0 <;> cases hk : k 0 <;> simp_all
  constructor
  · cases hx : x 1 <;> cases hy : y 1 <;> cases hk : k 1 <;> simp_all
  · cases hx : x 2 <;> cases hy : y 2 <;> cases hk : k 2 <;> simp_all

/-- Signs on a fixed XOR class, indexed by the explicit block enumeration. -/
def classSign (k : Blk) (a : Fin 8) : ℤ :=
  sgn (blockOfFin a) (bxor (blockOfFin a) k)

def classLeft (k : Blk) (a : Fin 7) (x y : Blk) : ℤ :=
  if bxor x y = k then
    (if x = blockOfFin a.succ then 1 else 0) -
      (if x = blockOfFin 0 then classSign k a.succ else 0)
  else 0

def classRight (k : Blk) (a : Fin 7) (u v : Blk) : ℤ :=
  if bxor v u = k then
    8 * (if blockOfFin a.succ = v then 1 else 0) -
      classSign k a.succ * sgn v u
  else 0

theorem bxor_zero_left (k : Blk) : bxor (zeroBits 3) k = k := by
  funext i
  cases h : k i <;> simp [bxor, zeroBits, h]

theorem classSign_zero (k : Blk) : classSign k 0 = 1 := by
  simp [classSign, blockOfFin_zero, bxor_zero_left, sign_zero_left]

/-- Each of the eight XOR classes has seven integer channels. -/
theorem class_factorization (k x y u v : Blk) :
    ∑ a : Fin 7, classLeft k a x y * classRight k a u v =
      if bxor x y = k ∧ bxor v u = k then
        8 * (if x = v then 1 else 0) - sgn x y * sgn v u
      else 0 := by
  by_cases hxy : bxor x y = k
  · by_cases hvu : bxor v u = k
    · obtain ⟨ix, rfl⟩ := blockOfFin_bijective.2 x
      obtain ⟨iv, rfl⟩ := blockOfFin_bijective.2 v
      have hinj (a b : Fin 8) : blockOfFin a = blockOfFin b ↔ a = b :=
        blockOfFin_bijective.1.eq_iff
      have hx : classSign k ix = sgn (blockOfFin ix) y := by
        rw [classSign, ← bxor_solve (blockOfFin ix) y k hxy]
      have hv : classSign k iv = sgn (blockOfFin iv) u := by
        rw [classSign, ← bxor_solve (blockOfFin iv) u k hvu]
      have h := signed_complement_factorization (classSign k) (classSign_zero k)
        (fun a => sgn_mul_self _ _) ix iv
      simp only [complementLeft, complementRight] at h
      simp only [classLeft, classRight, hxy, hvu, ite_true, hinj, ← hx, ← hv]
      exact h
    · simp [classRight, hvu]
  · simp [classLeft, hxy]

/-- Summing the eight classes gives an explicit 56-channel integer
factorization of the complement crossing tensor. -/
theorem complement_56_factorization (x y u v : Blk) :
    ∑ k : Blk, ∑ a : Fin 7, classLeft k a x y * classRight k a u v =
      8 * (if x = v ∧ y = u then 1 else 0) - Rx x y u v := by
  simp_rw [class_factorization]
  rw [Rx_eq]
  by_cases h : bxor x y = bxor v u
  · have heq : (x = v ∧ y = u) ↔ x = v := by
      constructor
      · exact And.left
      · intro hx
        subst v
        have hy := bxor_solve x y (bxor x y) rfl
        have hu := bxor_solve x u (bxor x y) h.symm
        exact ⟨rfl, hy.trans hu.symm⟩
    simp only [heq, h, ite_true]
    rw [← h]
    simp only [and_self, Finset.sum_ite_eq, mem_univ, ite_true]
  · have hnone : ∀ k : Blk, ¬(bxor x y = k ∧ bxor v u = k) := by
      intro k hk
      exact h (hk.1.trans hk.2.symm)
    have hneq : ¬(x = v ∧ y = u) := by
      rintro ⟨rfl, rfl⟩
      exact h rfl
    simp [hnone, h, hneq]

/-- A computable code for a block, in the same order as `blockOfFin`. -/
def blockCode (b : Blk) : Fin 8 :=
  finFunctionFinEquiv (fun i => finTwoEquiv.symm (b i))

theorem blockCode_blockOfFin : ∀ a : Fin 8, blockCode (blockOfFin a) = a := by
  decide +kernel

theorem blockOfFin_blockCode : ∀ b : Blk, blockOfFin (blockCode b) = b := by
  decide +kernel

def blockCodeEquiv : Blk ≃ Fin 8 where
  toFun := blockCode
  invFun := blockOfFin
  left_inv := blockOfFin_blockCode
  right_inv := blockCode_blockOfFin

/-- The eighth channel is zero, so seven channels fit in three bits. -/
def paddedClassLeft (k : Blk) (a : Fin 8) (x y : Blk) : ℤ :=
  if h : a.val = 0 then 0 else classLeft k ⟨a.val - 1, by omega⟩ x y

def paddedClassRight (k : Blk) (a : Fin 8) (u v : Blk) : ℤ :=
  if h : a.val = 0 then 0 else classRight k ⟨a.val - 1, by omega⟩ u v

/-- A six-bit `Fac`: three bits for the XOR class, three for the padded
seven-way channel. All factor computations are explicit. -/
def complementFac : Fac ℤ (Blk × Blk) (Blk × Blk) 6 where
  l p θ :=
    let channel := (Fin.appendEquiv 3 3).symm θ
    paddedClassLeft channel.1 (blockCode channel.2) p.1 p.2
  r p θ :=
    let channel := (Fin.appendEquiv 3 3).symm θ
    paddedClassRight channel.1 (blockCode channel.2) p.1 p.2

/-- The six-bit interface evaluates to the complement; eight of its 64
binary labels are zero, leaving exactly 56 potentially active channels. -/
theorem eval_complementFac (x y u v : Blk) :
    complementFac.eval (x, y) (u, v) =
      8 * (if x = v ∧ y = u then 1 else 0) - Rx x y u v := by
  unfold Fac.eval complementFac
  rw [← (Fin.appendEquiv 3 3).sum_comp, Fintype.sum_prod_type]
  simp only [Equiv.symm_apply_apply]
  calc
    (∑ k : Blk, ∑ b : Blk,
        paddedClassLeft k (blockCode b) x y * paddedClassRight k (blockCode b) u v) =
      ∑ k : Blk, ∑ a : Fin 7, classLeft k a x y * classRight k a u v := by
        apply Finset.sum_congr rfl
        intro k _
        rw [← blockCodeEquiv.symm.sum_comp]
        simp only [blockCodeEquiv]
        rw [Fin.sum_univ_succ]
        simp [paddedClassLeft, paddedClassRight, blockCode_blockOfFin]
    _ = _ := complement_56_factorization x y u v

theorem eval_complementFac_linear (x y u v : Blk) :
    complementFac.eval (x, y) (u, v) = 3 * Rx x y u v + 4 * Cx x y u v := by
  rw [eval_complementFac, crossing_complement_identity]

#print axioms signed_complement_factorization
#print axioms crossing_complement_identity
#print axioms seven_pow_26_lt_two_pow_73
#print axioms blockOfFin_bijective
#print axioms blockOfFin_zero
#print axioms sign_zero_left
#print axioms bxor_solve
#print axioms class_factorization
#print axioms complement_56_factorization
#print axioms blockCode_blockOfFin
#print axioms blockOfFin_blockCode
#print axioms eval_complementFac
#print axioms eval_complementFac_linear

end Factoring.Improvements
