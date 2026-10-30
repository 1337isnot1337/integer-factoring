/-
SPDX-License-Identifier: Apache-2.0

A bounded-index polynomial-time numeric implementation of a modular
transform stage. Gate-specific data sources and full messages remain separate.
-/
import GateFields
import FactoringAux.Transform

open Groundwork Factoring.Aux Finset

namespace IntegerFactoring.HalfAppendix

set_option maxHeartbeats 2000000

def numericHat (A M root k : ℕ) (f : ℕ → ℕ) : ℕ :=
  ((List.range A).map fun i => f i * root ^ (A * min i A * min k A)).sum % M

theorem numericHat_polyT {α : Type} [Enc α] {A : α → ℕ} {M root k : α → ℕ} {f : α → ℕ → ℕ}
    (hA : Unary A) (hM : PolyT M) (hroot : PolyT root) (hk : PolyT k)
    (hf : PolyT fun p : α × ℕ => f p.1 p.2) :
    PolyT fun a => numericHat (A a) (M a) (root a) (k a) (f a) := by
  let β := α × ℕ
  have hA' : Unary fun p : β => A p.1 := Unary.comp hA PolyT.id'.fst
  have hk' : PolyT fun p : β => k p.1 := PolyT.comp hk PolyT.id'.fst
  have hi : PolyT fun p : β => min p.2 (A p.1) := PolyT.min' PolyT.id'.snd hA'.poly
  have hiu : Unary fun p : β => min p.2 (A p.1) := Unary.of_le hA' hi (fun _ => Nat.min_le_right _ _)
  have hku : Unary fun p : β => min (k p.1) (A p.1) :=
    Unary.of_le hA' (PolyT.min' hk' hA'.poly) (fun _ => Nat.min_le_right _ _)
  have hexp : Unary fun p : β => A p.1 * min p.2 (A p.1) * min (k p.1) (A p.1) :=
    Unary.mul (Unary.mul hA' hiu) hku
  have hroot' : PolyT fun p : β => root p.1 := PolyT.comp hroot PolyT.id'.fst
  have hp := PolyT.pow_unary hroot' hexp
  have ht : PolyT fun p : β => f p.1 p.2 * root p.1 ^ (A p.1 * min p.2 (A p.1) * min (k p.1) (A p.1)) :=
    PolyT.mul hf hp
  have hl : PolyT fun a => (List.range (A a)).map (fun i => f a i * root a ^ (A a * min i (A a) * min (k a) (A a))) :=
    PolyT.map (f := fun a i => f a i * root a ^ (A a * min i (A a) * min (k a) (A a))) ht (PolyT.range hA)
  exact PolyT.mod (PolyT.sum hl) hM

theorem numericHat_eq (A M root k : ℕ) [NeZero M] (hk : k ≤ A)
    (f : ℕ → ℕ) (z : ℕ → ZMod M) (hf : ∀ i < A, (f i : ZMod M) = z i) :
    numericHat A M root k f = (hatOut A (root : ZMod M) z k).val := by
  unfold numericHat
  rw [← ZMod.val_natCast]
  congr 1
  rw [← List.sum_toFinset (fun i => f i * root ^ (A * min i A * min k A)) List.nodup_range,
    List.toFinset_range]
  push_cast
  unfold hatOut
  apply Finset.sum_congr rfl
  intro i hi
  have hi' := Finset.mem_range.mp hi
  simp [Nat.min_eq_left (Nat.le_of_lt hi'), Nat.min_eq_left hk, hf i hi']

def numericFa (A M root r k : ℕ) (f : ℕ → ℕ) : ℕ :=
  (root ^ (min r A * min k A) * numericHat A M root k f) % M

theorem numericFa_polyT {α : Type} [Enc α] {A : α → ℕ} {M root r k : α → ℕ} {f : α → ℕ → ℕ}
    (hA : Unary A) (hM : PolyT M) (hroot : PolyT root) (hr : PolyT r) (hk : PolyT k)
    (hf : PolyT fun p : α × ℕ => f p.1 p.2) :
    PolyT fun a => numericFa (A a) (M a) (root a) (r a) (k a) (f a) := by
  have hru : Unary fun a => min (r a) (A a) := Unary.of_le hA (PolyT.min' hr hA.poly) (fun _ => Nat.min_le_right _ _)
  have hku : Unary fun a => min (k a) (A a) := Unary.of_le hA (PolyT.min' hk hA.poly) (fun _ => Nat.min_le_right _ _)
  have hp := PolyT.pow_unary hroot (Unary.mul hru hku)
  have hh : PolyT fun a => numericHat (A a) (M a) (root a) (k a) (f a) :=
    numericHat_polyT (A := A) (M := M) (root := root) (k := k) (f := f) hA hM hroot hk hf
  exact PolyT.mod (PolyT.mul hp hh) hM

theorem numericFa_eq (A M root r k : ℕ) [NeZero M] (hr : r ≤ A) (hk : k ≤ A)
    (f : ℕ → ℕ) (z : ℕ → ZMod M) (hf : ∀ i < A, (f i : ZMod M) = z i) :
    numericFa A M root r k f = (faOut A (root : ZMod M) r z k).val := by
  unfold numericFa faOut
  rw [Nat.min_eq_left hr, Nat.min_eq_left hk, numericHat_eq A M root k hk f z hf, ZMod.val_mul]
  rw [← Nat.cast_pow, ZMod.val_natCast]
  simp [Nat.mul_mod, hatOut]

def numericMid (A M root inv k j : ℕ) (f g : ℕ → ℕ) : ℕ :=
  numericFa A M inv k j (fun t => (numericHat A M root t f * numericHat A M root t g) % M)

theorem numericMid_polyT {α : Type} [Enc α] {A : α → ℕ} {M root inv k j : α → ℕ} {f g : α → ℕ → ℕ}
    (hA : Unary A) (hM : PolyT M) (hroot : PolyT root) (hinv : PolyT inv) (hk : PolyT k) (hj : PolyT j)
    (hf : PolyT fun p : α × ℕ => f p.1 p.2) (hg : PolyT fun p : α × ℕ => g p.1 p.2) :
    PolyT fun a => numericMid (A a) (M a) (root a) (inv a) (k a) (j a) (f a) (g a) := by
  let β := α × ℕ
  have hA' : Unary fun p : β => A p.1 := Unary.comp hA PolyT.id'.fst
  have hM' : PolyT fun p : β => M p.1 := PolyT.comp hM PolyT.id'.fst
  have hroot' : PolyT fun p : β => root p.1 := PolyT.comp hroot PolyT.id'.fst
  have hf' : PolyT fun p : β × ℕ => f p.1.1 p.2 :=
    PolyT.comp (f := fun q : α × ℕ => f q.1 q.2)
      (g := fun p : β × ℕ => (p.1.1, p.2)) hf (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hg' : PolyT fun p : β × ℕ => g p.1.1 p.2 :=
    PolyT.comp (f := fun q : α × ℕ => g q.1 q.2)
      (g := fun p : β × ℕ => (p.1.1, p.2)) hg (PolyT.prodMk PolyT.id'.fst.fst PolyT.id'.snd)
  have hhatf : PolyT fun p : β => numericHat (A p.1) (M p.1) (root p.1) p.2 (f p.1) :=
    numericHat_polyT (A := fun p : β => A p.1) (M := fun p : β => M p.1)
      (root := fun p : β => root p.1) (k := fun p : β => p.2) (f := fun p i => f p.1 i)
      hA' hM' hroot' PolyT.id'.snd hf'
  have hhatg : PolyT fun p : β => numericHat (A p.1) (M p.1) (root p.1) p.2 (g p.1) :=
    numericHat_polyT (A := fun p : β => A p.1) (M := fun p : β => M p.1)
      (root := fun p : β => root p.1) (k := fun p : β => p.2) (f := fun p i => g p.1 i)
      hA' hM' hroot' PolyT.id'.snd hg'
  have hmid : PolyT fun p : α × ℕ => (numericHat (A p.1) (M p.1) (root p.1) p.2 (f p.1) *
      numericHat (A p.1) (M p.1) (root p.1) p.2 (g p.1)) % M p.1 :=
    PolyT.mod (PolyT.mul hhatf hhatg) hM'
  exact numericFa_polyT (A := A) (M := M) (root := inv) (r := k) (k := j)
    (f := fun a t => (numericHat (A a) (M a) (root a) t (f a) * numericHat (A a) (M a) (root a) t (g a)) % M a)
    hA hM hinv hk hj hmid

theorem numericMid_eq (A M root inv k j : ℕ) [NeZero M] (hk : k ≤ A) (hj : j ≤ A)
    (f g : ℕ → ℕ) (z w : ℕ → ZMod M)
    (hf : ∀ i < A, (f i : ZMod M) = z i) (hg : ∀ i < A, (g i : ZMod M) = w i) :
    numericMid A M root inv k j f g = (midOut A (root : ZMod M) (inv : ZMod M) k z w j).val := by
  unfold numericMid midOut
  apply numericFa_eq A M inv k j hk hj
  intro t ht
  rw [numericHat_eq A M root t (Nat.le_of_lt ht) f z hf,
    numericHat_eq A M root t (Nat.le_of_lt ht) g w hg]
  simp [ZMod.natCast_mod, ZMod.natCast_zmod_val]

def numericIb (A M inv tinv j : ℕ) (f : ℕ → ℕ) : ℕ :=
  (tinv * numericHat A M inv j f) % M

theorem numericIb_polyT {α : Type} [Enc α] {A : α → ℕ} {M inv tinv j : α → ℕ} {f : α → ℕ → ℕ}
    (hA : Unary A) (hM : PolyT M) (hinv : PolyT inv) (htinv : PolyT tinv) (hj : PolyT j)
    (hf : PolyT fun p : α × ℕ => f p.1 p.2) :
    PolyT fun a => numericIb (A a) (M a) (inv a) (tinv a) (j a) (f a) := by
  have hh := numericHat_polyT (A := A) (M := M) (root := inv) (k := j) (f := f) hA hM hinv hj hf
  exact PolyT.mod (PolyT.mul htinv hh) hM

theorem numericIb_eq (A M inv tinv j : ℕ) [NeZero M] (hj : j ≤ A)
    (f : ℕ → ℕ) (z : ℕ → ZMod M) (hf : ∀ i < A, (f i : ZMod M) = z i) :
    numericIb A M inv tinv j f = (ibOut A (inv : ZMod M) (tinv : ZMod M) z j).val := by
  unfold numericIb ibOut
  rw [numericHat_eq A M inv j hj f z hf, ZMod.val_mul, ZMod.val_natCast]
  simp [Nat.mul_mod]

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.numericHat_polyT
#print axioms IntegerFactoring.HalfAppendix.numericHat_eq
#print axioms IntegerFactoring.HalfAppendix.numericFa_polyT
#print axioms IntegerFactoring.HalfAppendix.numericFa_eq
#print axioms IntegerFactoring.HalfAppendix.numericMid_polyT
#print axioms IntegerFactoring.HalfAppendix.numericMid_eq
#print axioms IntegerFactoring.HalfAppendix.numericIb_polyT
#print axioms IntegerFactoring.HalfAppendix.numericIb_eq
