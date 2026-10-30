/-
SPDX-License-Identifier: Apache-2.0

Polynomial-time access to an Appendix B level given its unary index and an
input-width bound. This does not yet certify the selected-level search.
-/
import LevelSelect
import GeneralPieces.Num.PartA
import ToTuring.Lib.Nat2
import ToTuring.Lib.List2

open Groundwork Levels.Gen

namespace IntegerFactoring.HalfAppendix

private theorem level_x_closed (x l : ℕ) :
    ClauseBound.Gen.x x l + 6 * l + 2 * x + 6 = (3 * x + 6) * 2 ^ l := by
  induction l with
  | zero => simp [ClauseBound.Gen.x]; omega
  | succ l ih =>
      rw [GeneralPieces.PartA.x_succ, pow_succ]
      nlinarith

theorem level_parameter_polyT (x : ℕ) {α : Type} [Enc α] {l : α → ℕ} (hl : Unary l) :
    PolyT fun a => ClauseBound.Gen.x x (l a) := by
  have hp : PolyT fun a => 2 ^ l a := PolyT.two_pow hl
  have hlp : PolyT l := Unary.poly hl
  have hc : PolyT fun a => (3 * x + 6) * 2 ^ l a - (6 * l a + 2 * x + 6) := by fun_prop
  apply hc.congr
  intro a
  have h := level_x_closed x (l a)
  omega

theorem level_parameter_unary (x : ℕ) {α : Type} [Enc α] {l n : α → ℕ}
    (hl : Unary l) (hn : Unary n) (hbound : ∀ a, 2 ^ l a ≤ n a) :
    Unary fun a => ClauseBound.Gen.x x (l a) := by
  have hlim : Unary fun a => (3 * x + 6) * n a := by fun_prop
  apply Unary.of_le hlim (level_parameter_polyT x hl)
  intro a
  have h := level_x_closed x (l a)
  have h' := Nat.mul_le_mul_left (3 * x + 6) (hbound a)
  omega

theorem level_K_polyT (x : ℕ) {α : Type} [Enc α] {l n : α → ℕ}
    (hl : Unary l) (hn : Unary n) (hbound : ∀ a, 2 ^ l a ≤ n a) :
    PolyT fun a => K x (l a) := by
  have h := PolyT.two_pow (level_parameter_unary x hl hn hbound)
  exact h.congr fun a => (ClauseBound.Gen.K_eq x (l a)).symm

def greatestFold (p : ℕ → Bool) (n : ℕ) : ℕ :=
  (List.range n).foldl (fun s k => if p (k + 1) then k + 1 else s) 0

theorem greatestFold_eq (p : ℕ → Bool) : ∀ n,
    greatestFold p n = Nat.findGreatest (fun k => p k = true) n
  | 0 => rfl
  | n + 1 => by
      have ih := greatestFold_eq p n
      unfold greatestFold at ih ⊢
      rw [List.range_succ, List.foldl_append]
      simp only [List.foldl_cons, List.foldl_nil, Nat.findGreatest_succ]
      rw [ih]

theorem greatestFold_bound (p : ℕ → Bool) : ∀ (ks : List ℕ) (s m : ℕ),
    s ≤ m → (∀ k ∈ ks, k + 1 ≤ m) →
      ks.foldl (fun s k => if p (k + 1) then k + 1 else s) s ≤ m
  | [], s, m, hs, _ => hs
  | k :: ks, s, m, hs, hk => by
      apply greatestFold_bound p ks _ m ?_ (by intro j hj; exact hk j (by simp [hj]))
      dsimp
      split <;> simp_all

theorem findGreatest_polyT {α : Type} [Enc α] {m : α → ℕ} {p : α → ℕ → Bool}
    (hm : Unary m) (hp : PolyT fun q : α × ℕ => p q.1 q.2) :
    PolyT fun a => Nat.findGreatest (fun k => p a k = true) (m a) := by
  have hlist : PolyT fun a => List.range (m a) := PolyT.range hm
  have hstep : PolyT fun q : α × ℕ × ℕ =>
      if p q.1 (q.2.2 + 1) then q.2.2 + 1 else q.2.1 := by fun_prop
  obtain ⟨_, t, ht, _, _, hlen⟩ := hm.elim
  have hf : PolyT fun a => (List.range (m a)).foldl
      (fun s k => if p a (k + 1) then k + 1 else s) 0 := by
    refine PolyT.foldl hlist (PolyT.const 0) hstep ?_
    refine ⟨t, ht, ?_⟩
    intro a k
    let ks := (List.range (m a)).take k
    have hb := greatestFold_bound (p a) ks 0 (m a) (Nat.zero_le _) (by
      intro j hj
      have : j ∈ List.range (m a) := List.mem_of_mem_take hj
      exact List.mem_range.mp this)
    have hs := Nat.size_le_size hb
    have hsize : Nat.size (m a) ≤ m a := Nat.size_le.2 Nat.lt_two_pow_self
    have hlen' := hlen a
    rw [len_units, List.length_replicate] at hlen'
    change len (ks.foldl _ _) ≤ t (len a)
    rw [len_nat]
    omega
  exact hf.congr fun a => greatestFold_eq (p a) (m a)

private theorem findGreatest_bounded_congr (P Q : ℕ → Prop) [DecidablePred P] [DecidablePred Q] :
    ∀ n, (∀ k, k ≤ n → (P k ↔ Q k)) → Nat.findGreatest P n = Nat.findGreatest Q n
  | 0, _ => rfl
  | n + 1, h => by
      rw [Nat.findGreatest_succ, Nat.findGreatest_succ]
      have h' := h (n + 1) (le_refl _)
      by_cases hp : P (n + 1)
      · simp [hp, h'.mp hp]
      · have hq : ¬Q (n + 1) := fun hq => hp (h'.mpr hq)
        simp only [hp, hq, ↓reduceIte]
        exact findGreatest_bounded_congr P Q n (fun k hk => h k (by omega))

theorem selectedLevel_polyT (x : ℕ) {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    PolyT fun a => selectedLevel x (n a) := by
  have hnp : PolyT n := Unary.poly hn
  let m : α → ℕ := fun a => Nat.log2 (n a) + 1
  have hmp : PolyT m := by dsimp [m]; fun_prop
  have hm : Unary m := by
    have hlim : Unary fun a => n a + 1 := by fun_prop
    apply Unary.of_le hlim hmp
    intro a
    dsimp [m]
    rw [Nat.log2_eq_log_two]
    have h := Nat.log_le_self 2 (n a)
    omega
  let cap : α × ℕ → ℕ := fun p => min p.2 (m p.1)
  have hcap : PolyT cap := by dsimp [cap]; fun_prop
  have hcapu : Unary cap := by
    have hlim : Unary fun p : α × ℕ => m p.1 := Unary.comp hm (by fun_prop)
    exact Unary.of_le hlim hcap (fun p => Nat.min_le_right _ _)
  have hn2 : Unary fun p : α × ℕ => 2 * (n p.1 + 1) := by fun_prop
  have hpow : ∀ p : α × ℕ, 2 ^ cap p ≤ 2 * (n p.1 + 1) := by
    intro p
    calc
      2 ^ cap p ≤ 2 ^ m p.1 := Nat.pow_le_pow_right (by norm_num) (Nat.min_le_right _ _)
      _ = 2 * 2 ^ Nat.log 2 (n p.1) := by dsimp [m]; rw [Nat.log2_eq_log_two, pow_succ]; ring
      _ ≤ 2 * (n p.1 + 1) := Nat.mul_le_mul_left _ (Nat.pow_log_le_add_one 2 _)
  have hK : PolyT fun p : α × ℕ => K x (cap p) := level_K_polyT x hcapu hn2 hpow
  have hp : PolyT fun p : α × ℕ => decide ((K x (cap p)) ^ 40 ≤ n p.1) := by fun_prop
  have hfind := findGreatest_polyT (p := fun a k => decide ((K x (cap (a, k))) ^ 40 ≤ n a)) hm hp
  apply hfind.congr
  intro a
  unfold selectedLevel
  rw [← Nat.log2_eq_log_two]
  apply findGreatest_bounded_congr
  intro k hk
  dsimp [cap, m] at *
  simp [Nat.min_eq_left hk]

theorem selectedLevel_unary (x : ℕ) {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    Unary fun a => selectedLevel x (n a) := by
  have hlim : Unary fun a => n a + 1 := by fun_prop
  apply Unary.of_le hlim (selectedLevel_polyT x hn)
  intro a
  have h := Nat.findGreatest_le (P := fun l => K x l ^ 40 ≤ n a) (Nat.log 2 (n a) + 1)
  have hlog := Nat.log_le_self 2 (n a)
  change Nat.findGreatest _ _ ≤ _
  omega

theorem selected_K_polyT (x : ℕ) {α : Type} [Enc α] {n : α → ℕ} (hn : Unary n) :
    PolyT fun a => K x (selectedLevel x (n a)) := by
  have hlim : Unary fun a => 2 * (n a + 1) := by fun_prop
  apply level_K_polyT x (selectedLevel_unary x hn) hlim
  intro a
  have h := Nat.findGreatest_le (P := fun l => K x l ^ 40 ≤ n a) (Nat.log 2 (n a) + 1)
  have hp : 2 ^ selectedLevel x (n a) ≤ 2 ^ (Nat.log 2 (n a) + 1) :=
    Nat.pow_le_pow_right (by norm_num) h
  calc
    2 ^ selectedLevel x (n a) ≤ 2 ^ (Nat.log 2 (n a) + 1) := hp
    _ = 2 * 2 ^ Nat.log 2 (n a) := by rw [pow_succ]; ring
    _ ≤ 2 * (n a + 1) := Nat.mul_le_mul_left _ (Nat.pow_log_le_add_one 2 _)

theorem selected_K_unary (x : ℕ) (hx : 1 ≤ x) {α : Type} [Enc α] {n : α → ℕ}
    (hn : Unary n) : Unary fun a => K x (selectedLevel x (n a)) := by
  have hlim : Unary fun a => n a + K x 0 := by fun_prop
  apply Unary.of_le hlim (selected_K_polyT x hn)
  intro a
  by_cases hlarge : K x 0 ^ 40 ≤ n a
  · have hs := (selectedLevel_spec x (n a) hx hlarge).1
    have hk : K x (selectedLevel x (n a)) ≤ K x (selectedLevel x (n a)) ^ 40 :=
      Nat.le_self_pow (by norm_num) _
    omega
  · have hmono : ∀ l, K x 0 ≤ K x l := by
      intro l
      induction l with
      | zero => exact le_rfl
      | succ l ih => exact ih.trans (GeneralPieces.PartA.K_lt_le x l hx).1.le
    have hzero : selectedLevel x (n a) = 0 := by
      unfold selectedLevel
      apply Nat.findGreatest_eq_zero_iff.mpr
      intro k _ _ hk
      have h := Nat.pow_le_pow_left (hmono k) 40
      omega
    rw [hzero]
    omega

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.level_K_polyT
#print axioms IntegerFactoring.HalfAppendix.selectedLevel_polyT
#print axioms IntegerFactoring.HalfAppendix.selected_K_polyT
#print axioms IntegerFactoring.HalfAppendix.selected_K_unary
