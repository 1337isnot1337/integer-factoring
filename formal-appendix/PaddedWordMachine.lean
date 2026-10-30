/-
SPDX-License-Identifier: Apache-2.0

Conditional finite-alphabet machine bridge for a polynomial-time padded
bit-string algorithm. The factoring-specific padded algorithm is not supplied
by this file.
-/
import PaddedDescent.Pad
import CountKSat.Host
import Challenge

open Groundwork Turing PaddedDescent

namespace IntegerFactoring.HalfAppendix

def wordPadM (e : List Bool → ℕ) (A : List Bool × List Unit → List Bool) :
    Stepped (List Bool) (List Bool) where
  σ := (List Bool × List Unit) ⊕ List Bool
  init x := .inl (x, [()])
  step s := s.elim (fun p => bif decide (e p.1 ≤ Nat.log2 p.2.length) then .inr (A p)
    else .inl (p.1, p.2 ++ p.2)) fun b => .inr b
  out s := s.elim (fun _ => none) fun b => some b

@[instance_reducible] def wordPadM_isPolyT {e : List Bool → ℕ} {A : List Bool × List Unit → List Bool}
    (he : PolyT e) (hA : PolyT A) : IsPolyT (wordPadM e A) where
  enc := (inferInstance : Enc ((List Bool × List Unit) ⊕ List Bool))
  carried := (inferInstance : Carried ((List Bool × List Unit) ⊕ List Bool))
  init := (by fun_prop : PolyT fun x : List Bool => (.inl (x, [()]) : (List Bool × List Unit) ⊕ List Bool))
  step := (by fun_prop : PolyT fun s : (List Bool × List Unit) ⊕ List Bool =>
    s.elim (fun p => bif decide (e p.1 ≤ Nat.log2 p.2.length) then
      (.inr (A p) : (List Bool × List Unit) ⊕ List Bool) else .inl (p.1, p.2 ++ p.2)) fun b => .inr b)
  out := (by fun_prop : PolyT fun s : (List Bool × List Unit) ⊕ List Bool => s.elim (fun _ => none) fun b => some b)

private theorem size_bools (x : List Bool) : size x = 2 * x.length + 1 := by
  induction x with
  | nil => rfl
  | cons b x ih => simp only [size_list, List.map_cons, List.sum_cons, size_bool, List.length_cons] at ih ⊢; omega

private theorem size_pad (k : ℕ) : size (pad k) = 2 * 2 ^ k + 1 := by
  simp [pad, size_list, Nat.mul_comm]

private theorem pad_append (k : ℕ) : pad k ++ pad k = pad (k + 1) := by
  simp only [pad, List.replicate_append_replicate, pow_succ, Nat.mul_two]

private theorem length_pad (k : ℕ) : (pad k).length = 2 ^ k := List.length_replicate

theorem wordPadM_reaches (e : List Bool → ℕ) (A : List Bool × List Unit → List Bool)
    (x : List Bool) : ∀ (d k : ℕ), k + d = e x →
    (wordPadM e A).Reaches (2 * x.length + 2 * 2 ^ e x + size (A (x, pad (e x))) + 4)
      (A (x, pad (e x))) (.inl (x, pad k)) (d + 1)
  | 0, k, h => by
      obtain rfl : k = e x := by omega
      refine .step rfl ?_ ?_
      · change size x + size (pad (e x)) + 1 ≤ _
        rw [size_bools, size_pad]; omega
      · have hs : (wordPadM e A).step (.inl (x, pad (e x))) = .inr (A (x, pad (e x))) := by
          simp [wordPadM, length_pad, Nat.log2_two_pow]
        rw [hs]
        exact .halt rfl (by change size (A (x, pad (e x))) + 1 ≤ _; omega)
  | d + 1, k, h => by
      refine .step rfl ?_ ?_
      · change size x + size (pad k) + 1 ≤ _
        have hp : 2 ^ k ≤ 2 ^ e x := Nat.pow_le_pow_right (by norm_num) (by omega)
        rw [size_bools, size_pad]; omega
      · have hs : (wordPadM e A).step (.inl (x, pad k)) = .inl (x, pad (k + 1)) := by
          have hlt : ¬ e x ≤ k := by omega
          simp [wordPadM, length_pad, Nat.log2_two_pow, hlt, pad_append]
        rw [hs]
        exact wordPadM_reaches e A x d (k + 1) (by omega)

theorem wordPadM_runs (e : List Bool → ℕ) (A : List Bool × List Unit → List Bool) (x : List Bool) :
    (wordPadM e A).Runs x (e x + 1)
      (2 * x.length + 2 * 2 ^ e x + size (A (x, pad (e x))) + 4) (A (x, pad (e x))) :=
  wordPadM_reaches e A x (e x) 0 (by omega)

theorem exists_wordPad_machine {e : List Bool → ℕ} {A : List Bool × List Unit → List Bool}
    (he : PolyT e) (hA : PolyT A) :
    ∃ (M : TM2ComputableAux Bool Bool) (q : ℕ), (∀ k, Finite (M.tm.Γ k)) ∧
      ∀ x : List Bool, Nonempty (TM2OutputsInTime M.tm (x.map M.inputAlphabet.invFun)
        (some ((A (x, pad (e x))).map M.outputAlphabet.invFun))
        (q * (e x + 2) * (x.length + (2 * x.length + 2 * 2 ^ e x + size (A (x, pad (e x))) + 4) + 1) ^ q)) := by
  letI : IsPolyT (wordPadM e A) := wordPadM_isPolyT he hA
  obtain ⟨M, q, hfin, hrun⟩ := CountKSat.exists_tmW (wordPadM e A)
  refine ⟨M, q, hfin, ?_⟩
  intro x
  simpa [Nat.add_assoc] using hrun x (e x + 1)
    (2 * x.length + 2 * 2 ^ e x + size (A (x, pad (e x))) + 4)
    (A (x, pad (e x))) (wordPadM_runs e A x)

/-- Conditional assembly into the reference factoring-machine type. The
    hypotheses still require a concrete, correct polynomial-time algorithm. -/
theorem padded_factoring_machine {e : List Bool → ℕ} {A : List Bool × List Unit → List Bool}
    (E W : ℕ → ℕ) (he : PolyT e) (hA : PolyT A)
    (hE : ∀ N, e (Computability.encodeNat N) ≤ E (Computability.encodeNat N).length)
    (hW : ∀ N, size (A (Computability.encodeNat N, pad (e (Computability.encodeNat N)))) ≤
      W (Computability.encodeNat N).length)
    (hc : ∀ N, A (Computability.encodeNat N, pad (e (Computability.encodeNat N))) =
      Computability.encodeNat (Factoring.factor N)) :
    ∃ (M : Factoring.FactoringMachine) (q : ℕ), M.FiniteAlphabets ∧
      ∀ n, M.time n ≤ q * (E n + 2) * (3 * n + 2 * 2 ^ E n + W n + 5) ^ q := by
  obtain ⟨M, q, hfin, hrun⟩ := exists_wordPad_machine he hA
  let tm : Factoring.FactoringMachine :=
    { toTM2ComputableAux := M
      time := fun n => q * (E n + 2) * (3 * n + 2 * 2 ^ E n + W n + 5) ^ q
      outputsFun := by
        intro N
        let ev := Classical.choice (hrun (Computability.encodeNat N))
        rw [hc N] at ev
        refine ⟨ev.toEvalsTo, ev.steps_le_m.trans ?_⟩
        have hp : 2 ^ e (Computability.encodeNat N) ≤ 2 ^ E (Computability.encodeNat N).length :=
          Nat.pow_le_pow_right (by norm_num) (hE N)
        have heN := hE N
        have hwN := hW N
        rw [hc N] at hwN
        apply Nat.mul_le_mul (Nat.mul_le_mul_left q (by omega))
        apply Nat.pow_le_pow_left
        omega }
  exact ⟨tm, q, hfin, fun _ => le_rfl⟩

end IntegerFactoring.HalfAppendix

#print axioms IntegerFactoring.HalfAppendix.exists_wordPad_machine
#print axioms IntegerFactoring.HalfAppendix.padded_factoring_machine
