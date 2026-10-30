# Modular positivity with a scaled count

The candidate evaluator is intended to return `H * (q + b)` modulo each
small modulus, where `q` is the number of desired divisors, `b` is the known
contribution of the divisor 1, and `H` is a positive product of level
multipliers. A modulus may divide `H`, so division by `H` is unavailable.

[ModularPositivity.lean](../formal/ModularPositivity.lean) proves that, if
`H > 0`, `M > 0`, and `(2M+1) Hq < 4^M`, then `q > 0` exactly when a
modulus `2 <= m <= 2M` gives different residues for `H(q+b)` and `Hb`.
The proof uses the pinned upstream `eq_zero_of_all_dvd` theorem and cancels
the *added baseline* in a congruence; it never cancels or inverts `H`.

The same file proves that the smallness condition follows from
`M >= 4` and `Hq <= 2^(M-4)`. In particular, if `H <= 2^S` and
`q <= 2^n`, the choice `M = n + S + 4` is sufficient. This formalizes the
arithmetic interface used in [integration.md](integration.md); it does not
prove that the proposed evaluator computes these residues, that its concrete
multiplier obeys the assumed bound, or that its execution has the claimed
cost. Those remain separate obligations.

Checked against the pinned upstream commit in `upstream.lock.json` with Lean
4.35.0 using `lake env lean` on the file. The [axiom log](../evidence/lean-modular-positivity.log)
contains only `propext`, `Classical.choice`, and `Quot.sound`. This is a
kernel-checked interface lemma, not a verified factoring improvement.
