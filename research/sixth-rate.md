# A checked rate bridge for the sixth-power schedule

## Scope and sources

[`formal/SixthRate.lean`](../formal/SixthRate.lean) proves the arithmetic
conversion needed by the sixth-power integration route. It does not construct
the candidate evaluator, bound that evaluator's costs, or assert that a faster
factoring machine exists. The ordinary comparison rate for the unchanged
reference is audited in [baseline-audit.md](baseline-audit.md).

The reference is Isaac Rudich and Louis-Martin Rousseau's
`SubexponentialDeterministicFactoring`, commit
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`, Apache-2.0. Its README
credits Claude for proof/code generation. This adaptation was prepared with
OpenAI Codex assistance. Mathlib is pinned by that checkout to
`56a005855ea17e0b0958f73ec801ae140b0a1679`, with Lean 4.35.0.

The new file imports Mathlib directly; it does not import or assume any
candidate factoring-machine theorem. Its statements expose every assumption
needed from a future implementation.

## Growth, cutoff, and volume

For any bottom count `K0 <= 1024`, define

\[
K_0=\texttt{K0},\qquad K_{d+1}=25K_d^{12}.
\]

This covers bottom counts 48, 600, and 1024 used by the original mathematical
chain and current integration work. The rate theorem does not identify these
different schedules. `groups_eq_of_recurrence` connects a separately defined
implementation with the same bottom count and recurrence.

The checked theorem `groups_scaled_le` gives

\[
2K_d\le 2^{11\cdot12^d}.
\]

Consequently, `log_le_of_cutoff` proves

\[
n<K_{d+1}^{40}
\quad\Longrightarrow\quad
\lfloor\log_2 n\rfloor+1\le5280\cdot12^d.
\tag{1}
\]

The **next-level** cutoff matters: an upper bound on `K_d` alone would not
give the needed lower bound on the selected depth. The runtime depth selector
must be proved to satisfy the cutoff for the actual schedule. A lower bound
on `K0` needed for the volume contraction is a separate obligation.

For natural numbers `A` and `v`, the theorem `volume_rate_of_cutoff` combines
(1) with a supplied volume bound to prove

\[
25^d v\le A23^d n
\quad\Longrightarrow\quad
v^{30}(\lfloor\log_2 n\rfloor+1)
\le5280A^{30}n^{30}.
\tag{2}
\]

The proof raises the volume estimate to the thirtieth power, multiplies by
the logarithm bound, and uses the exact integer certificate

\[
12\cdot23^{30}<25^{30}.
\]

It cancels the positive integer factor `25^(30*d)`. There is no floating-point
calculation or asymptotic estimate in this argument.

## Padding and the real-valued bound

Write `L = floor(log_2(n)) + 1`. The theorem `padded_rate` requires explicit
budgets for both a main term and a padding term:

\[
E\le v+p,\qquad v^{30}L\le Bm^{30},\qquad
p^{30}L\le Cm^{30}.
\]

It concludes

\[
E^{30}L\le2^{29}(B+C)m^{30}.
\tag{3}
\]

Here `m` can be `n+1`, so a positive bookkeeping cost at `n=0` is compatible
with the hypotheses. Source bounds with numerator `n` can be weakened to
numerator `n+1` by monotonicity. No padding, top-level cost, or compiler
overhead is dismissed as negligible without a bound.

`real_rate_of_power` proves, for an arbitrary numerator `m`,

\[
v^{30}L\le Bm^{30}
\quad\Longrightarrow\quad
v\le B^{1/30}\frac{m}{L^{1/30}}.
\tag{4}
\]

The denominator is positive even at `n=0`. With `m=n+1`,
`time_rate_of_power` carries (4) through any **separately proved** complete
time estimate `t <= C * 2^(C*(E+1))`. Thus both the root extraction and the
exponential time conversion are checked, conditional on their explicit
power-budget and time hypotheses.

The exact theorem `exponent_comparison` proves

\[
\frac{\log(12/11)}{\log16}
<\frac1{30}
<\frac{\log(25/23)}{\log12}.
\tag{5}
\]

The left comparison uses `12^30 < 16*11^30`; the right comparison uses the
contraction certificate above. This compares exponents. It is not a lower
bound on the reference machine's running time or a measured speedup.

## Remaining machine obligations

To use this bridge for a factoring theorem, supply the actual machine and
establish all of the following for its runtime-selected depth and costs:

1. Its group-count recurrence and bottom count agree with the schedule used
   in the cutoff theorem.
2. Its depth selector satisfies the next-level cutoff.
3. Its contracted main volume satisfies the hypothesis of (2), with all
   startup or exceptional levels accounted for in the fixed constant `A`.
4. Its full exponent is bounded by main volume plus padding, and the padding
   satisfies the second budget in (3), uniformly over all inputs.
5. Its correct, uniform factoring machine satisfies the complete time
   hypothesis, including reduction, parameter selection, table operations,
   moduli, compiler overhead, and finite-alphabet simulation.

The rate bridge does not supply any of these machine-specific facts by
assuming them. Independent review of the resulting full theorem is also
separate from kernel checking this file.

## Reproduction

Run `make verify` and `make verify-formal` from the repository root. The
formal target rebuilds the pinned reference and checks every file under
`formal/`, including this one. A focused check in the pinned reference's
Lake environment is:

```sh
cd .lake/upstream-reference
lake env lean ../../formal/SixthRate.lean
```

The source hash, exact checks run, axiom reports, and full-gate status are
recorded in [sixth-rate-verification.json](../evidence/sixth-rate-verification.json).
