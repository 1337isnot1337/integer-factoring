# Integer factoring bounds

Research on deterministic integer factoring, targeting the [Big O Speedrun
integer-factoring problem](https://bigospeedrun.com/integer-factoring).
The input is an integer `N` written in `n` bits; the output is its smallest
prime factor. Running time counts bit operations for every input.

## Baseline

The starting point is Isaac Rudich and Louis-Martin Rousseau's
[SubexponentialDeterministicFactoring](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring),
pinned to commit `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
Its published Lean theorem establishes one deterministic machine with time
`2^{o(n)} = N^{o(1)}`. Its separate square-root-logarithm conjecture is unproved.

## Research results

The notes contain ordinary mathematical proofs and exact finite checks.
They do **not** constitute a Lean verification of an improved factoring machine.

| Result | Value | Evidence and scope |
| --- | --- | --- |
| Explicit analysis of the unchanged baseline | `alpha = log(12/11)/log(16) = 0.0313827205...` | Handwritten derivation from the pinned source's inequalities |
| Grouping with independent accuracy | `g = O(D^2 K^2 log K)` | Matching decomposition, Holder's inequality, and an explicit deterministic search |
| Smaller crossing factorization | 56 integer channels instead of a trivial 64-channel factor | Explicit integer factors and exhaustive finite verification |
| Grouping-only rate target | `alpha = log(12/11)/log(4) = 0.0627654410...` | New parameter schedule and runtime analysis; modified machine not Lean verified |
| Combined rate target | `alpha = log(1/rho)/log(4) = 0.0755104681...`, `rho = (8+log_2(7))/12` | Combined mathematical construction; full formal integration remains open |

Here `alpha` denotes the bound

\[
T(n)=\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right).
\]

Larger `alpha` is better. These are asymptotic bounds with very large constants;
the finite checkers are not practical factoring implementations or runtime
benchmarks. The combined target remains much weaker than the upstream conjecture
`alpha = 1/2`.

## Proofs and certificates

- [Baseline audit](research/baseline-audit.md): exact depth recurrence,
  contraction, all-level cost bound, and machine-cost bridge.
- [Grouping theorem](research/grouping-improvement.md): complete concentration
  proof, integer tests, deterministic construction, and integration obligations.
- [Crossing compression](research/crossing-compression.md): explicit factors,
  coefficient identities, changed multiplier, and a limitation of this coefficient
  family.
- [Exact crossing checker](tools/check_crossing.py): finite identities, channel
  encodings, and rank/width comparisons.

Further executable certificates, recurrence analysis, and checked Lean lemmas
are being added alongside these notes. No external leaderboard submission is
part of this repository update.

## Attribution

The factoring reduction and counting engine are due to the pinned upstream
development. Its README credits Isaac Rudich and Louis-Martin Rousseau and
describes its proof/code generation with Claude. Adapted source is identified
in file headers. This repository uses the upstream-compatible Apache-2.0 license.
