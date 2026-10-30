# Integer factoring bounds

[![Verify research certificates](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![Verify complete k6 theorem](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify-k6.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify-k6.yml)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-59636e)](LICENSE)

An open research effort on deterministic integer factoring, targeting the
[Big O Speedrun problem](https://bigospeedrun.com/integer-factoring). Given an
`n`-bit integer `N`, the task is to output its smallest prime factor. The time
bound counts bit operations for every input.

**Contributing research or asking an agent to help?** Start with
[autoresearch.md](autoresearch.md), [AGENTS.md](AGENTS.md), and
[CONTRIBUTING.md](CONTRIBUTING.md). Check the
[coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8)
and [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls)
for overlap; post your scope and handoff on the board.

Search the live board and open work alongside this checkout's research notes:

```sh
python3 tools/find_overlap.py grouping crossing
```

This read-only helper requires GitHub CLI access. Results are leads, not a
reservation or a novelty check; it reports which sources it did not search.

**Keep this README current.** Changes to results, evidence, verification, or
open obligations must update the relevant status, links, and reproduction
instructions here in the same change. Check the underlying evidence before
promoting a claim.

## Closed machine theorem: exponent 1/30

The degree-twelve variant has a closed Lean theorem asserting one deterministic,
finite-alphabet factoring machine with

```math
T(n)\le C\exp\!\left(\frac{Cn}{(\log n)^{1/30}}\right)\qquad(n\ge2).
```

It computes the smallest prime factor using the reference's original binary
input/output and machine-time definitions. The actual grouping program changes
from `G = 3K^8`, `K_next = 9K^16` to `G = 5K^6`, `K_next = 25K^12`.
The starting group count remains 48.

**Verified:** a fresh-source Lean 4.35.0 build completed all 2,252 jobs,
including the modified `Solution` and all eight new proof modules. A separate
[GitHub Actions rebuild](https://github.com/1337isnot1337/integer-factoring/actions/runs/41000040607)
also passed all 2,252 jobs. The three final axiom reports contain only
`propext`, `Classical.choice`, and `Quot.sound`.
The [verification record](evidence/k6/README.md) includes source hashes,
the complete compressed local build log, and the fresh endpoint reports.

This meets the official scoring criterion for `alpha = 1/30`, above the
reference's listed score 0. The same machine works for every positive scoring
tolerance; only the eventual input threshold changes. The
[rate and score proof](research/certified-k6.md) gives the exact statement.
No leaderboard submission or acceptance is claimed.

The exact comparison `log(12/11)/log(16) < 1/30` is also kernel checked.
This compares asymptotic upper-bound exponents, not measured runtimes or a
lower bound on the reference machine. The square-root-logarithm conjecture
remains open here.

[Closed theorem](variants/k6/FactoringK6/Certified.lean) ·
[Proof and score](research/certified-k6.md) ·
[Semantic and source audit](research/k6-proof-audit.md) ·
[Variant and reproduction](variants/k6/README.md) ·
[Fresh-source verifier](tools/check_k6.sh)

Run `make verify-k6` with Lean 4.35.0 and Lake on `PATH`. The verifier rebuilds
the modified reference program and the complete rate proof in a new source
tree, checks protected definitions and pinned dependencies, and prints the
endpoint axioms again. It rehashes the copied sources before issuing a success
receipt; a failed rerun invalidates old receipts. The separate
[k6 workflow](.github/workflows/verify-k6.yml) runs on proof, patch, verifier,
and lock-file pushes without duplicating pull-request builds.

## Stronger candidate: grouping and crossing compression

The proposed combination of a stronger grouping construction and a 56-channel
crossing factorization targets

```math
T(n)=\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right),
\qquad \alpha=\frac{\log(1/\rho)}{\log 4}=0.0755104681\ldots,
\quad \rho=\frac{8+\log_2 7}{12}.
```

Larger `alpha` is better. The comparison is between asymptotic bounds, not a
measured speedup. The target remains much weaker than the upstream conjecture
`alpha = 1/2`.

**Status: research candidate, not a completed machine-level theorem.** The
notes describe a uniform counting construction and its table-operation costs;
the generalized evaluator, complete machine-level proof, and Lean verification
of the improved factoring machine remain open. Exact finite certificates and
checked local lemmas establish only their stated scope.

[Integration and open obligations](research/integration.md) ·
[Recurrence and costs](research/recurrence.md) ·
[Grouping proof](research/grouping-improvement.md) ·
[Crossing construction](research/crossing-compression.md) ·
[Formal verification evidence](evidence/formal-verification.json)

Two further directions have separate scope:

- [Bounded-edge grouping](research/bounded-edge-grouping.md), based on the
  [dense-grouping argument](research/dense-grouping.md), gives smaller ordinary
  mathematical grouping bounds and exact finite checks. Its degree-two
  schedule still needs a uniform evaluator and a complete machine-time proof.
- [Appendix B integration](research/appendix-b-factoring-route.md) identifies
  a source-pinned route toward exponent `1/2`. The compatible dependencies,
  table adapters, top constructor, and factoring-machine composition remain open.

Neither direction is included in the verified `1/30` result.

## Baseline

The starting point is Isaac Rudich and Louis-Martin Rousseau's
[SubexponentialDeterministicFactoring](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring),
pinned to commit `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
Its published Lean theorem establishes one deterministic machine with time
`2^{o(n)} = N^{o(1)}`. Its separate square-root-logarithm conjecture is unproved.

## Research results

The repository contains ordinary mathematical proofs, exact finite checks,
Lean-checked local improvements, a quantitative bound for the unchanged reference
machine, and the new k6 machine theorem above. The k6 theorem has passed local
and hosted fresh-source builds. The stronger compressed candidates still lack
a complete machine theorem.

| Result | Value | Evidence and scope |
| --- | --- | --- |
| Degree-twelve modified factoring machine | `alpha = 1/30 = 0.033333...` | Closed machine and real-rate theorem; fresh-source local and hosted builds passed |
| Explicit analysis of the unchanged baseline | `alpha = log(12/11)/log(16) = 0.0313827205...` | Handwritten derivation from the pinned source's inequalities |
| Explicit time bound for the reference machine | `C * 2^(C * (explicitZu(n) + 1))` | Lean-checked bound for one finite-alphabet factoring machine; logarithmic-power conversion remains handwritten |
| Grouping with independent accuracy | `g = O(D^2 K^2 log K)` | Matching decomposition, Holder's inequality, and an explicit deterministic search |
| Formally checked grouping improvement | `g > K^2 + 4D^2K^4`; original accuracy works at `g = 5K^6` | Lean kernel checked, including network value preservation |
| Smaller crossing factorization | 56 integer channels instead of a trivial 64-channel factor | Explicit integer factors and exhaustive finite verification |
| Grouping-only rate target | `alpha = log(12/11)/log(4) = 0.0627654410...` | New parameter schedule and runtime analysis; modified machine not Lean verified |
| Combined rate target | `alpha = log(1/rho)/log(4) = 0.0755104681...`, `rho = (8+log_2(7))/12` | Combined mathematical construction; full formal integration remains open |

These bounds have very large constants; the finite checkers are not practical
factoring implementations or runtime benchmarks.

The combined target uses exact integer widths: `ceil(log_2(7^i))` is
`(7^i - 1).bit_length()`. A simpler rational width formula gives the slightly
weaker target `alpha = log(312/281)/log(4) = 0.0754879493...`. Neither formula
requires real-number approximation in the algorithm.

## Proofs and certificates

- [Certified k6 bound](research/certified-k6.md): the complete theorem, exact
  exponent comparison, official score implication, and machine-cost proof.
- [K6 audit](research/k6-proof-audit.md) and [build evidence](evidence/k6/README.md):
  preserved problem definitions, actual changed program, source checks, and
  completed verification results.
- [Baseline audit](research/baseline-audit.md): exact depth recurrence,
  contraction, all-level cost bound, and machine-cost bridge.
- [Grouping theorem](research/grouping-improvement.md): complete concentration
  proof, integer tests, deterministic construction, and integration obligations.
- [Crossing compression](research/crossing-compression.md): explicit factors,
  coefficient identities, changed multiplier, and a limitation of this coefficient
  family.
- [Exact crossing checker](tools/check_crossing.py): finite identities, channel
  encodings, and rank/width comparisons.
- [Combined integration](research/integration.md) and [recurrence proof](research/recurrence.md):
  explicit integer parameters, common widths, all-level costs, and remaining
  formalization obligations.
- [Finite grid step](research/grid-level-validation.md): actual row and column
  contractions, signed tables, composite moduli, and the required `8^t` multiplier.
- [Coefficient-basis obstruction](research/coefficient-basis.md): an exact minimax
  result ruling out further gains from linear recombinations in the same coefficient
  span, with an exact small-instance search.
- [GroupingThreshold.lean](formal/GroupingThreshold.lean) and
  [CrossingComplement.lean](formal/CrossingComplement.lean): checked finite lemmas.
- [ExplicitChain.lean](formal/ExplicitChain.lean): a checked bound on every
  chain exponent and a complete time bound for the unchanged reference machine.
  This formal envelope retains the source's `23/24` contraction; the sharper
  `11/12` baseline analysis remains handwritten. The k6 machine has the separate
  closed proof described above.
  [Verification evidence](evidence/formal-verification.json) records source hashes
  and the allowed axiom list.

## Reproduce the checks

The Python checks require only Python 3 and the standard library:

```sh
git clone https://github.com/1337isnot1337/integer-factoring.git
cd integer-factoring
make verify
```

Without Make, run the commands in the `verify` target of the [Makefile](Makefile).

With Lean 4.35.0 and Lake on the executable path, run:

```sh
make verify-formal
```

The script fetches the pinned public reference, rebuilds its `Solution`, and
checks every new formal file. It rejects axiom reports containing anything
outside `propext`, `Classical.choice`, and `Quot.sound`. The GitHub Actions workflow
specifies the same checks with a pinned Lean archive and checkout action. It runs
for pull requests, pushes to `main`, and manual dispatch; branch pushes without a
pull request do not start CI. No leaderboard submission has been made.

The [workflow badge](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
shows current GitHub Actions results. For the recorded local reference build and
axiom scope, see the [build record](evidence/upstream-build.json) and
[formal verification evidence](evidence/formal-verification.json). Historical
CI observations are kept in [evidence/ci-run.json](evidence/ci-run.json).

## Attribution

The factoring reduction and counting engine are due to the pinned upstream
development. Its README credits Isaac Rudich and Louis-Martin Rousseau and
describes its proof/code generation with Claude. Adapted source is identified
in file headers. This repository uses the upstream-compatible Apache-2.0 license.
