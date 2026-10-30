# Deterministic integer factoring bounds

[![Verification](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

This repository studies worst-case deterministic factoring. Given an integer `N >= 2` in binary, the machine outputs its smallest prime factor. Runtime is measured in finite-alphabet Turing-machine steps, with `n` equal to the input bit length.

**Want to contribute?** [Start here](docs/START_HERE.md) for setup, a reading path,
and small tasks. Reproductions, mathematical corrections, documentation, independent
review, and scoped negative results are welcome.

## Main theorem: logarithmic saving exponent 1/2

**Lean proves the exact pinned `Factoring.SqrtLogLogConjecture`:** one deterministic finite-alphabet machine computes the smallest prime factor for every input in time

```text
T(n) <= C * 2^(C * n / (Nat.sqrt (Nat.log 2 n) + 1))
     = exp(O(n / sqrt(log n))).
```

The constant `C` does not depend on the input. This gives score **`alpha = 1/2`** under the [Big O Speedrun rule](https://bigospeedrun.com/integer-factoring). The bound is asymptotic and does not provide a practical factoring speedup. Big O Speedrun independently rebuilt commit `e91442d3` and accepted the theorem with both Lean and nanoda; see the [reviewer's report](https://github.com/1337isnot1337/integer-factoring/issues/28#issuecomment-6533697488) and [successful verifier run](https://github.com/bigospeedrun/verifier/actions/runs/41108289491). Broader mathematical review is welcome.

The complete Lean endpoint is [`IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`](formal-appendix/HalfTheorem.lean). It proves the source statement `Factoring.SqrtLogLogConjecture` and depends only on `propext`, `Classical.choice`, and `Quot.sound`. The proof covers the factoring-table constructor, deterministic grouping, value preservation, padding, prefix recovery, and the finite-alphabet machine bound.

The composition uses the existing clique evaluator from [Anthropic's Appendix B](https://github.com/anthropics/formal-math) and the multiplication-network reduction by [Isaac Rudich and Louis-Martin Rousseau](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring). See the [written proof](research/half-rate-factoring.md), the [constructor](research/half-top-constructor.md), and the [evaluator audit](research/half-evaluator-audit.md).

## Results and evidence

| Claim | Status | Details |
| --- | --- | --- |
| Saving exponent `1/2` | **Lean-verified** | [Final theorem](formal-appendix/HalfTheorem.lean), [proof modules](formal-appendix/README.md) |
| Saving exponent `1/30` for the k6 machine | **Lean-verified** | [Proof](research/verified-k6.md), [machine theorem](variants/k6/FactoringK6/Certified.lean) |
| Sharper k6 saving `log(12/11)/log(12) > 1/29` | Ordinary proof, not Lean-verified at this rate | [Rate derivation](research/k6-sharp-rate.md) |
| Saving exponent above `1/2` | **Not established** | [Research index](research/INDEX.md) |

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
make verify-docs          # Local Markdown-link check
make verify-schur-fermion # Complete Schur entry formulas and exact comparisons
make verify-cancellation-modules # Exact module states and word-rate certificates
```

The Python checks use Python 3 and its standard library. Some optional checks use NumPy for speed.

For the Lean builds, install the stated Lean release, Lake, Git, and the prerequisites in the corresponding scripts. The scripts fetch pinned upstream revisions.

```sh
make verify-half-appendix # Lean 4.33.1; complete 1/2 theorem
make verify-k6           # Lean 4.35.0; complete 1/30 theorem
make verify-formal       # Lean 4.35.0; original local lemmas
```

The [GitHub Actions workflow](.github/workflows/verify.yml) runs these four groups on PRs and on `main`: exact Python checks, original Lean lemmas, k6, and the half-rate theorem. The half-rate script audits the axioms of its final theorem and supporting declarations.

## Repository map

| Directory | Contents |
| --- | --- |
| [`formal-appendix/`](formal-appendix/README.md) | Lean proof of the `1/2` theorem and its dependencies |
| [`variants/k6/`](variants/k6/FactoringK6/Certified.lean) | Earlier complete `1/30` proof and source patch |
| [`formal/`](formal/) | Supporting Lean lemmas and local checks |
| [`research/`](research/INDEX.md) | Mathematical proofs, candidate constructions, counterexamples, and open questions |
| [`evidence/`](evidence/) | Finite certificates and recorded verification outputs |
| [`tools/`](tools/) and [`tests/`](tests/) | Exact checkers, reproduction scripts, and unit tests |

## Current research frontier

No uniform factoring exponent above `1/2` is proved here. One open approach seeks effective crossing pieces with total rank width `4*rho*lambda + O(log lambda)` for `rho < 1/sqrt(2)`. A suitable family could improve the evaluator bound, but finite rank improvements alone do not establish that inequality.

The [complete Schur entry algorithm](research/schur-fermion-queries.md) computes
individual modes with determinants of order at most `3 sqrt(H)`, and full
crossing coefficients with determinants of order at most `4 sqrt(h)`. It works
directly in finite characteristic and removes factorial-precision conversion
from these queries. Its ordinary proof and independent exact checks improve
the available construction tools; no stronger rank bound is claimed.

[Finite module automata](research/cancellation-module-smith-automata.md)
retain cancellation between auxiliary paths. Their bounds improve monotonically
with precision, dominate the corresponding capped entry-order bounds, and
converge to the exact full-matrix word rate when both boundaries are free.
The checked crossing controls do not yet meet the strict improvement threshold.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with the [contributor quickstart](docs/START_HERE.md) and
[CONTRIBUTING.md](CONTRIBUTING.md). You can report a correction or reproduction
without proposing a new bound. [Open a focused issue](https://github.com/1337isnot1337/integer-factoring/issues/new/choose)
or choose an [available starter task](https://github.com/1337isnot1337/integer-factoring/issues?q=is%3Aissue%20is%3Aopen%20label%3A%22good%20first%20issue%22).

For substantial research, read [autoresearch.md](autoresearch.md) and coordinate
through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8).
Reviewers can use the [review guide](docs/REVIEWING.md) and
[exact-head review ledger](docs/REVIEW_LEDGER.md). Document source revisions,
proof status, reproduction commands, and remaining gaps. The
[contribution record](CONTRIBUTORS.md) retains credit for useful partial and
superseded work as well as incorporated results.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
