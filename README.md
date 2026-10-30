# Deterministic integer factoring bounds

[![Verification](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

This repository studies worst-case deterministic factoring. Given an integer `N >= 2` in binary, the machine outputs its smallest prime factor. Runtime is measured in finite-alphabet Turing-machine steps, with `n` equal to the input bit length.

## Main theorem: logarithmic saving exponent 1/2

**Lean proves the exact pinned `Factoring.SqrtLogLogConjecture`:** one deterministic finite-alphabet machine computes the smallest prime factor for every input in time

```text
T(n) <= C * 2^(C * n / (Nat.sqrt (Nat.log 2 n) + 1))
     = exp(O(n / sqrt(log n))).
```

The constant `C` does not depend on the input. This gives score **`alpha = 1/2`** under the [Big O Speedrun rule](https://bigospeedrun.com/integer-factoring). The bound is asymptotic. It is not a practical factoring speedup, and independent review remains pending.

The complete Lean endpoint is [`IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`](formal-appendix/HalfTheorem.lean). It proves the source statement `Factoring.SqrtLogLogConjecture` and depends only on `propext`, `Classical.choice`, and `Quot.sound`. The proof covers the factoring-table constructor, deterministic grouping, value preservation, padding, prefix recovery, and the finite-alphabet machine bound.

The composition uses the existing clique evaluator from [Anthropic's Appendix B](https://github.com/anthropics/formal-math) and the multiplication-network reduction by [Isaac Rudich and Louis-Martin Rousseau](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring). See the [written proof](research/half-rate-factoring.md), the [constructor](research/half-top-constructor.md), and the [evaluator audit](research/half-evaluator-audit.md).

## Results and evidence

| Claim | Status | Details |
| --- | --- | --- |
| Saving exponent `1/2` | **Lean-verified** | [Final theorem](formal-appendix/HalfTheorem.lean), [proof modules](formal-appendix/README.md) |
| Saving exponent `1/30` for the k6 machine | **Lean-verified** | [Proof](research/verified-k6.md), [machine theorem](variants/k6/FactoringK6/Certified.lean) |
| Sharper k6 saving `log(12/11)/log(12) > 1/29` | Ordinary proof, not Lean-verified at this rate | [Rate derivation](research/k6-sharp-rate.md) |
| Saving exponent above `1/2` | **Not established** | [Research index](research/INDEX.md) |
| Score transfer with arbitrary vanishing width error | Ordinary conditional theorem | [Proof and effective hypotheses](research/vanishing-slack-score-transfer.md) |
| Uniform endpoint from arbitrarily sparse finite seeds | Ordinary conditional theorem | [Bounded search and type amplification](research/gap-free-seed-score-transfer.md) |
| One uniform finite-seed endpoint over all coefficient fields | Ordinary conditional theorem | [Input-wise selection](research/universal-field-seed-selection.md), [exact rank specialization](research/finite-field-rank-specialization.md) |
| Auxiliary divisor-count exponent with leading coefficient one | Ordinary unconditional bound; overall score unchanged | [Exact integer parameter](research/optimized-divisor-count-bound.md) |
| All-field binary-word and middle-layer obstruction | Ordinary theorem; restricted family only | [Uniform rank floor `rho >= 3/4`](research/all-field-binary-word-barrier.md) |
| Matrix-entry-oracle rank and affine queries | Quasipolynomial time and polylogarithmic auxiliary space; local algorithm only | [Polynomial weighting and recomputed bases](research/polylog-space-rank-recovery.md), [scalar Berkowitz foundation](research/small-space-berkowitz-rank.md) |
| Border and finite-precision seed extraction | Ordinary conditional theorem; no stronger seed supplied | [Coefficient extraction and nonunit budgets](research/border-and-nonunit-seed-transfer.md) |

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
make verify-further       # New transfers, arithmetic bounds, and word certificates
make verify-energy-full   # Optional larger weighted-energy certificate
make verify-docs          # Local Markdown-link check
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

No uniform factoring exponent above `1/2` is proved here. One open approach seeks effective crossing pieces with total rank width `4*rho*lambda + o(lambda)` for `rho < 1/sqrt(2)`. The [score-transfer theorem](research/vanishing-slack-score-transfer.md) permits this arbitrary vanishing relative error, together with piece entropy `lambda/(log lambda)^(delta-o(1))`, and gives score `min(log2(1/rho), delta)`. Its uniform construction and factor-query hypotheses are part of the theorem. It supplies no new crossing family; finite rank improvements alone do not establish a stronger factoring bound.

The [finite-seed endpoint theorem](research/gap-free-seed-score-transfer.md) removes any gap or generation-speed requirement from an enumerable improving sequence. Its [input-wise field extension](research/universal-field-seed-selection.md) gives one machine attaining the best finite rank-sum score over all coefficient fields. The [specialization proof](research/finite-field-rank-specialization.md) preserves every finite rank profile over a finite field, and type amplification makes prime fields sufficient for the same infimum. Finding a better finite certificate remains the mathematical requirement.

Existence of one machine attaining the supremum of achievable scores also follows from [universal verified search](research/universal-verified-factoring-search.md), using polynomial-time primality testing. That is a standard search principle; the explicit tensor selector additionally describes the construction and its local resource costs.

The [all-field binary-word theorem](research/all-field-binary-word-barrier.md) rules out crossing pieces invariant under simultaneous permutation of all binary positions, and includes middle Hamming layers and products of such layers. One piece must have both ranks at least `2^(3t/2-O(log^2 t))`, uniformly in the field, so this family cannot reach the half-rate threshold. This does not restrict arbitrary crossing constructions or factoring algorithms.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

The same word certificate applies to [symmetric binary-energy shells](research/binary-energy-symmetry-barrier.md), giving `rho >= 3/4` for all fields. For the [specific reciprocal pencil](research/weighted-energy-reciprocal-pencil.md), a separate exact block argument proves the sharper optimum `rho = 1`. Other constructions on these low-growth spaces remain open. The [border-transfer theorem](research/border-and-nonunit-seed-transfer.md) makes polynomial or sufficiently precise truncated-parameter identities usable as exact finite seeds, with their coefficient-extraction costs included.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
