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

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
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

No uniform factoring exponent above `1/2` is proved here. One open approach seeks effective crossing pieces with total rank width `4*rho*lambda + O(log lambda)` for `rho < 1/sqrt(2)`. A suitable family could improve the evaluator bound, but finite rank improvements alone do not establish that inequality.

The [weighted cut-rank theorem](research/weighted-cut-rank-budget.md) proves
`sum_i sqrt(r_i*c_i) >= d*(2*H_(2d)-H_d-1)` for crossing decompositions over
every field. An exact strengthened certificate excludes the finite rank-sum
amplification criterion through 46 labels. The
[residual-budget refinement](research/weighted-cut-residual-budget.md)
raises the exclusion to 50 labels, and six
[positive rational duals](research/weighted-cut-dual-certificate.md)
extend it through 52 labels. These restrictions concern that criterion;
the [two-product refinement](research/two-product-seed-cutoff.md)
extends it further through **54 labels**. These restrictions leave other
constructions open and do not change the factoring exponent. The exact
replays are included in `make verify-frontier`.

The [two-product complement theorem](research/two-product-complement.md)
proves `rank(E_d - P tensor Q - R tensor S) >= d*(d-1)/2` over every
field, including singular matrix pencils. A separate
[scalar-relaxation example](research/weighted-cut-relaxation-limit.md)
shows why these rank budgets alone cannot settle every dimension:
abstract rank lists attain `O(d log d)` cost while satisfying them.
Realizing those lists by compatible tensor pieces remains unresolved.
The same theorem gives a sharp common-block-flag bound, ruling out
nested projection realizations of the dyadic rank lists and quantifying
the rank cost of a perturbation.
The [formal-degeneration lemma](research/weighted-cut-border-stability.md)
also preserves the weighted quotient budgets when individual pieces
have Laurent-series poles; the harmonic exclusion through 41 labels
survives that operation. The stronger finite cutoffs require separate
complement arguments and are not asserted for degenerations.

The [quadratic Boolean-code oracle](research/quadratic-code-oracle.md)
computes crossing coefficients in polynomial time without expanding
the evaluation words. Its `2^O(sqrt(lambda))` coefficient count and
effective factor construction meet a broader transfer interface, but
the required rank contraction is still open. The note also proves a
quadratic-annihilator theorem and reports the unsuccessful literal
32- and 64-label rank profiles.

The [dual code-section normal forms](research/code-section-normal-forms.md)
reduce arbitrary ambient code lengths to finite affine-selector data.
Complete screens cover all 67,108,864 coordinate-section classes at
logical width five and all 65,536 general-section classes at width three.
Neither screen improves the elementary pencil's maximum rank product;
these finite results do not exclude larger widths or other families.

The [global-block transfer](research/global-block-field-transfer.md)
permits coefficient fields to vary with the raw block length, with
explicit construction cost `2^O(t^h)` for any fixed `h>=1` and width
error `O(t/log(t+2))`. It retains the exact simultaneous time and space
bound `exp(O(n/log(n)^alpha))` and `O(n/log(n)^alpha)`, respectively,
conditional on the supplied contraction. A
[verified recovery extension](research/characteristic-independent-recovery.md)
covers finite fields of varying characteristic. These results broaden
the admissible construction interface without supplying a stronger
contraction.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
