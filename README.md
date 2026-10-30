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
| Factoring from counts in every fixed positive characteristic | Ordinary conditional reduction; no score change | [Unitary-orbit recovery](research/unitary-orbit-factor-recovery.md) uses coprime tuples and perfect-power preprocessing, including composite characteristic |
| Factoring from quadratic parity counts | Ordinary conditional reduction; no score change | [Square roots of one](research/quadratic-congruence-factor-recovery.md) give a scalar witness and a fused quadratic-plus-linear arithmetic network |
| Fixed-precision modular counts | Ordinary conditional reduction and algebraic lift; no score change | [Prime-power lifting](research/fixed-precision-count-lifting.md) preserves the leading crossing exponent; [characteristic selection](research/high-order-characteristic-selection.md) avoids quadratic-count cancellation |
| Joined multiplication networks | Ordinary compiler theorem | [Weighted network compiler](research/weighted-network-top-compiler.md) preserves linear volume and exact witness counts |
| Local polynomial crossing certificates | Ordinary conditional transfer; no improving certificate yet | [Smith profiles](research/smith-profile-crossing-transfer.md), [auxiliary gauges](research/auxiliary-gauge-smith-transfer.md), [exact lattice enumeration](research/auxiliary-lattice-gauge-search.md), and [balanced physical gauges](research/balanced-physical-gauge-transfer.md) give exact finite construction criteria |
| Explicit auxiliary entanglers | Exact algebraic certificates; no score change | [Product-vector exclusion and its limits](research/explicit-entanglers-and-density-scalarization.md), including a complete base-field density screen |
| Exterior crossing boundaries | Ordinary reduction and scoped lower bounds | [Common projection and fixed-occupation bounds](research/coherent-exterior-boundary-reduction.md); a complete binary screen at four positions |

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
make verify-recovery      # Fixed-characteristic recovery and multiplication ports
make verify-constructive  # Entanglers, local Smith profiles, and auxiliary gauges
make verify-native        # Optional complete native screens (C++17/C++20)
make verify-docs          # Local Markdown-link check
```

The default checks use Python 3 and its standard library. Some optional checks use NumPy for speed. The optional `verify-native` target requires `g++`: the coherent-boundary replay uses C++17, and the complete product-density screen uses C++20. The latter examines all `2^30` normalized product densities for the displayed binary matrix.

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

The [unitary-orbit reduction](research/unitary-orbit-factor-recovery.md) is the shortest general recovery interface here. Pairwise coprime tuples have a free cyclic action, and their nontrivial orbit count is `k^(omega(N)-1)-1`. Taking the fixed tuple length equal to the coefficient characteristic, with exact perfect-power preprocessing, gives complete factoring with polynomial query overhead. Unique bounded inverse witnesses compile coprimality without changing the count. The earlier [prime-characteristic construction](research/fixed-characteristic-factor-recovery.md) and [parity power-divisor construction](research/parity-power-divisor-recovery.md) remain valid alternatives.

In characteristic two, [quadratic recovery](research/quadratic-congruence-factor-recovery.md) needs only roots of `z^2=1 mod N` in the canonical half interval. A recovered root gives a proper gcd factor. Its fused verifier uses one source-shaped network with pointwise map `hat(z)^2+hat(N)*hat(u)` and only a prefix restriction on `z`. The [weighted compiler](research/weighted-network-top-compiler.md) preserves linear total bit volume and charges construction costs. These reductions remove separate sparse-detector overhead, but do not supply a faster evaluator or higher unconditional exponent.

The [local Smith-profile criterion](research/smith-profile-crossing-transfer.md) converts two finite valuation spectra into an all-length rank bound for ordered polynomial matrix products. [Independent rational auxiliary gauges](research/auxiliary-gauge-smith-transfer.md) change only a constant boundary cost. [Balanced physical gauges](research/balanced-physical-gauge-transfer.md) add independent weights on all four external legs, paying only a polynomial loss in label-set size and explicitly correcting mean shifts. The elementary binary control recovers `rho=H_2(1/4)`; no screened certificate reaches `rho<1/sqrt(2)`. A [sum-preserving conjugation bound](research/sum-preserving-conjugation-smith-barrier.md) excludes another unweighted Smith family, without asserting a lower bound for its actual ranks.

The [explicit entangler certificates](research/explicit-entanglers-and-density-scalarization.md) separate product-vector exclusion from the stronger absence of a factorizing functional. The [exterior-boundary reduction](research/coherent-exterior-boundary-reduction.md) preserves the leading exponent under a common histogram projection and proves all-width restrictions for fixed occupations. Its remaining coherent family is open. An exact [131-bit counterexample](research/root-threshold-parity-counterexample.md) also shows why testing only whole root intervals cannot replace power-divisibility queries in parity recovery.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
