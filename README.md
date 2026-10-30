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
every field. The stronger [geometric cut budgets](research/geometric-cut-rank-budget.md)
now exclude the finite rank-sum amplification criterion through **95 labels**.
[Explicit integer Jacobian certificates](research/geometric-subspace-compression.md)
force common input and output subspaces for arbitrary matrix tuples.
They strengthen the earlier flag quotients, and the weighted budgets are
conditioned on an actual piece of largest imbalance. Exact rational
verification covers every required class. These restrictions concern that
criterion and do not change the factoring exponent. Their replays are
included in `make verify-frontier`.

The [two-product complement theorem](research/two-product-complement.md)
proves `rank(E_d - P tensor Q - R tensor S) >= d*(d-1)/2` over every
field, including singular matrix pencils. A separate
[scalar-relaxation example](research/weighted-cut-relaxation-limit.md)
shows why these rank budgets alone cannot settle every dimension:
abstract rank lists attain `O(d log d)` cost while satisfying them.
Realizing those lists by compatible tensor pieces remains unresolved.
The [geometric extension](research/geometric-cut-relaxation-limit.md)
proves that the same lists also satisfy every new geometric subset budget.
The same theorem gives a sharp common-block-flag bound, ruling out
nested projection realizations of the dyadic rank lists and quantifying
the rank cost of a perturbation.
The [formal-degeneration lemma](research/weighted-cut-border-stability.md)
also preserves the weighted quotient budgets when individual pieces
have Laurent-series poles. Moving geometric subspaces and a
[factor-kernel refinement](research/rank-one-formal-cutoff.md)
now exclude formal rank-sum seeds through **71 labels**. The exact
95-label cutoff uses stronger small-product arguments and remains a
separate theorem.

The [quadratic Boolean-code oracle](research/quadratic-code-oracle.md)
computes crossing coefficients in polynomial time without expanding
the evaluation words. Its `2^O(sqrt(lambda))` coefficient count and
effective factor construction meet a broader transfer interface, but
the [word-restriction theorem](research/word-restriction-barrier.md)
now excludes a contraction from this fixed 23-term coefficient template.
For arbitrary distinct binary words of subexponential length, one piece
in every coefficient-span decomposition has both widths
`2*lambda-o(lambda)`, regardless of the number of pieces. The obstruction
also covers coordinate weights and suitable unions of label embeddings.
The oracle and its quadratic-annihilator theorem remain valid.

The [dual code-section normal forms](research/code-section-normal-forms.md)
reduce arbitrary ambient code lengths to finite affine-selector data.
Complete screens cover all 67,108,864 coordinate-section classes at
logical width five and all 65,536 general-section classes at width three.
Neither screen improves the elementary pencil's maximum rank product;
the subsequent [frame-filtration theorem](research/frame-filtration-barrier.md)
now proves the sharp leading obstruction `rho >= H_2(1/4)` for all
affine-selector frames and arbitrary mixtures of their coefficient spans.
Its common projection preserves the crossing and never increases either
cut rank. A separate [diagonal-word theorem](research/diagonal-word-rank-sum-barrier.md)
covers every unsymmetrized tensor-word combination in characteristic two:
any crossing decomposition in that space has rank sum at least `d^(3/2)`,
regardless of its number of pieces. Off-diagonal nonlinear selectors and
nonlinear label maps require separate analysis; the full conservation
result below now closes the former whenever they retain that identity.

The broader [Boolean-incidence bound](research/incidence-kernel-rank-sum-barrier.md)
forces the same `d^(3/2)` cost over every field for all conservation
kernels supported on `w subset S`, including arbitrary row controls.
Every nonzero unrestricted conservation kernel still costs at least `d`,
so a passing 128-label seed in that larger space has at most seven pieces.
In characteristic two, the new
[full conservation-space theorem](research/full-conservation-rank-sum-barrier.md)
now closes that entire larger space at every size: arbitrary kernels,
including off-incidence and nilpotent corrections, satisfy
`sum sqrt(r_i*c_i)>=d^(3/2)`. The proof uses matched row and column
bases from a double Cauchy–Binet expansion. It covers the full
four-feature finite-state tensor space and supplied Laurent limits.
Constructions using incompatible physical label bases remain outside it.
A [fixed-target Feistel quotient](research/fixed-target-feistel-barrier.md)
also excludes arbitrary controlled-translation coefficient spans with
one common target set. Two nonlinear rounds with different target sets
escape this quotient: a checked sixteen-label tensor has ranks
`(254,18)`, while the proposed quotient has ranks `(256,1)`.

An [exact exterior-algebra construction](research/bi-exterior-crossing-algebra.md)
now realizes both cut ranks as multiplication ranks in one square-zero
algebra. This permits arbitrary position-dependent and mixed-degree
polynomials, with an explicit crossing endpoint and exact rank invariance
under all invertible linear changes of the generators. One sixteen-label
direction has actual ranks `(36,36)`, compared with `(46,46)` for the
earlier Boolean-selector realization. It is not a crossing decomposition.
The accompanying [degree and quadratic-rank budgets](research/bi-exterior-degree-budget.md)
exclude every seed in this entire model through **256 labels**, including
supplied formal limits. A binary seed in this model must therefore have
at least 512 labels to pass the rank-sum criterion.

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

The [raw rank-sum extension](research/rank-sum-global-transfer.md)
also permits exponentially many pieces when their total displayed
cost satisfies `sum sqrt(r_i*c_i)<=2^(2*rho*t+O(t/log(t+2)))`.
Keeping exact channel cardinalities until after type amplification
preserves the same quantitative time and space endpoint.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
