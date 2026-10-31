# Deterministic integer factoring bounds

[![Verification](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

This repository studies worst-case deterministic factoring. Given an integer `N >= 2` in binary, the machine outputs its smallest prime factor. Runtime is measured in finite-alphabet Turing-machine steps, with `n` equal to the input bit length.

**Want to contribute?** [Start here](docs/START_HERE.md) for setup, a reading path,
and small tasks. Work directly on `main` by default and push useful checkpoints
frequently. Reproductions, mathematical corrections, documentation, independent
review, and scoped negative or partial research results are welcome.

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

Arithmetic research now uses [inverse sections](research/dyadic-inverse-sections.md)
and an [interval orbit-gcd interface](research/interval-orbit-gcd.md).
It supplies [exact complete-orbit aggregation](research/dyadic-orbit-aggregation.md)
and [global modular input circuits](research/interval-packet-circuits.md)
with `O(n+p)` arithmetic nodes and polynomial generation cost.
Recovering the genuine rational common factor within the desired time
bound remains open; these results do not improve the factoring score.
[Normalized integer-center jets](research/normalized-interval-jets.md)
also separate huge local powers of two from polynomial-cost unit data.
The [research index](research/INDEX.md) links proofs, exact evidence, and replay commands.

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
make verify-docs          # Local Markdown-link check
make verify-arithmetic    # All arithmetic-framework checkers
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

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work is committed directly to `main` by default; useful partial research can be integrated as a checkpoint without being promoted to a stronger theorem.

An [arithmetic framework program](research/arithmetic-framework.md) studies
deterministic class-group relation generation outside the crossing evaluator.
It proves a positive Fourier hitting lemma and deterministic extraction of
two-torsion from unsaturated relation data. A powering seed with published
expander amplification supplies biased words from generators without a group
decomposition. The fast generator and smooth-pair
construction steps remain open, so the score is still `1/2`. Reproduce the
current exact controls with `make verify-arithmetic-framework`.
The [actual quadratic-form controls](research/arithmetic-class-group-controls.md)
also recover proper factors for 18 finite targets and check complete
two-torsion recovery; their generator and order searches enumerate classes.
A [density/index transfer](research/arithmetic-density-lift.md) also turns
relations on a possibly proper small-prime subgroup into a global
class-number multiple. Fast pair extraction and deterministic coverage of
the missing ambiguous directions remain unresolved.
The [bounded smooth-norm interface](research/bounded-smooth-norm.md)
characterizes the required pairs and converts supplied representations at
`2^O(b/log b)` divisor-search cost; locating those representations efficiently
is still open.

The [Bell-diagonal seed construction](research/bell-diagonal-crossings.md)
gives exact two-cut transforms and excludes its complete symmetric
coefficient spans at local dimensions two and four. Dimension eight
remains open under that argument. Replay its physical tensors and pencil
ranks with `make verify-bell`; no stronger factoring exponent is claimed.
The [linear-pencil search](research/bell-linear-search.md) includes two
complete small-domain calibrations and explicit larger finite searches,
with no improvement in the recorded schedules.
The [polynomial synthesis engine](research/bell-polynomial-synthesis.md)
certifies uniform rational rank bounds and led to an explicit
[folded-cube crossing family](research/folded-bell-crossings.md).
Its optimal leading contraction within its complete coefficient span
is `H_2(1/4)=0.811278...`, above the critical threshold. The factoring
score is unchanged.
The [code-quotient analysis](research/code-quotient-crossings.md) preserves
the crossing duality, pays the maximal common polynomial division, and
identifies nonuniform finite-field cancellations as a remaining route;
ordinary covering-radius gains alone do not improve its pressure criterion.
The [signed-boundary extension](research/signed-boundary-crossings.md)
preserves duality under more general intertwiners and isotropic reductions,
with exact common-factor costs; its strict rank improvement remains open.
An [endpoint pressure theorem](research/bell-cubic-outer-rank-floor.md)
excludes the stated outer-rank-one cubics in every odd characteristic,
and quartics at local dimensions 8 and 16. Its first open power-of-two
quartic alphabet at least eight is 32.
For [self-reciprocal quartics](research/bell-quartic-rigidity.md), an exact
matrix identity relates a pressure witness to sparse Hadamard rank reduction;
it is a necessary condition and supplies no improved factoring bound.


## Contribution and attribution

Start with the [contributor quickstart](docs/START_HERE.md) and
[CONTRIBUTING.md](CONTRIBUTING.md). Work directly on `main` when you have
push access. Pull requests are mainly for external contributors or when
isolation and separate integration review are useful; every useful PR is
intended to merge, including partial and negative research. A merge records
integration, not mathematical acceptance.

For substantial research, read [autoresearch.md](autoresearch.md) and coordinate
through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8).
Reviewers can use the [review guide](docs/REVIEWING.md) and
[exact-head review ledger](docs/REVIEW_LEDGER.md). Document source revisions,
proof status, reproduction commands, and remaining gaps. The
[contribution record](CONTRIBUTORS.md) retains credit for useful partial and
superseded work as well as incorporated results.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
