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
| Simultaneous time and workspace at an enumerable seed endpoint | Ordinary conditional theorem; no larger numerical exponent | [Bounded selection and direct type channels](research/simultaneous-seed-endpoint.md) |
| Polynomial-time modular angle-factor queries | Ordinary exact construction; channel ranks unchanged | [Three-block Pfaffian and transpose identity](research/homogeneous-central-modular-symbols.md) |
| Uniform eighth-order residual identity at even state energy | Ordinary exact theorem; first contact term identified | [Determinant-minor pairing and residue charge](research/homogeneous-central-eighth-jet.md) |
| Homogeneous modular duals at every label width | Ordinary uniform construction; retains exponent `1/2` | [Rectangular vacuum minor selection](research/homogeneous-central-selection.md) |
| Noncommutative 64-label constant-cut rank reduced from 31 to 28 | Exact finite construction and exhaustive search; complete profile does not improve the score | [Group-word crossings and all central shears](research/noncommutative-crossing-orbits.md) |
| Cumulative bosonic crossing profile | Ordinary exact identities and finite construction; no stronger exponent | [Moment pairing, exact rank, and both channel cuts](research/bosonic-hermite-profile.md) |
| Crossing families from binary codes | Ordinary all-size construction and one-cut bounds; complete column profile open | [Coded central quotients](research/coded-central-crossings.md) |
| Explicit complement-shell crossings | Ordinary all-size construction and exact four-state tensor reduction; tested seeds do not improve the score | [Identity Gram, sparse support, and local tensor](research/complement-shell-crossings.md) |
| Graded noncommutative Hopf crossings | Ordinary all-size construction with explicit bounds for both cuts; no stronger exponent | [Generator substitution, scalar factorizations, and complete profiles](research/noncommutative-crossing-graded-hopf.md) |

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-constructions # New crossing and resource-interface certificates
make verify-native-crossings # Optional C++ screen with independent Python audit
make verify-post-half     # Larger optional rank tests
make verify-docs          # Local Markdown-link check
```

The Python checks use Python 3 and its standard library. Some optional checks use NumPy for speed. Native crossing checks also require a C++17 compiler.

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

The [simultaneous endpoint theorem](research/simultaneous-seed-endpoint.md) preserves the limiting seed score in both logarithmic running time and working space, even when useful certificates have arbitrary gaps, take arbitrarily long to generate, or use different coefficient fields. It uses bounded certificate selection and explicit amplified factor channels. It does not supply a better seed or an exact constant-factor bound at an endpoint approached without a convergence rate.

A [three-block Pfaffian formula](research/homogeneous-central-modular-symbols.md) evaluates the existing characteristic-two angle factors in `O(h^4)` bit operations at partition weight `h`. An independent determinant-sum formula gives their exact reciprocal relation and shows that all odd angle and crossing-shift indices vanish. This complements the earlier four-block formula for complete pieces; it does not decrease their ranks or improve the factoring score.

The [eighth-order residual theorem](research/homogeneous-central-eighth-jet.md) identifies the first possible nonconstant Taylor coefficient of the homogeneous modular profile at even state energy. It factors into an explicit parameter pairing and a diagonal residue charge. A complete weight-32 example makes this coefficient nonzero, so the order-16 behavior of smaller examples does not extend uniformly. [Complete native profiles](research/homogeneous-central-profiles.md) record both cut ranks and independent operator comparisons.

[Complement shells](research/complement-shell-crossings.md) give explicit homogeneous duals with an identity vacuum Gram. An equality argument reduces their full crossing tensors to products of one universal four-state local tensor and leaves at most `2L+1` pieces for labels of length `L`. This proves the reciprocal-pencil reduction of the paired subfamily at every size and exposes a larger family for rank analysis. [Quadratic charge characters](research/binary-quadratic-charge-shells.md) describe precisely how imposing all low-weight charges limits the available label dimension. Neither result establishes a stronger factoring exponent.

The [homogeneous selection theorem](research/homogeneous-central-selection.md) constructs duals of one weight at every binary label width. A large nonsingular minor exists in a homogeneous block of the proved source pairing, even when the same-weight Gram matrix is singular. The selected operator and state weights differ by `O(sqrt(h) log h)`, so the construction retains the uniform half-rate bound. Its raw-parameter Pfaffian makes complete finite profiles accessible without listing intermediate partition spaces.

[Noncommuting group words](research/noncommutative-crossing-orbits.md) give a different exact crossing family. An exhaustive screen of all `33554432` normalized functions on `F_4^2` finds the optimal constant-cut rank 28 in that central-shear family. The other pieces remain too expensive, including after all binary recombinations of the displayed winning profile. This is a finite rank improvement and a reproducible search result; the factoring score remains `1/2`.

The [cumulative bosonic profile](research/bosonic-hermite-profile.md) uses a vanishing first moment and a scalar second moment to build new exact pieces. Its scalar parameter pairing has rank `p(h)-p(h-2)`, proved using an explicit Virasoro operator. Prefix and suffix formulas give both channel factorizations without dropping contact terms. These pieces retain the leading half-rate channel scale; their complete ranks are measured separately in the finite examples.

[Coded central quotients](research/coded-central-crossings.md) let a binary code control the noncommuting group relations. A column-matroid formula bounds the constant cut, and a weighted filtration bounds every row coefficient. Some parameter ranges survive both the constant-cut test and the known quotient obstruction. A matching bound for the complete column profile is still needed.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
