# Deterministic integer factoring bounds

[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

This repository studies worst-case deterministic factoring. Given an integer `N >= 2` in binary, the algorithm outputs its smallest prime factor. Let `n` be the input bit length. Runtime is measured in finite-alphabet Turing-machine steps.

## Main theorem: logarithmic saving exponent 1/2

The strongest complete result on `main` is a deterministic algorithm with runtime

```text
T(n) <= C * 2^(C*n/(sqrt(floor(log_2 n))+1))
     = exp(O(n/sqrt(log n))).
```

Here `C` is an absolute constant. This gives saving exponent **`alpha = 1/2`** under the [Big O Speedrun scoring rule](https://bigospeedrun.com/integer-factoring).

The Lean endpoint is [`IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`](formal-appendix/HalfTheorem.lean), which proves the pinned `Factoring.SqrtLogLogConjecture`. Big O Speedrun independently rebuilt commit [`e91442d3`](https://github.com/1337isnot1337/integer-factoring/commit/e91442d3) and accepted the result with Lean and nanoda. See the [review report](https://github.com/1337isnot1337/integer-factoring/issues/28#issuecomment-6533697488) and [verifier run](https://github.com/bigospeedrun/verifier/actions/runs/41108289491).

The construction combines the clique evaluator from [Anthropic's Appendix B](https://github.com/anthropics/formal-math) with the multiplication-network reduction by [Isaac Rudich and Louis-Martin Rousseau](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring). Read the [written proof](research/half-rate-factoring.md), [constructor](research/half-top-constructor.md), [evaluator audit](research/half-evaluator-audit.md), and [Lean proof map](formal-appendix/README.md).

This is an asymptotic worst-case bound, not a practical factoring speedup.

## Results and evidence

| Claim | Status | Reference |
| --- | --- | --- |
| Saving exponent `1/2` | **Lean-verified; accepted by Big O Speedrun** | [Theorem](formal-appendix/HalfTheorem.lean) |
| Saving exponent above `1/2` | **Not established** | [Per-cube audit](research/two-axis-cube-regrouping.md); [global-pairing audit](research/two-axis-global-pairing-repair.md) |
| Finite factorization experiments | Evidence for the stated inputs only; no uniform bound follows | [Latest checkpoint](research/beyond-half-research-checkpoint.md) |

## Current research frontier

**Arithmetic route.** The project has deterministic prime-form generation, smooth-norm existence results, relation recovery, and determinant-free factor extraction. An adaptive selector factored a fresh 144-bit semiprime after 8,232,730 norm tests. The open step is a uniform bound on the number of candidate powers needed. Without that bound, the finite result does not improve the worst-case exponent. See the [research checkpoint](research/beyond-half-research-checkpoint.md) and [adaptive selector report](research/adaptive-power-window-selection.md).

The [Fourier window analysis](research/window-orbit-normalization.md) separates the actual target-cycle mean from nonzero-phase error. Exact [cycle and power controls](research/entry-moment-cycles.md) distinguish smooth density from clustered hits. Counting [changes of the relation row](research/relation-change-event-density.md) retains event density at least `delta/(2K_p-1)` uniformly, removing a separate raw-cluster-length assumption. Density for arbitrary target subgroups and uniform event-moment estimates remain open.

An [even-product density theorem](research/small-index-smooth-density.md) now proves original smooth density at least `exp(-O(b/log b))` on every subgroup of index at most `D^kappa`, for fixed `kappa<1/4`, under the principal split-prime estimate. With the full pinned character estimate, every coset is covered for `kappa<1/64`. Covered prime orbits inherit row-change density at the same scale. Actual target indices and short-window moments remain open.

A [sparse faithful-character reduction](research/sparse-faithful-character.md) replaces the full character-order family by one coefficient vector at each eligible prime level. Exact code distance, polynomial sparse decoding and classical lower-sieve prime acquisition preserve every fixed logarithmic-saving exponent, conditionally on a faster composite-index coefficient engine. That engine remains open.

A [connected-event criterion](research/connected-event-window-criterion.md) gives an explicit window linear in inverse event density, a cumulant scale, and the logarithm of a period bound, conditional on factorial cumulant estimates through that order. The arithmetic estimates and uncovered target indices remain open; saved-cycle controls do not certify the criterion or a stronger complete bound.

A [small-norm event theorem](research/small-norm-event-separation.md) preserves the covered-subgroup density scale while forcing local event separation and explicit connected correlations. An [exact word-boundary representation](research/event-word-boundary.md) reduces window estimates to a static target-prime-free smooth set by bounded phase transport. It also gives annular character sums and a contributed finite recurrence/Mellin multiplier. Long-range correlations and actual target indices remain open.

A [pooled-rank alternative](research/clustered-hits-and-pooled-rank.md) can certify enough relations before every target prime has a smooth hit. Fixed integer rows repeat at most `2K_p-1` times, but different rows need not be independent. The finite pool uses more total probes in the checked example; a uniform rank-window bound remains open.

[Character-lift guards](research/power-character-lift-guards.md) explain one finite rank stall exactly and certify long rough gaps. The existing constant-bias hypothesis defeats this specific obstruction on sufficiently large valid bases; it still supplies no quantitative rank window.

A [near-identity positivity criterion](research/near-identity-positive-minorants.md) would bound the raw selector from a positive Fourier minorant. A [uniform arithmetic norm gap](research/uniform-ambient-minorant-obstruction.md) now excludes a subexponential inverse ambient mean for every sufficiently large discriminant, under the same pinned analytic input. It covers subexponential cofactor caps and multiplier choices, and excludes the one-sided negative-spectrum repair. Subgroup-specific certificates and the original target-cycle questions remain open.

**Tensor and crossing route.** This work seeks a complete crossing below the required pressure threshold or a different valid transfer. The proposed two-axis cube regroupings fail their adaptive relabeling steps: the [per-cube audit](research/two-axis-cube-regrouping.md) and [global-pairing audit](research/two-axis-global-pairing-repair.md) retain exact counterexamples, including coefficients of the actual Gaussian family. The local Gaussian profiles remain valid; no stronger complete exponent follows. See [directions beyond the half-rate result](research/post-half-directions.md).

The [box-minor theorem](research/bare-seven-to-twelve-letters.md) already excludes the bare, ungauged two-place Smith criterion for every alphabet `q<=12`, at every degree and over every field. The old [seven-letter asymmetric profile](research/seven-letter-asymmetric-profile.md) is excluded by the later adjugate condition at every integer ramification. Physical gauges, collective recombinations, actual coefficient ranks and other transfer criteria remain separate questions.

A [confluent seven-letter quantum-spin construction](research/confluent-quantum-spin.md) has a rank-one specialization, but an exact determinant budget excludes its bare two-place target over every finite coefficient extension. That fixed degeneration is stopped.

A [noncommuting nilpotent Gram test](research/nilpotent-gram-crossing.md) also preserves the endpoint and both cold ranks one, but its row Smith profile forces pressure above the required threshold over every finite coefficient extension. This fixed construction is stopped.

A [color-parity sector theorem](research/color-distinct-sector-obstruction.md) excludes every homogeneous antipodal cube retaining separate color parity at `q>=34`, over every field: its common-scale, unweighted three-axis profile pressure is at least `3/4>1/sqrt(3)` for every common physical type. This includes the one-hot Gaussian. Closed odd-color sectors and the existing scalar polynomial certificate prove the bound without a generic profile computation; changing Gaussian weights while preserving color parity cannot escape it.

For the quartic charge-fusion construction in [PR #97](https://github.com/1337isnot1337/integer-factoring/pull/97), an [exact global creation identity](research/charge-creation-square.md) removes the entire quadratic gate term at `q=8`. A binary relation-code formula controls that term at larger alphabets; higher cut ranks and favorable complete pressure remain open.

For the fixed two-word audit family, an [all-block functional profile](research/cube-global-functional-profile.md) also excludes repairing globally favorable common-orientation products by scalar reweighting. A [recombination-before-orientation check](research/recombination-before-orientation.md) finds that the unique recombination of the 22 pieces sent to E4 by a fixed pairing has width sum 12, above the 5/8 budget; endpoint restoration is possible only at effective rho=1. Larger native alphabets, other families and coupled blocks remain open.

An [exact residue-system backend](research/precision-efficient-residue-indicators.md#5-implemented-compilers-and-the-system-count-bound) retains counts at growing 2-adic precision and handles collected unit inverses. A complete integer factoring query within that polynomial interface remains open.

A [Witt-trace word interface](research/witt-trace-native-word-circuits.md) gives a complete two-input factor/prefix predicate using explicitly charged exact divisions. Fast contraction of that extended circuit family remains open; the polynomial carry backend does not supply it.

A [one-trace inverse predicate](research/one-trace-inverse-words.md) determines the cofactor from one half-width candidate, with a division-free inverse polynomial and a certified no-wrap interval test. It preserves the complete smallest-factor reduction and charges all precision; a fast contractor remains open.

A [divided-power extension](research/divided-power-trace-lifting.md) handles polynomial trace constraints at growing precision, with its coefficient-algebra size charged explicitly. An [Artin--Hasse alternative](research/artin-hasse-trace-projection.md) uses classical binary characters and charges extra precision for Fourier projection. Both provide exact counts for their represented input class; a complete factoring representation remains open.

**Semiprime transform route.** Sanftenberg's D1 work offers a deterministic boundary contractor and a packed factor-moment scalar as candidate evaluator components. Its stated runtime target is also `alpha = 1/2`, and the hand-proved, not independently reviewed multiplier theorem in PR #73 maps D1's output-pinned bit graph to linear treewidth. Dense elimination is therefore blocked under that theorem; structured contraction and the full primitive-machine cost remain open. See the [source audit](research/riemann-gaussian-semiprime-source-audit.md) and [multiplier treewidth result](research/multiplier-treewidth-and-relation-cycles.md).

The [research index](research/INDEX.md) is the full catalogue of proofs, constructions, finite searches, obstructions, and open questions. It contains the detailed research log; this README is only a summary.

## Verification

From the repository root:

```sh
make verify               # Default Python tests, exact certificates, and Markdown links
make verify-docs          # Markdown-link check only
make verify-half-appendix # Lean 4.33.1; complete 1/2 theorem
make verify-k6            # Lean 4.35.0; separate k6 theorem
make verify-formal        # Original supporting Lean lemmas
```

The default `make verify` gate does not rebuild the Lean proofs. Research notes list their own focused replay commands and dependencies.

## Repository map

| Path | Contents |
| --- | --- |
| [`research/`](research/INDEX.md) | Research notes, full index, and open questions |
| [`formal-appendix/`](formal-appendix/README.md) | Lean proof of the `1/2` theorem |
| [`variants/k6/`](variants/k6/FactoringK6/Certified.lean) | Separate earlier k6 proof |
| [`formal/`](formal/) | Supporting Lean lemmas and checks |
| [`evidence/`](evidence/) | Finite certificates and verification records |
| [`tools/`](tools/) and [`tests/`](tests/) | Checkers, reproduction scripts, and tests |
| [`docs/`](docs/START_HERE.md) | Onboarding and review guides |

## Contribution and attribution

Start with [docs/START_HERE.md](docs/START_HERE.md) and [CONTRIBUTING.md](CONTRIBUTING.md). Useful contributions include proof corrections, independent replays, partial or negative results, finite certificates, checkers, tests, documentation, and maintenance. A stronger factoring exponent is not required.

Work directly on `main` by default. Make focused commits and push useful checkpoints frequently. Merge useful pull requests by default. A merge records integration; it does not imply independent mathematical review or acceptance of a stronger bound.

Use [GitHub Discussions](https://github.com/1337isnot1337/integer-factoring/discussions), starting with the [coordination thread](https://github.com/1337isnot1337/integer-factoring/discussions/51), for overlapping research and agent handoffs. Issue #8 is a historical archive. Preserve source attribution and contributor credit, and disclose substantial AI assistance.
