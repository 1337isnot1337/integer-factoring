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

An [entry-moment sign correction](research/entry-moment-tail-sign-correction.md) recomputes the old tables from saved flags: their earlier lower-tail labels actually measured upper deviations. The claimed general raw-versus-entry lower-tail contrast is withdrawn; periods, gaps and the complete bound are unchanged.

An [even-product density theorem](research/small-index-smooth-density.md) now proves original smooth density at least `exp(-O(b/log b))` on every subgroup of index at most `D^kappa`, for fixed `kappa<1/4`, under the principal split-prime estimate. With the full pinned character estimate, every coset is covered for `kappa<1/64`. Covered prime orbits inherit row-change density at the same scale. Actual target indices and short-window moments remain open.

A [sparse faithful-character reduction](research/sparse-faithful-character.md) replaces the full character-order family by one coefficient vector at each eligible prime level. Exact code distance, polynomial sparse decoding and classical lower-sieve prime acquisition preserve every fixed logarithmic-saving exponent, conditionally on a faster composite-index coefficient engine. Actual divisor reciprocity also gives an inversion-fixed weight-one target and recovers paired field components under the stated condition. The engine remains open.

[One chosen 2-adic completion](research/one-cyclotomic-completion.md) at explicitly polynomial precision also suffices for all inputs: algebraic norm separation and exact LLL decoding recover the whole integral coefficient from one fixed-character local value. All level, Hensel, precision and reconstruction costs preserve a hypothetical evaluator's fixed saving exponent. Its [source-pinned full default replay](evidence/one-completion-full-verification.json) passed in 1385.48 seconds. This is a conditional interface result; the original-index evaluator remains missing and the complete bound stays `1/2`.

For [promised semiprimes](research/semiprime-single-field.md), one primitive field component at every coprime q>=7, q=3 mod 4, determines unordered signed factor residues through a quadratic equation. Symmetric CRT recovers the factors with polynomial additional cost. The coefficient engine and an arbitrary-input reduction remain open.

A [connected-event criterion](research/connected-event-window-criterion.md) gives an explicit window linear in inverse event density, a cumulant scale, and the logarithm of a period bound, conditional on factorial cumulant estimates through that order. The arithmetic estimates and uncovered target indices remain open; saved-cycle controls do not certify the criterion or a stronger complete bound.

A [CM mixing audit](research/cm-mixing-event-scope.md) distinguishes published whole-packet geometric estimates from the required arithmetic prime-cycle correlations. An explicit pair of nearby CM points has opposite truncated smooth flags and forces an exponential Lipschitz interpolation cost on a family of order discriminants. Arithmetic weighted estimates and orbit-specific representations remain possible; no event-moment theorem is imported from geometric equidistribution.

A [small-norm event theorem](research/small-norm-event-separation.md) preserves the covered-subgroup density scale while forcing local event separation and explicit connected correlations. An [exact word-boundary representation](research/event-word-boundary.md) reduces window estimates to a static target-prime-free smooth set by bounded phase transport. It also gives annular character sums and a contributed finite recurrence/Mellin multiplier. Long-range correlations and actual target indices remain open.

A [near-height enumeration theorem](research/near-height-smooth-descent.md) makes the complete smaller-norm ellipse enumerable in `O(1+p^k/sqrt(D))` candidates, up to polynomial costs, after a prime-power target's leading smoothness test fails. A [rough composite-input counterfamily](research/prime-power-descent-obstruction.md) disproves strict smooth descent for ordinary prescribed scalar thresholds at every fixed logarithmic saving exponent. A [classical-sieve extension](research/full-trial-sieve-obstruction.md) makes the failure survive every subexponential trial-division cutoff. The extra ellipse points can all be imprimitive; adaptive or collective targets remain separate. No stronger complete bound follows.

A [pooled-rank alternative](research/clustered-hits-and-pooled-rank.md) can certify enough relations before every target prime has a smooth hit. Fixed integer rows repeat at most `2K_p-1` times, but different rows need not be independent. The finite pool uses more total probes in the checked example; a uniform rank-window bound remains open.

An [odd-quotient relation relaxation](research/odd-quotient-relations.md) accepts smooth-pair equality modulo a computable odd kernel, preserving all two-part factoring directions and charging the multiplier. A certified two-group image permits direct extraction; general quotient-pair extraction remains open, with generation retaining the existing analytic premises.

The odd-quotient [full default replay](evidence/odd-quotient-full-verification.json) passes at frozen `002bcf3`; later documentation and PR104/105 integrations have focused checks. This is implementation evidence, with the complete bound unchanged.

[High two-power roots of row residues](research/high-root-projection.md) give an alternative verified projection to the full supplied two-primary subgroup. The ambient bound, determinant-dependent precision and numerical root depth are charged; efficient root acquisition remains open and the complete bound is unchanged.
Its [full default replay](evidence/high-root-full-verification.json) passes at frozen `be4aae1`; the later guaranteed-positive specialization and archive prose have focused checks.

An [always-positive deep-root reduction](research/deep-class-root-projection.md) removes relation and smooth-density acquisition: any root z with z^(2^(2L))=g^(2^L) yields g's exact two-primary projection. Existing generation and early-stop extraction give poly(b)(F(b)+1) complete conditional cost; a fast root engine remains missing. Its [full default replay](evidence/deep-root-full-verification.json) passes at frozen `9877f64`, with later integration/prose checked separately and a scoped ordinary inspection of the all-input assembly.

[Character-lift guards](research/power-character-lift-guards.md) explain one finite rank stall exactly and certify long rough gaps. The existing constant-bias hypothesis defeats this specific obstruction on sufficiently large valid bases; it still supplies no quantitative rank window.

A [near-identity positivity criterion](research/near-identity-positive-minorants.md) would bound the raw selector from a positive Fourier minorant. A [uniform arithmetic norm gap](research/uniform-ambient-minorant-obstruction.md) now excludes a subexponential inverse ambient mean for every sufficiently large discriminant, under the same pinned analytic input. It covers subexponential cofactor caps and multiplier choices, and excludes the one-sided negative-spectrum repair. Subgroup-specific certificates and the original target-cycle questions remain open.

**Tensor and crossing route.** This work seeks a complete crossing below the required pressure threshold or a different valid transfer. The proposed two-axis cube regroupings fail their adaptive relabeling steps: the [per-cube audit](research/two-axis-cube-regrouping.md) and [global-pairing audit](research/two-axis-global-pairing-repair.md) retain exact counterexamples, including coefficients of the actual Gaussian family. The local Gaussian profiles remain valid; no stronger complete exponent follows. See [directions beyond the half-rate result](research/post-half-directions.md).

The [box-minor theorem](research/bare-seven-to-twelve-letters.md) already excludes the bare, ungauged two-place Smith criterion for every alphabet `q<=12`, at every degree and over every field. The old [seven-letter asymmetric profile](research/seven-letter-asymmetric-profile.md) is excluded by the later adjugate condition at every integer ramification. Physical gauges, collective recombinations, actual coefficient ranks and other transfer criteria remain separate questions.

A [confluent seven-letter quantum-spin construction](research/confluent-quantum-spin.md) has a rank-one specialization, but an exact determinant budget excludes its bare two-place target over every finite coefficient extension. That fixed degeneration is stopped.

A [noncommuting nilpotent Gram test](research/nilpotent-gram-crossing.md) also preserves the endpoint and both cold ranks one, but its row Smith profile forces pressure above the required threshold over every finite coefficient extension. This fixed construction is stopped.

A [color-parity sector theorem](research/color-distinct-sector-obstruction.md) excludes every homogeneous antipodal cube retaining separate color parity at `q>=34`, over every field: its common-scale, unweighted three-axis profile pressure is at least `3/4>1/sqrt(3)` for every common physical type. This includes the one-hot Gaussian. Closed odd-color sectors and the existing scalar polynomial certificate prove the bound without a generic profile computation; changing Gaussian weights while preserving color parity cannot escape it.

A [quartic additive-charge construction](research/quartic-charge-fusion.md) breaks separate color parity while preserving the endpoint. Its generic order-four multiplicity falls from binomial(q,2) to `2(q-1)` on additive color alphabets; the eight-color control changes 28 channels to 14. It pays degree four per occupied vertex instead of three. A [returning order-five layer](research/quartic-order-five-obstruction.md) forces restored pressure at least `7/12>1/sqrt(3)` for every q>=4096 in this construction and its stated spatial-quartic extensions. The complete score remains `1/2`.

A [fixed q8 two-axis projection](research/fixed-curve-quartic-floor.md) uses (X,Y,Z)=(x^2,y^2,xy) with one fixed leg pairing. Exact binary jets and a degree/entropy certificate give a profile floor 3/4 for every common paired type after restoring the doubled degree. This fixed curve is stopped under that criterion. Its [source-pinned full default replay](evidence/fixed-curve-full-verification.json) passed in 1981.39 seconds; later integrated contributions have separate focused checks. This is implementation verification, without independent mathematical review or a stronger complete bound.

For the quartic charge-fusion construction in [PR #97](https://github.com/1337isnot1337/integer-factoring/pull/97), an [exact global creation identity](research/charge-creation-square.md) removes the entire quadratic gate term at `q=8`. A binary relation-code formula controls that term at larger alphabets; higher cut ranks and favorable complete pressure remain open.

On Sol worker 2's fixed two-axis curve `(X,Y,Z)=(x^2,y^2,xy)`, [colour-column alternation](research/colour-alternation-joint-cold-curve.md) cancels symmetric quartic terms and preserves Gaussian exterior-minor channels at every alphabet size. Their exact minors force common-scale, unweighted paired-type pressure at least `3/4>1/sqrt(2)` for the conservative threshold `q-1>=2^90`. This result does not settle `q=8`, weighted gauges or unequal scales; increasing the alphabet in this fixed symmetric family cannot cross the target. Its [full default replay](evidence/colour-full-verification.json) passes at frozen `54acb0e`, excluding later integrations. The complete bound remains `1/2`.

A [spatial-mode colour fusion](research/spatial-mode-colour-fusion.md) changes that symmetry and fuses the whole leading colour-pair block at both cuts, with unchanged endpoint and degree costs. For `q=2^s>=4`, its exact prefix has only `5(q-1)` order-twelve channels; at q8 this lowers 91 to 35. Returning orders still force a `3/4` profile floor for this fixed q8 tensor, so that core is stopped. This does not exclude every growing alphabet in the new family or give a stronger factoring bound.

For the fixed two-word audit family, an [all-block functional profile](research/cube-global-functional-profile.md) also excludes repairing globally favorable common-orientation products by scalar reweighting. A [recombination-before-orientation check](research/recombination-before-orientation.md) finds that the unique recombination of the 22 pieces sent to E4 by a fixed pairing has width sum 12, above the 5/8 budget; endpoint restoration is possible only at effective rho=1. Larger native alphabets, other families and coupled blocks remain open.

An [exact residue-system backend](research/precision-efficient-residue-indicators.md#5-implemented-compilers-and-the-system-count-bound) retains counts at growing 2-adic precision and handles collected unit inverses. A complete integer factoring query within that polynomial interface remains open.

A [Witt-trace word interface](research/witt-trace-native-word-circuits.md) gives a complete two-input factor/prefix predicate using explicitly charged exact divisions. Fast contraction of that extended circuit family remains open; the polynomial carry backend does not supply it.

A [one-trace inverse predicate](research/one-trace-inverse-words.md) determines the cofactor from one half-width candidate, with a division-free inverse polynomial and a certified no-wrap interval test. It preserves the complete smallest-factor reduction and charges all precision; a fast contractor remains open.

A [trace-dual coordinate encoding](research/trace-dual-native-counts.md) represents each odd half-width candidate exactly once. After complete trial division through `2^t`, unweighted prefix counts fit `min(ceil(b/2)-1,ceil(b/t))` output bits, allowing polylogarithmic output precision below any fixed positive-score budget. It uses several coordinate traces of one Teichmuller input and polynomial Boolean gates for odd division, with the same small starting precision and no exact halvings. A fast sum of the unrestricted-degree circuit remains open.

With [direct Boolean hash masks](research/native-parity-hashes.md), one-bit scalar parities of that native circuit suffice for arbitrary inputs after the charged pretest. The existing polynomial hash isolation and verified bit tags recover the smallest factor with polynomial cost in the divisor bound. This precision-one provider is still open; the finite replay uses exhaustive field states.

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

The [combined PR96–PR99 replay](evidence/charge-reciprocity-full-verification.json) passed at clean frozen source `a21c955` in 888.40 seconds; its [complete log](evidence/charge-reciprocity-full-verification.stdout.txt.xz) and source/tree hashes are archived. Later changes have their separate focused scope in the [review ledger](docs/REVIEW_LEDGER.md). This verifies implementations and certificates, without establishing a stronger factoring bound.

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
