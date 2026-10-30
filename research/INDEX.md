# Research index

This index covers research notes on `main`. Open PRs can contain later work. Each note states its own assumptions, proof scope, and evidence. An obstruction to a specific construction is not a factoring lower bound.

## Proved machine bounds and their foundations

- [Half-rate factoring composition](half-rate-factoring.md): ordinary proof of the `1/2` machine bound; its final Lean theorem is in [`formal-appendix/HalfTheorem.lean`](../formal-appendix/HalfTheorem.lean).
- [Top constructor and codec](half-top-constructor.md): deterministic table construction and exact input translation.
- [Evaluator and encoding audit](half-evaluator-audit.md): recursive evaluator, computable pad, and charged bit costs.
- [Verified k6 machine](verified-k6.md): Lean proof of the earlier `1/30` bound.
- [Sharper k6 rate](k6-sharp-rate.md): ordinary rate derivation above `1/29` for the k6 program; this sharper exponent is not Lean-certified.
- [k6 source audit](k6-source-audit.md): hashes, modified source scope, and compiler preservation.
- [Baseline audit](baseline-audit.md): source recurrence and the unchanged reference's machine-cost analysis.

## Transfers and proposed stronger evaluators

- [Integral Schur crossing](integral-schur-crossing.md): the complete partition construction has multiplier one over every characteristic; full homogeneous duals exist at prime-power weights.
- [Effective integral construction](effective-integral-schur-transfer.md): division-free Schur determinants, memoized modes, and polynomial coefficient sizes give a uniform `2^O(lambda)` construction.
- [Finite min-plus channel automata](tropical-smith-channel-words.md): normalized capped cost vectors count each channel word once, with exact rational potential certificates.
- [Ternary nilpotent resolvents](ternary-nilpotent-resolvents.md): exact characteristic-three endpoints, a small determinant nilpotence test, and independent Toeplitz verification of native profiles.
- [Degree-sensitive auxiliary paths](degree-sensitive-smith-paths.md): entrywise series orders constrain channel walks; rational positive potentials certify a stronger finite criterion, with exact zero tests from bordered minors.
- [General-characteristic trace-shear parameter](trace-shear-general-characteristic.md): extends the existing Fierz common-boundary identity, checks both cut signs, and distinguishes identity normalization from preservation of the endpoint. The profiles remain noncompetitive.
- [Fixed-precision count lifting](fixed-precision-count-lifting.md): sorted witness tuples recover prime-power count digits; local nilpotent corrections preserve every fixed field seed and Smith-profile leading exponent asymptotically.
- [High-order characteristic selection](high-order-characteristic-selection.md): an elementary polynomial prime bound prevents cancellation of the quadratic root detector, with an explicit uniform-evaluator requirement.
- [Exhaustive auxiliary lattice gauges](auxiliary-lattice-gauge-search.md): canonical two-dimensional lattice representatives cover every rational gauge; neighbor inequalities certify global minima at a fixed evaluation point.
- [Unitary-orbit factor recovery](unitary-orbit-factor-recovery.md): coprime tuples and perfect-power preprocessing give polynomial-overhead recovery in every fixed positive characteristic, with unique inverse witnesses.
- [Quadratic-congruence parity recovery](quadratic-congruence-factor-recovery.md): a scalar square root of one yields a proper gcd factor; a fused quadratic-plus-linear verifier has the source's original graph shape.
- [Balanced physical degree gauges](balanced-physical-gauge-transfer.md): independent four-leg weights become known coefficient shifts on a balanced label set of size `q^t/poly(t)`.
- [Local Smith-profile crossing certificate](smith-profile-crossing-transfer.md): two exact local valuation spectra certify all-length rank bounds for ordered polynomial matrix products.
- [Independent auxiliary rational gauges](auxiliary-gauge-smith-transfer.md): different telescoping basis changes in the two cuts enlarge the Smith search with constant boundary costs.
- [Recovery in every fixed prime characteristic](fixed-characteristic-factor-recovery.md): higher exponent digits and unique cyclic-orbit representatives give complete factorization with polynomial query overhead and linear total arithmetic length.
- [Parity and power-divisor recovery](parity-power-divisor-recovery.md): characteristic-two counts detect exponent bits and recover a complete factorization with polynomial query overhead, without a sparse-recovery square-root floor.
- [General weighted network compiler](weighted-network-top-compiler.md): deterministic-message networks with linear total weight and sublinear maximum gate weight compile to linear-volume top cliques; variable multiplication output ports preserve unique witnesses.
- [Pieces-to-factoring transfer](pieces-to-factoring.md): conditional crossing-parameter transfer, with full bit-cost accounting.
- [Scale-one modular transfer](scale-one-modular-transfer.md): fixed-field transfer and machine costs.
- [Sparse divisor nonzeroness](sparse-divisor-nonzeroness.md): detection of divisors in weighted polynomial counts.
- [Modular crossing](modular-crossing.md): normalized finite-field pieces and quotient construction.
- [Homogeneous duals](homogeneous-duals.md): alternative algebraic pieces and their checked identities.
- [Near-ones label families](near-ones-labels.md): label construction and rank experiments.
- [Binary-prefix top channel](binary-prefix-top-channel.md): structural analysis of channel widths.
- [Primary central Gram matrix](primary-central-gram.md): source-linked rank calculations.
- [Half-piece ranks](half-piece-rank.md): finite crossing-piece rank evidence.

## Obstructions and work beyond 1/2

- [Complete weight-four Schur spans](full-schur-span-recombination.md): exact integer shortest-path bounds cover every binary linear recombination with any number of summands; both spans remain above the needed finite-seed threshold.
- [Sum-preserving conjugation Smith floor](sum-preserving-conjugation-smith-barrier.md): arbitrary constant exchange-commuting block conjugations have certificate parameter at least `log_2(64/9)/4`, just above the half-rate threshold; actual rank sums remain outside the claim.
- [Explicit entanglers and scalarization controls](explicit-entanglers-and-density-scalarization.md): product-ket exclusions over algebraic closures, a surviving full-rank product functional, and a complete binary product-density screen.
- [Coherent exterior-boundary reduction](coherent-exterior-boundary-reduction.md): a common histogram projection, all-width fixed-occupation lower bounds, and a complete four-position binary coefficient-span screen.
- [Whole-root-interval parity counterexample](root-threshold-parity-counterexample.md): an exact non-perfect-power composite defeats every whole root-threshold parity count; general interval recovery remains unrestricted.
- [Directions beyond half](post-half-directions.md): experiments, scope, and unresolved questions.
- [Network tradeoff](half-top-structure.md): restriction for the current local-enumeration and wire-channel constructor.
- [Restricted-label counterexample](half-restricted-labels.md): failure of a proposed low-mode dimension reduction.
- [Restricted intermediate spans](half-restricted-span.md): quantitative obstruction, with precise selected-label conditions.
- [Half-rate recurrence barrier](half-recurrence-barrier.md): limits of the existing operators and recurrence.
- [Graded-space barrier](graded-space-barrier.md): lower bound under the stated graded-space assumptions.
- [Rectangular-grid frontier](rectangular-grid-frontier.md): limits of unequal splits under the declared-profile recurrence.
- [Top-channel rank](top-channel-rank.md): exact factor rank on the full partition domain.
- [Coefficient basis](coefficient-basis.md): scoped minimax calculation for the chosen coefficient span.

## Earlier grouping and crossing research

- [Integration notes](integration.md): earlier combined grouping/crossing proposal and open proof obligations.
- [Recurrence](recurrence.md): the original recurrence and cost accounting.
- [Grouping improvement](grouping-improvement.md): concentration argument and deterministic construction.
- [Grouping validation](grouping_validation.md): finite grouping examples.
- [Crossing compression](crossing-compression.md): explicit integer factorization and channel counts.
- [Grid-level validation](grid-level-validation.md): exact small-grid tests.

## Verification and proof artifacts

- [`../formal-appendix/README.md`](../formal-appendix/README.md) maps the Lean `1/2` proof.
- [`../variants/k6/FactoringK6/Certified.lean`](../variants/k6/FactoringK6/Certified.lean) gives the Lean `1/30` theorem.
- [`../evidence/`](../evidence/) stores finite certificates, source audits, and recorded CI output.
- [`../tools/`](../tools/) contains the checkers. Run `make verify` from the repository root.

No finite certificate alone proves a new asymptotic factoring exponent. Distinguish an exact finite identity from a completed uniform machine theorem.
