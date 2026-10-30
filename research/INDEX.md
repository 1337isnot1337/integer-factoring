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

- [Simultaneous finite-seed endpoint](simultaneous-seed-endpoint.md): one machine attains an enumerable seed family's limiting time and workspace exponent, with no generation-speed or field-uniformity assumption.
- [Pieces-to-factoring transfer](pieces-to-factoring.md): conditional crossing-parameter transfer, with full bit-cost accounting.
- [Scale-one modular transfer](scale-one-modular-transfer.md): fixed-field transfer and machine costs.
- [Sparse divisor nonzeroness](sparse-divisor-nonzeroness.md): detection of divisors in weighted polynomial counts.
- [Modular crossing](modular-crossing.md): normalized finite-field pieces and quotient construction.
- [Modular angle factors](homogeneous-central-modular-symbols.md): Boolean-zeta symbols, a three-block polynomial-time Pfaffian, reciprocal coefficients, and exact parity support.
- [Eighth-order modular residual](homogeneous-central-eighth-jet.md): an all-size even-energy identity, explicit determinant-minor pairing, diagonal state charge, and a sharp weight-32 witness.
- [Complete homogeneous profiles](homogeneous-central-profiles.md): exact larger cut ranks, Taylor coefficients, and independently checked native polynomial Pfaffians.
- [Complement-shell crossings](complement-shell-crossings.md): explicit identity-Gram labels, exact sparse support, a universal four-state tensor, and an all-size paired reciprocal-pencil identity.
- [Complete shell rank blocks](complement-shell-rank-blocks.md): exact syndrome direct sums through 496 labels and protected central minors for fixed-length full shells.
- [Three-parameter Pfaffian crossings](pfaffian-ternary-edge-family.md): arbitrary fixed-length codes, a complete uniform edge space, exact label retracts, and binary/four-element-field seed exclusions.
- [Pfaffian boundary sections](noncommutative-crossing-pfaffian-boundaries.md): all-size section and Gram formulas, the universal correction factor, and an all-field rank-cost certificate for four pairs.
- [Sparse isotropic Pfaffian sections](noncommutative-crossing-pfaffian-isotropic.md): explicit endpoint cancellation at every size from six pairs, four-term lifts from seven pairs, exact active correction ranks, and a six-piece crossing curve.
- [High-mode residue folding](residue-folded-channels.md): exact ambient dimensions, finite incoming and outgoing factors, both complete channel counts, and eventual energy periodicity.
- [Complete folded Gram channels](folded-gram-channels.md): all-angle cancellation inside one Gram form and a constructive polynomial-time exact Johnson-rank algorithm.
- [Weighted Gram aggregation](folded-gram-aggregation.md): exact energy-preserving multiplicities, Gaussian-binomial evaluation, and 542 complete finite energy-bound calculations.
- [Homogeneous modular selection](homogeneous-central-selection.md): one operator weight at every label width, selected rectangular vacuum minors, and an effective uniform half-rate construction.
- [Noncommutative group-word crossings](noncommutative-crossing-orbits.md): exact all-size construction, compressed constant cut, optimal rank-28 central shears on 64 labels, and complete-profile limitations.
- [Coded central quotients](coded-central-crossings.md): binary-code group construction, exact constant-cut compression, matroid bounds, and an all-coefficient row filtration estimate.
- [Graded Hopf crossings](noncommutative-crossing-graded-hopf.md): exact noncommutative families, homogeneous row bounds, compressed translation cuts, and complete finite recombination and gauge audits.
- [Cumulative bosonic crossing profile](bosonic-hermite-profile.md): exact moment identities, a Virasoro rank formula, and complete pieces with independently collapsed prefix and suffix channels.
- [Homogeneous duals](homogeneous-duals.md): alternative algebraic pieces and their checked identities.
- [Near-ones label families](near-ones-labels.md): label construction and rank experiments.
- [Binary-prefix top channel](binary-prefix-top-channel.md): structural analysis of channel widths.
- [Primary central Gram matrix](primary-central-gram.md): source-linked rank calculations.
- [Half-piece ranks](half-piece-rank.md): finite crossing-piece rank evidence.

## Obstructions and work beyond 1/2

- [Directions beyond half](post-half-directions.md): experiments, scope, and unresolved questions.
- [Quadratic charges and label entropy](binary-quadratic-charge-shells.md): exact residue-syndrome characters and dimension bounds for joint eigenspaces of all low quadratic charges.
- [PBW degeneration of Hopf crossings](noncommutative-crossing-pbw-degeneration.md): a common four-leg specialization proves a finite rank-sum obstruction for the full PBW-domain coefficient spaces, including field extensions and filtered local gauges.
- [Two-stage group degeneration](noncommutative-crossing-group-degeneration.md): augmentation weight followed by ordinary PBW length excludes full-domain class-two word crossings with quadratic controls; cubic controls have an explicit failure of the first filtration.
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
