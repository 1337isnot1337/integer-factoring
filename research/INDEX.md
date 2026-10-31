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

- [Complete Schur queries](schur-fermion-queries.md): division-free small determinant formulas for modes and complete crossing coefficients; polynomial query cost with exact independent checks, but no new rank bound.
- [Schur hook residues and central ranks](schur-hook-residue-central.md): all-weight hook identities give rank at most `p(2h)-3`, attained at weights 8 and 16 by complete exact computations.
- [Schur digit reduction and a fixed Gaussian seed](schur-gaussian-digit-reduction.md): exact local maps, retained energy projections, explicit finite coefficients, and a four-label occupation-sector crossing.
- [Sparse singleton holonomy](sparse-holonomy-singleton-profile.md): a complete all-size Smith profile for zero-and-singleton indices, with a Gaussian encoding for general down-sets.
- [Sparse down-set Schur crossings](schur-downset-crossing.md): an exact polynomial kernel, energy-slice endpoint proof, reciprocal identities, and rigorous aggregate profile limits.
- [Offset interval reduction](schur-offset-interval-reduction.md): two independent interval scales reduce to a four-label Gaussian kernel with explicit finite coefficient pieces and polynomial overhead.
- [Two-block holonomy decomposition](schur-two-block-holonomy.md): two unequal Boolean intervals split into scalar Gaussian sectors, with direct polynomial seeds beyond power-of-two parameters.
- [Fixed-mode crossings and independent occupations](fixed-mode-crossings.md): arbitrary index sets give exact typed crossings, including unequal species counts, with arbitrary mode heights and a rational product-pressure criterion.
- [Certified three-mode pressure improvement](three-mode-pressure-improvement.md): asymmetric species occupations give the exact bound `rho < 741/1000`, with all-minor replay and outward integer interval arithmetic; the post-half target remains open.
- [One-species interval spectra](schur-one-species-interval.md): exact odd orders, rational certificates through 63 modes, and an analytic description of approach to the half-rate threshold from above in the displayed gauges.
- [Quantitative interval frontier](quantitative-interval-frontier.md): exact certificates through 16,383 modes and a sharp all-parity asymptotic for the gap above the half-rate threshold, with bounded-mantissa integer arithmetic.
- [Greedy Pascal complements](schur-pascal-complement.md): an exact binary-linear-algebra formula for the infinity spectrum of every finite mode set, with a traversal that skips large index gaps and an all-length interval corollary.
- [Constant-flag crossings](constant-flag-crossings.md): polynomial Gaussian crossings from two arbitrary constant flags and a unitriangular endpoint; an exact four-mode, fixed-gauge screen covers all 4,096 lower flags over the four-element field.
- [Positive-energy perturbations](positive-energy-quartic-endpoints.md): an index-sum criterion for valid nonvanishing endpoint changes and rational local coefficients whose poles cancel after every fixed-type projection.
- [Localized quartic profiles](localized-quartic-profiles.md): exact creation-annihilation determinant identities, complete zero spectra, and a checked rank-deficient residual endpoint for a non-Gaussian perturbation family.
- [Polynomial residual codes](residual-quartic-codes.md): a polynomial eight-letter candidate with a surviving quartic correction, an all-power typed endpoint proof, and an exact finite classification of constant binary endpoint restrictions.
- [Coded arbitrary-mode families](coded-quartic-mode-families.md): arbitrary finite mode sets and five independently chosen Laurent coefficient functions preserve the polynomial typed endpoint; a reusable constructor covers sixteen-letter examples.
- [Conjugate phase codes](conjugate-residual-phase-codes.md): two conjugate copies cancel a residual phase, with original-variable polynomiality after common types and three explicit eight-letter polynomial crossings that descend to the binary field.
- [Fixed field descent](fixed-field-descent.md): a linear projection after the complete typed tensor power transfers fixed-extension coefficient certificates to the base field with unchanged physical alphabet and only a constant cost factor.
- [Residual three-letter obstruction](residual-quartic-three-letter-obstruction.md): a matrix-space compression proof excludes single-copy constant three-letter crossings over every characteristic-two field; identical codes on both crossing pairs always give a singular three-letter transfer.
- [Signed residual codes](signed-residual-codes.md): a signed exterior model extends the constant three-letter obstruction to every field and admits an exact three-conjugate residual phase repair over the 27-element field; the claims stop at finite endpoint identities.
- [Asymmetric involutive core](involutive-core-crossings.md): a general cancellation identity gives polynomial candidates on `3 * 2^(n-2)` letters with a common-type crossing endpoint, including six-letter and twelve-letter examples and four integer deformation parameters.
- [Coordinate core crossings](coordinate-core-crossings.md): coordinate restrictions retain a filled low-mode background; independent matching types on each crossing pair cancel genuine local poles and give polynomial crossings for arbitrary finite mode sets, with exact rational profile normalization.
- [Independent paired-type pressure](paired-type-crossing-pressure.md): rectangular alphabets and independent type distributions transfer to equal-alphabet crossings, with an exact rational Gibbs criterion and a product cost using separate coefficient indices.
- [Rectangular residual repair](rectangular-residual-repair.md): four-by-three extreme core alphabets admit explicit rational repair, an unrestricted two-parameter inverse, finite polynomial mixed examples, and a smaller height-triangular repair with a checked non-Gaussian witness.
- [Even-core postselection](even-core-postselection.md): an explicit four-letter polynomial family has a strict endpoint, a Gaussian zero-parameter control, and two Wick identities that exclude Gaussianity for every nonzero rational parameter path.
- [Two-output Gaussian locus](two-parameter-even-core.md): independent output parameters enlarge the four-letter family; exact Wick identities and an explicit factorization classify its Gaussian members over every characteristic-two coefficient field.
- [Four-boundary Gaussian surface](four-parameter-even-core.md): the complete four-parameter Gaussian locus consists of the coordinate hyperplanes and a scalar copy of the same exceptional tensor, proved by exact multivariate coefficient identities.
- [Optional native F4 search](native-weighted-smith-f4.md): a standalone compiled profile engine uses exact field arithmetic and support-block decomposition, with independent rational controls and explicit safeguards for exploratory floating weights.
- [Optional rectangular Smith search](native-rectangular-smith-search.md): a standalone binary engine supports unequal crossing-pair alphabets, independent type optimization, exact denominator restoration, and independent rational-weight controls.
- [Exact type-optimized pressure](type-optimized-crossing-pressure.md): arbitrary rational letter frequencies reduce to one rational certificate inequality; a Gaussian physical gauge gives `rho < 61/80`, still above the post-half target.
- [Polynomial crossing interpolation](polynomial-crossing-interpolation.md): prescribed local row lattices give a complete finite linear synthesis problem, including a classified 64-element unrestricted fiber and a non-Gaussian control.
- [Tensorized Smith sections](tensorized-smith-sections.md): an explicit 220-dimensional affine family contains new two-site directions; common nilpotent residues certify identical profiles for every member at the stated gauges.
- [Weighted Smith separation oracle](weighted-smith-separation-oracle.md): exact minimum-valuation pivoting finds every determinantal minimum and supplies the constraints for convex fixed-type gauge optimization.
- [Cancellation module automata](cancellation-module-smith-automata.md): finite-state certificates that retain auxiliary-product cancellation, a monotone precision hierarchy, and a complete certificate limit for free-boundary word rates.
- [Minimal representations and boundary completeness](minimal-boundary-module-completeness.md): exact reachable and observable reduction makes the finite module hierarchy complete for fixed-boundary grouped-word rates.
- [Exact rational module orbits](exact-rational-module-orbits.md): finite exact lattice orbits compute grouped-channel word pressure, including rank-deficient states; a companion example has explicit small transition matrices.
- [A complete crossing with a strict module gain](cancellation-module-crossing-example.md): exact rational upper and lower potentials separate the module certificate from its matching capped min-plus certificate.
- [Fixed block Smith profiles](block-smith-crossings.md): submultiplicative profile bounds and an exact three-position improvement for the ordered-pairing control.
- [Uniform degree-one Smith obstruction](linear-smith-barrier.md): coupled endpoint ranks force `B >= q^3` over every field and alphabet; sharp binary bounds and all 65,536 scalar binary pencils are checked.
- [Full varying-image cubic spaces](cubic-full-image-span.md): exact 786- and 1,225-dimensional coefficient spaces over characteristic two, and smaller matrices for exact column ranks; no strict seed is supplied.

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

- [One-sided charge pruning](one-sided-charge-pruning.md): common-type restriction removes every strictly one-directional charge perturbation; both-cut Smith majorization gives a matching untyped obstruction.
- [Typed Gaussian generic blocks](typed-gaussian-generic-blocks.md): arbitrary coefficient recombinations have cost `Omega(D^2)` at fixed nonvacuum occupations, with exact finite controls and a separate growing-density bound.
- [Exact typed charge-block compression](typed-gaussian-charge-blocks.md): polynomially many orbit blocks compute every coefficient-combination rank, including complete binary searches through 1,260 physical labels and independent dense controls.
- [Sparse bulk species-gauge floor](sparse-bulk-species-gauge-floor.md): an exact entropy dual proves `rho >= 3/4` in the stated common-density limit model; no exclusion of general sparse families or unbounded energy slopes is claimed.
- [Controlled monomial Smith floor](controlled-permutation-smith-floor.md): a fiber-span argument extends the nilpotent-path obstruction to dependent, noncommuting, and weighted monomial controls.
- [Qutrit Weyl crossing](qutrit-weyl-crossing.md): a complete auxiliary crossing with an exact all-length constant-coefficient rank formula and a scoped raw-piece obstruction.
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
