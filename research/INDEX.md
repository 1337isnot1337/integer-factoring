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

- [Arithmetic Cauchy factorial law and Smith profiles](arithmetic-cauchy-factorial-law.md): a finite-sampling moment bound proves the permanent congruence at every size, and an integral orthogonal basis proves every squared-Cauchy invariant factor over each prime-adic coefficient ring.
- [All arithmetic Cauchy powers](arithmetic-cauchy-powers.md): a factorial-normalized alternant bound proves the determinant congruence and every Smith order for each positive Hadamard power, including an explicit tunable first spectral gap.

- [Mixed-field rectangular repair](rectangular-f4-repair.md): an explicit two-parameter height-triangular correction over every characteristic-two field, a uniform simple-pole bound, and four exact mixed-F4 constructors.
- [Full-alphabet quartic repair](full-alphabet-quartic-repair.md): a rank-four projector identifies every invertible coefficient, a sparse correction retains all four states on both pairs, and fixed-size profiles support 23 or 29 parameters.
- [Bounded rectangular defect blocks](rectangular-defect-factorization.md): full-interval repairs reduce to seven and nine active fermion modes at the stated Gaussian gauges, with structurally proved generic ranks and exact rational support checks.
- [All four extreme cores](rectangular-full-core-factorization.md): ten active modes at both places give fixed 768-order blocks and 25 continuous search parameters, with exact regular-basis and full-matrix controls.
- [Coupled core transports](rectangular-core-transports.md): arbitrary invertible constant core changes preserve a certified endpoint through an exact matrix identity; cached coordinate updates and two-base scalar interpolation avoid rebuilding the Fock program.
- [Critical rectangular pressure](rectangular-critical-pressure.md): a target-sign discriminator, explicit Gaussian control, frozen-boundary limit, and reproducible starts that activate the extreme core states.
- [Bivariate core limits](bivariate-core-pressure-limit.md): eventual affine Smith orders for fixed rational two-variable matrices and a finite limiting pressure problem; four-size profile fits remain separate from an exported symbolic certificate.
- [Exact leading-face LP certificates](rectangular-leading-face-lp.md): unbounded rational primal/dual certificates isolate the all-low face in fitted 19- and 25-parameter models; uniform entry bounds remain a separate obligation.
- [Uniform entry and ray obstruction](rectangular-uniform-ray-obstruction.md): 28 symbolic coefficient witnesses replace the fitted entry bounds and support exact dual certificates for the baseline family.
- [Moving balanced-gauge obstruction](rectangular-balanced-gauge-obstruction.md): every balanced gauge vector has positive critical pressure eventually when M log s exceeds c log M with c greater than 1/sqrt(2), uniformly over unbounded weights and scales.
- [Full-alphabet leading obstruction](full-alphabet-leading-obstruction.md): 34 symbolic witnesses and 22 rational duals extend the moving-gauge theorem to both conjugate full-alphabet partial repairs and all 29 parameters.
- [All-scale balanced obstruction](all-scale-balanced-obstruction.md): the positive finite Gaussian tail removes the logarithmic scale restriction for both partial-repair families, uniformly over all balanced weights and all positive scales at sufficiently large intervals.
- [All-low Gaussian background](rectangular-low-background.md): a uniform two-cut calculation identifies the preserved conditional tensor with the interval of length M-2, including bulk charges and scalar cancellation.
- [Polynomial frozen-core rows](polynomial-frozen-core-crossings.md): one affine endpoint equation gives eleven free parameters, with locally polynomial coefficients and an exact non-Gaussian member outside the preserved-background family.
- [Polynomial core shears](rectangular-polynomial-shears.md): Laurent core changes equal to the identity at the endpoint permit every height direction, with exact cached updates and original-tensor controls.
- [Four additional quartics](rectangular-extra-quartics.md): positive-energy extreme monomials preserve the existing repair, with exact Laurent constructors and seven/nine- or ten-mode support.
- [Core monomials and opposite cubic charges](rectangular-core-monomials.md): arbitrary positive-energy monomials preserve the endpoint, while a uniform full-alphabet obstruction reduces the useful cubic search to opposite-charge pairs.
- [Critical cubic residual repairs](full-alphabet-monomial-residuals.md): four cubic pairs restore all sixteen states at the singular canonical quartic, with symbolic determinant (1+gamma+ab)^4 and a uniform height-block correction.
- [Critical-cubic pressure screens](critical-cubic-pressure-screen.md): a complete two-mode floating survey, saved paired-type starts, and reproducible drivers for cubic repairs and ramified polynomial perturbations; no threshold improvement.
- [Secondary-height repairs](rectangular-refined-repairs.md): 231 individually admissible F4 shears, dense additive flags, and Laurent quartic freedom preserve the all-length crossing and bounded support.
- [Rectangular convex pressure search](rectangular-convex-pressure.md): paired Gibbs derivatives, attaining-minor separation, generic-rank interfaces, and independent exact and numerical implementation controls.
- [Retained profile precision](retained-profile-precision.md): exported coefficient arrays carry their actual retained cutoff, with an adversarial regression for reuse in a larger elimination window.
- [Cauchy packet polynomial realization](cauchy-packet-realization.md): every degree-at-most-wire-count polynomial integral on the interpolation grid has an explicit affine packet realization, with separate canonical-lift controls.
- [Canonical alternating packets](canonical-alternating-packet.md): an all-size binary-input family has exact minimum odd-center content `L(1+log_2 L)-1`, with creation precision one.
- [Canonical packet upper bounds](canonical-packet-upper-frontier.md): exact exhaustive four-wire visibility, adaptive positive endpoint bounds, and an all-size limitation of low-order endpoint jets.
- [Finite canonical decoder thresholds](canonical-packet-decoder-frontier.md): exact parity and full-count separation for two and four wires, with an all-size binary full-count collision.
- [Full Boolean Cauchy operator](boolean-cauchy-smith-profile.md): an all-size binomial Smith profile, a one-bit inverse at the nearest centers, uniform implicit inverse entries, and an exact grouping cost for direct insertion.
- [Grouped Cauchy inverse profile](grouped-cauchy-inverse-profile.md): a weighted Boolean substitution lemma and graded elimination give the exact grouped Smith orders and normalization cost of the direct inverse.
- [Walsh twirling and Cauchy permanents](walsh-principal-permanents.md): exact discrepancy valuations, a power-of-two factorial law for even-spacing Cauchy permanents, and source-oriented grouped-profile bounds separate scalar recovery from full-operator contraction.
- [Eight-wire polynomial lattice](eight-wire-cauchy-lattice.md): an exact finite module certificate realizes every integral degree-at-most-120 polynomial and three full wire-scale tapered examples, with canonical precision still separate.
- [Dyadic affine mode reduction](one-species-dyadic-affine-reduction.md): appending common trailing binary ones gives exactly the smaller typed one-species tensor after ramification, including independent-gauge transport.
- [Pascal complement convolution](pascal-complement-convolution.md): characteristic-independent complement-order inequalities and a field-free pressure envelope for arbitrary mode sets.
- [All-field Pascal crossings](all-field-pascal-crossings.md): signed global Gaussian endpoints and Pascal-complement profiles over every field, including independent matching types.
- [Weighted Lucas-Pascal profiles](weighted-lucas-pascal-profiles.md): a unit factorization gives exact independent-height profiles by constant-field elimination, and prime-radix block traversal skips arbitrary index gaps.
- [Graded even tensor powers](graded-even-power-transfer.md): both cut profiles tensorize with explicit signs, admitting non-Gaussian endpoint-vanishing perturbations in odd characteristic.
- [Truncated symmetric Pascal crossings](truncated-symmetric-pascal-crossings.md): p-state modes in characteristic p give ordinary tensor powers, exact digit-sum spectra, and arbitrary-degree endpoint-vanishing perturbations.
- [Every-digit interval frontier](truncated-symmetric-interval-frontier.md): the displayed interval pressure stays above the half-rate threshold for every digit count, every interval length, and every scale, with a sharp fixed-digit Lambert-W gap.
- [Fixed-prime recovery source audit](fixed-prime-transfer-source-audit.md): pinned scalar factoring recovery, direct primary-source network accounting, and a byte-identical finite replay of the external dependency.

- [31 October constructive frontier checkpoint](constructive-frontier-checkpoint.md): exact finite bounds, proved construction families, reproducible floating starts, and the remaining threshold obligation after a three-hour research session.
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
- [Accelerated exact interval certificates](accelerated-interval-certificates.md): a positive Euler-series remainder removes the scan through modes and certifies an interval of 2^209-1 modes, with a rational upper bound less than 5e-123 above the unchanged threshold.
- [Greedy Pascal complements](schur-pascal-complement.md): an exact binary-linear-algebra formula for the infinity spectrum of every finite mode set, with a traversal that skips large index gaps and an all-length interval corollary.
- [Constant-flag crossings](constant-flag-crossings.md): polynomial Gaussian crossings from two arbitrary constant flags and a unitriangular endpoint; an exact four-mode, fixed-gauge screen covers all 4,096 lower flags over the four-element field.
- [Strict Gaussian endpoints and deformations](strict-gaussian-relabeling-crossings.md): a profile-preserving Pascal repair enables arbitrary affine relabelings, binary charge projections with exact blocks, and nonseparable entrywise Laurent shifts.
- [Compact polynomial frozen-row profiles](frozen-core-row-profiles.md): twelve compiled source matrices give exact order-64 profiles at every interval length, with two parameter orbits and recorded finite cancellation searches.
- [Slow corrections after ramification](slow-ramification-crossings.md): an explicit separated-degree threshold proves exact affine profiles, and a rational qutrit example exhibits a favorable finite correction at one scale.
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

- [Prime-power vacuum locality](prime-power-vacuum-locality.md): Artin--Hasse transcendence excludes the stated scalar finite-locality replacement, with independent finite coefficient controls and no factoring lower-bound claim.
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
