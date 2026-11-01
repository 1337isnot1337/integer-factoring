# Research index

Related checkpoint: [31 October ordinary-proof research session](frontier-session-20261031.md) — complete score unchanged at 1/2; exact constructors, global charge obstructions, completed finite searches, and final validation.

This index covers research notes on `main`. Open PRs can contain later work. Each note states its own assumptions, proof scope, and evidence. An obstruction to a specific construction is not a factoring lower bound.

<!-- beyond-half-index:start -->
## November 1 arithmetic and tensor checkpoint

[Overview and exact selector contract](beyond-half-research-checkpoint.md). Ordinary proofs and executable evidence; complete score remains 1/2.

- [Prove polynomial-size prime-form generation and exact biased samples](polylog-prime-form-generators.md).
- [Make prime-form generator and bias cutoffs explicit](explicit-prime-form-cutoffs.md).
- [Record prime-translation and smooth-set Fourier obstructions](prime-translates-miss-smooth-forms.md).
- [Smooth set fourier norm](smooth-set-fourier-norm.md).
- [Extract useful ambiguous forms in polynomial time from an annihilator](polynomial-ambiguous-extraction.md).
- [Decode primitive factored norms by exact orientation matching](factored-norm-orientation-search.md).
- [Prove explicit multiband smooth-ideal density](smooth-band-density.md).
- [Construct the prime-form base deterministically](prime-form-root-construction.md).
- [Amplify prime bias to obtain short smooth norms in every class](powered-prime-smooth-norms.md).
- [Replay powered-prime smooth-norm existence on finite class groups](powered-prime-finite-replay.md).
- [Complete the selected-integer factoring consumer](selected-norm-factoring-interface.md).
- [Remove the multiplicity promise from the unsigned consumer](height-only-selected-integers.md).
- [Decode arbitrary factored norms and minimize primitive content](arbitrary-factored-norm-representation.md).
- [Bound coverage of one target-independent norm list](norm-selector-coverage.md).
- [Prove the elementary crossing floor after all binary type restrictions](direct-type-elementary-pressure-floor.md).
- [Classify all rank-seven constant restrictions of three elementary sites](linear-physical-restriction-crossings.md).
- [Prove Clifford compression layers and near-full pressure floors](clifford-compression-pressure-floor.md).
- [Pin the scope of the recent rigorous ideal sampler](rigorous-ideal-sampling-scope.md).
- [Show the row-frame threshold is sharp and exhibit the cubic induction obstruction](clifford-threshold-sharpness.md).
- [Sharpen supplied-norm processing to a uniform one-half leading coefficient](quasipolynomial-base-norm-processing.md).
- [Direct supplied smooth norms give signed short prime words](smooth-norm-short-words.md).
- [Polynomial ideal HNF and uniform Delta=-4N extraction structure](direct-discriminant-short-norms.md).
- [Complete finite norm sieve and original 60/92-bit factoring certificates](complete-smooth-norm-sieve.md).
- [Determinant-free full and partial relation halving](determinant-free-ambiguous-extraction.md).
- [Finite partial norm search and 124-bit certificate](partial-norm-search.md).
- [Height-only norms and logarithm-weighted relation independence](height-only-smooth-norm-relations.md).
- [Adaptive target exponents improve the finite height search](adaptive-height-norm-selection.md).
- [Polynomial-size complete subspace Smith-channel oracle](subspace-smith-channel-oracle.md).
- [Exact all-state channel-pressure separation](channel-pressure-separation.md).
- [Exact block and precision hierarchy for complete subspace channels](block-subspace-channel-hierarchy.md).
- [Exact modular higher-spin crossing exclusions](modular-spin-crossings.md).
- [Exact homogeneous physical restrictions of the two-site reflection seed](reflection-homogeneous-restrictions.md).
- [All-scale 5/7 floor for every homogeneous middle-fiber restriction](homogeneous-restriction-pressure-floor.md).
- [Complete joint physical conjugations cannot improve the four-letter reflection seed](joint-height-conjugation-floor.md).
- [State-adaptive physical Smith channels](state-adaptive-smith-channels.md).
- [Matroid optimization within equal-order channel blocks](matroid-adaptive-smith-channels.md).
- [Complete regular channel bases through finite jets](dvr-adaptive-smith-channels.md).
- [Complete triangular auxiliary row fiber and word obstruction](triangular-auxiliary-row-fiber.md).
- [Protected subtensor and maximal regular left projection](protected-subtensor-pressure.md).
- [Bounded boundary extraction of protected scalar lines](boundary-extracted-scalar-pressure.md).
- [All binary code embeddings with arbitrary diagonal gauges](arbitrary-binary-code-pressure-floor.md).
- [Integer shared-minor floors and the universal affine exclusion](integer-shared-minor-pressure-floor.md).
- [All-degree seven-letter adjugate floor and capped jets](adjugate-seven-letter-pressure-floor.md).
- [All coordinate codes of two reflection sites have an exact pressure floor](reflection-coordinate-code-pressure-floor.md).
- [Smith orders bounded by degree close alphabets through ten](bounded-order-smith-pressure-floor.md).

<!-- beyond-half-index:end -->

## Proved machine bounds and their foundations

- [Half-rate factoring composition](half-rate-factoring.md): ordinary proof of the `1/2` machine bound; its final Lean theorem is in [`formal-appendix/HalfTheorem.lean`](../formal-appendix/HalfTheorem.lean).
- [Top constructor and codec](half-top-constructor.md): deterministic table construction and exact input translation.
- [Evaluator and encoding audit](half-evaluator-audit.md): recursive evaluator, computable pad, and charged bit costs.
- [Verified k6 machine](verified-k6.md): Lean proof of the earlier `1/30` bound.
- [Sharper k6 rate](k6-sharp-rate.md): ordinary rate derivation above `1/29` for the k6 program; this sharper exponent is not Lean-certified.
- [k6 source audit](k6-source-audit.md): hashes, modified source scope, and compiler preservation.
- [Baseline audit](baseline-audit.md): source recurrence and the unchanged reference's machine-cost analysis.

## Transfers and proposed stronger evaluators

- [Weighted fixed-characteristic transfer](fixed-characteristic-transfer.md): uniform factoring from a modular evaluator.

- [Weighted prefix constructor](weighted-prefix-constructor.md): exact weighted tables and their construction cost.

- [Modular crossing pieces](past-half-pieces.md): effective characteristic-two pieces and the half-rate bound.

- [Polynomial Pfaffian formula](modular-pfaffian.md): polynomial-time computation of complete piece profiles.

- [Shared-mode Schur reduction](fermi-sea-schur.md): exact common-block convolution and rank restrictions.

- [Sparse-grid geometry](past-half-sparse-grids.md): limitations of regular-grid density changes.

- [Coarsened grid check](coarsened-grid-validation.md): exact finite integration checks for signed network grouping, orientation, composite moduli, and padded clique tables.

- [Bare criterion through twelve letters](bare-seven-to-twelve-letters.md): proved box-minor Kronecker budgets couple the two cuts; exact entropy-dual certificates give rho_bare>1/sqrt(2) for 7<=q<=12 at every degree and field, with fixed auxiliaries; an exact thirteen-letter profile satisfies every box budget and adjugate inequality and is the remaining construction target.

- [Auxiliary binary gauge floor](auxiliary-binary-gauge-floor.md): the gauged binary floor extends to matrix-valued Laurent seeds of every fixed auxiliary dimension with independent rational auxiliary conjugation gauges at the two places, by dominating the rectangular Smith orders coordinatewise with a twisted scalar boundary contraction; an explicit seed shows that no constant contraction dominates at both places and that a boundary contraction `f^T K e` can have strictly larger pressure.

- [Cell-set minor caps](cell-minor-caps.md): every k-cell endpoint minor bounds each normalized Smith profile prefix by kd plus the k smallest cell energies; this exactly excludes the four-symbol shared-minor target, leaves certified sharper negative relaxation targets at five and six symbols, and proves a compact ternary region; the finite floating searches found no negative four-symbol value, but global four-symbol closure remains open. No realizing tensor or factoring improvement.

- [Joint determinant actual-rank barrier](joint-determinant-barrier.md): all-field Weyl filtrations close every binary pure-value quadratic and cubic at leading `rho>=4/5`; the weaker rank-two determinant floor is `2/3`, and the common-cofactor result above closes the ternary quadratic class.

- [Constant transpose-intertwiner barrier](transpose-intertwiner-barrier.md): equal actual cut ranks at every tensor length preclude contraction for polynomially many coefficient pieces; complete tests exclude 32,454,657 ternary cubic coefficient seeds, with actual-jet witnesses checked through three copies.

- [Universal multipoint Smith floor](multipoint-smith-floor.md): determinant budgets and equal local minimum orders force `rho > 3/4` for the unweighted local-jet criterion over every field, alphabet and degree; a complete 1,579,014-seed cubic catalogue provides independent controls.

- [Low-rank ternary quadratic multipoint catalogue](multipoint-low-rank-seeds.md): all binary seeds with cut-rank product at most three at both finite places have an exact `rho >= 3/4` pressure floor; nine stabilizer representatives, 32,265 seeds, 321 signatures, complete local-rank coverage and independent tensor controls.

- [Complete ternary quadratic multipoint catalogue](multipoint-rank-one-seeds.md): all binary seeds with pure values at both finite places have an exact `rho >= 3/4` pressure floor, with 3,078 representatives, 23 signatures, independent Toeplitz profiles and full tensor-power reconstructions.

- [All-code elementary coefficient-span barrier](difference-support-rank-floor.md): weighted difference/orientation blocks give rank-sum at least (|A||B|)^(3/4) for arbitrary paired binary codes over every field; equal alphabets force rho >= 3/4, without symmetry or downset assumptions.

- [Difference-unique code rank obstruction](sidon-code-rank-obstruction.md): a literal two-cut block proof excludes the elementary crossing on all cubic binary-field graphs; 16 exact finite controls distinguish the separate multiplication-graph candidate.

- [Dyadic inverse sections](dyadic-inverse-sections.md): exact arithmetic interval recursion, affine-fibre floor sums, integer-factor recovery, and the unresolved aggregation cost.
- [Dyadic orbit aggregation](dyadic-orbit-aggregation.md): polynomial bit time for complete fractional-linear cyclic families with interval selectors and arbitrary affine changes.
- [Interval orbit gcd](interval-orbit-gcd.md): correlated-root and rational-part interfaces, exact specialization caveat, output-sensitive algorithm target, and source audit; its stronger factoring consequence is conditional.
- [Interval unit jets](interval-unit-jets.md): polynomial-bit-time local Taylor data for exponentially long interval products using parity, power sums, and exact Newton identities; global elimination remains open.
- [Interval packet circuits](interval-packet-circuits.md): global modular polynomial circuits with `O(n+p)` nodes and polynomial generation cost; exact artificial-common-factor controls state the remaining elimination obligation.
- [Normalized interval jets](normalized-interval-jets.md): exact binary content exponents, zero orders, and polynomial-cost unit jets at scaled integer centers; global root search remains open.
- [Global interval pole trace](interval-pole-trace.md): an exact finite-precision factor-pair projection, closed boundary and bilinear section rules, and a charged enumerative baseline; efficient global state reduction remains open.
- [Cauchy determinant slope](pole-cauchy-content.md): a determinant slope or normalized trace from two Cauchy powers recovers an arbitrary rectangle's factor count; padding and Vandermonde normalization are cheap, while global scalar evaluation remains open.
- [Arithmetic framework shift](arithmetic-framework.md): primary-source randomness audit and a new generator/relation construction contract, outside the crossing evaluator.
- [Uniformity and the half frontier](uniformity-and-half-frontier.md): standard padding shows that `P=BPP` would beat every fixed saving exponent; small universal random-tape families exist nonuniformly, while efficient uniform construction remains open.
- [Deterministic torsion extraction](arithmetic-torsion-extraction.md): a black-box two-group basis algorithm from an annihilator, with no saturated relation lattice or random group elements.
- [Positive relation kernel](positive-relation-kernel.md): a Fourier-positive smooth-pair hitting lemma and diagonally dominant relation construction; efficient arithmetic pair extraction is open.
- [Biased words without coordinates](arithmetic-biased-words.md): a powering seed and published derandomized squaring supply biased samples from generators and a group-size upper bound.
- [Actual class-group controls](arithmetic-class-group-controls.md): exact composition, smooth prime-form words, unsaturated relations and complete two-torsion recovery; exhaustive searches produce factors for 18 finite targets.
- [Density/index transfer](arithmetic-density-lift.md): dense smooth forms and relations on a possibly proper factor-base subgroup give a global class-number multiple; missing ambiguous directions remain an explicit obligation.
- [Bounded smooth-norm interface](bounded-smooth-norm.md): exact equivalence between smooth-form quotients and balanced smooth represented norms; supplied witnesses convert at `2^O(b/log b)` cost, while fast witness search remains open.
- [Deterministic norm sieve](arithmetic-norm-sieve.md): pair and direct-word certificates factor 60/80-bit inputs without class enumeration; polynomial composition and bounded diagonals reduce determinant encodings to 520/500 bits, while failed boxes retain the search limits.
- [Quantitative smooth-norm localization](smooth-norm-localization.md): the exact biased generator-power family, an explicit unproved short-box hypothesis and a conditional transfer above half; splitting coverage remains a separate hypothesis.
- [Known-root conic transport](class-square-conic.md): an integral matrix of determinant |Delta| parametrizes every rational square-norm point; three explicit points solve every modular instance without factorization, while a negative integer example isolates the unresolved global denominator-divisor search.
- [Conic denominator audit](conic-denominator-descent.md): for discriminant -4N, cancellation is one explicit linear gcd and a parity factor; nonsplitting cases preserve the root class up to the exceptional ambiguous class, with exact failures of monotone descent.
- [Divisor moments and Eisenstein coefficients](divisor-moment-recovery.md): positive Hankel recovery gives a deterministic smallest-factor reduction to a uniform family of level-one coefficients; a trial/weight tradeoff makes the weights and recovery polynomial in b for target scores through one, while the fast coefficient engine remains open.
- [Global sparse divisor residues](sparse-divisor-residues.md): Garg--Schost interpolation and a chosen trial budget preserve every fixed-characteristic saving exponent, with polynomial cyclic degrees through score one and a global-only query contract; a faster evaluator remains open.
- [Fixed weight-one coefficient target](weight-one-divisor-coefficients.md): one bit per cyclotomic coordinate supplies signed divisor residues at every prime level, by an injective binary decoder retaining nilpotent coordinates; paired sparse recovery preserves every fixed saving exponent, conditional on a faster composite-index weight-one coefficient engine.
- [Polynomial divisor-gap extraction](smooth-norm-gap.md): integer multiplicative trimming converts promised inner-radius balanced norms at a specified outer radius, without enumerating all divisors.
- [Relations from small smooth norms](smooth-norm-relations.md): primitive smooth norms yield bounded prime-form words directly, removing the balanced-factor condition and strictly enlarging the accepted target set.
- [Exact search for a factored norm](factored-norm-search.md): ideal-lattice composition and deterministic meet-in-the-middle decide one primitive representation and recover coordinates in `Y poly(b)+2^O(b/log b)` bit time; selecting candidate norms remains open.
- [Reflected duals and a neutral crossing compiler](neutral-reflected-crossings.md): Lucas duality in every prime characteristic and a characteristic-two Pfaffian formula give complete balanced and unbalanced fixed-sector crossings.
- [All homogeneous dual completions](homogeneous-dual-completion.md): the complete right-inverse family at weight 32 has 120 dual parameters and exactly 88 independent tensor directions, with exact source and endpoint controls.
- [Crossing synthesis from arbitrary F4 frames](f4-recentered-crossing-lattices.md): full-entry affine coefficient search can recenter on nonbinary tensors; finite polynomial unit frames certify complete orders, with independent binary-jet replays.
- [Odd-cyclic matrix crossings](odd-cyclic-matrix-crossings.md): explicit matrix/Fourier cut normal forms, an all-degree unweighted `rho >= 3/4` floor, and a completed five-symbol weighted search with recorded numerical failures.
- [Bell-diagonal crossings](bell-diagonal-crossings.md): exact two-cut transforms, effective polynomial paths, and a complete symmetric-span restriction at local dimensions two and four; dimension eight remains open.
- [Bell linear-pencil search](bell-linear-search.md): exact coefficient ranks, complete small-domain calibrations, and retained higher-dimension finite schedules.
- [Bell polynomial synthesis](bell-polynomial-synthesis.md): exact affine two-cut construction and rational pressure bounds uniform in tensor length.
- [Folded-cube Bell family](folded-bell-crossings.md): explicit all-dimension crossing in every odd characteristic, with optimal coefficient-span limit `H_2(1/4)`.
- [Code-quotient crossings](code-quotient-crossings.md): exact projection and maximal common division, a characteristic-zero pressure restriction, and a precise modular-cancellation target.
- [Signed boundary crossings](signed-boundary-crossings.md): exact intertwiners and isotropic code reduction, with complete common-factor normalization and retained explicit paths.
- [Bell outer-rank endpoint floors](bell-cubic-outer-rank-floor.md): pressure exceeds `5/7` for the stated cubics at all local dimensions at least eight and for quartics at 8 and 16; the refined bound leaves 32 open.
- [Quartic Bell paths and matrix rigidity](bell-quartic-rigidity.md): a chirped reshape turns a self-reciprocal quartic into an explicit Hadamard sparse modification, with necessary critical parameters and a stated missing converse.

- [Joint residues reduce to finite templates and convex offsets](joint-syndrome-residue-retiming.md): at fixed normalized degree, arbitrary coefficient fields and ramification reduce exactly to finitely many binary endpoint covers and compact convex offset programs; a transport certificate bounds the entire offset box, with exact color and Hall controls.

- [Typed rational kernels have pressure-preserving strict forms](typed-rational-crossing-strictification.md): balanced valuation circuits yield paired potentials and positive-characteristic strictification; the scalar binary floor extends to all stated rational endpoints, with a sharp four-copy endpoint test and uniform eighth-power gauge in characteristic two.

- [Mixed forward/reverse paired quartics](reflection-mixed-paired-quartics.md): a common charge permits both directions on disjoint logical pairs, with asymmetric pole budgets M+7 and M-1; same-pair directions are allowed as an additive correction, while unproved overlapping products are rejected.

- [Coefficient/height alternation at exact integer chambers](coefficient-height-alternation.md): two completed thirty-minute budgets perform 61,041 target steps and produce 1,432 full-profile classes; all eighty saved tensors replay, with no gain, and the unfinished 50,000-step requests remain explicit.

- [All extreme charges are cube height functions](cube-charge-extreme-points.md): a connected-tight-graph proof classifies every charge vertex at every complete alphabet; exact catalogues give 1,980 vertices at sixteen symbols, and independent odd-parity signs give 2^(q/2) vertices with common endpoint pole at most three.

- [Amplified removable poles in growing quartic words](reflection-paired-quartic-poles.md): simultaneous Pfaffian row congruences permit M/4 disjoint gates with pole M+7 each, unchanged effective degree in the stated numerator window, and actual common pole at most M-1; an exact three-copy counterexample checks the first excluded boundary.

- [Exact nonlinear symbol-charge domain](nonlinear-symbol-charge-lattices.md): arbitrary physical-symbol potentials satisfy the full reflection coefficient inequalities exactly when their Hamming-cube Lipschitz constants sum to at most two; the four-symbol polytope has twelve vertices and sixty-four facets, including nonlinear parity charges.

- [Finite fused six-vertex crossings](finite-fused-six-vertex-crossings.md): literal three- and four-symbol F16 fusion, exact common-factor cancellation, 24 independent binary profiles, and 144 replayed full-gauge starts; the completed finite searches remain uncompetitive.

- [Exact F4 regular row fibers](colored-regular-row-fibers.md): nonbinary polynomial row operations preserve a common binary affine floor, allowing exact compressed F4 synthesis with complete coefficient windows and literal untruncated tensor replays.

- [Joint filters beyond the complete scalar multiplier family](joint-boolean-multiplier-separation.md): a literal three-parity polynomial filter beats every scalar Boolean support and arbitrary-degree multiplier on that quotient; 10,668 exact envelope certificates give a uniform separation at every balanced physical length at least sixteen, with positive full critical pressure.

- [Complete-entry typed rational coefficient lattices](typed-charge-lattice-synthesis.md): explicit entrywise endpoint valuations give full affine spaces with local poles and all-length polynomial type restrictions; four-symbol families have 6,144 or 6,237 coefficients at effective degree 24, with complete two-cut synthesis and independent coefficient controls.

- [Joint global-syndrome polynomial synthesis](joint-syndrome-polynomial-filters.md): independent syndrome polynomials retain all endpoint colors and full zero/infinity costs; exact residue compression removes ramification-size growth, and coefficient-space releases escape a demonstrated greedy-search trap.

- [Coherent physical codes with full letter gauges](physical-code-gauge-searches.md): disjoint odd fibers give exact coherent split maps; completed suites contain thirty coordinate, 158 coherent five-letter, and 57 compressed fourteen-letter seeds, with full type costs and independent contraction/profile controls.

- [Paired-type rational endpoints force full local rank](paired-type-endpoint-full-rank.md): crossing specializations of largest type classes force both scalar cut matrices to have full generic rank, including local poles and unequal alphabets; smaller-entropy hypotheses give the precise weaker rank bound.

- [Complete regular row-frame cells](formal-row-flag-cells.md): every formal frame with a fixed constant flag has a unique bounded lower polynomial representative; the degree-gap formula covers all jets over every field.

- [Finite-precision regular Smith frames](finite-precision-smith-frames.md): actual polynomial unit frames retain all coefficient equations below explicit row caps, avoiding irrelevant fraction-free tails; full products, independent profiles, and 316 Toeplitz ranks verify the implementation.

- [Fast matrix multiplication and the current recursion](fast-matrix-multiplication-recursion-audit.md): pinned source accounting isolates the child-volume term; changes to a fixed terminal solver or summable construction overhead do not change the certified logarithmic saving exponent, while a new multi-output primitive remains outside the audit.

- [Larger-alphabet shared-minor targets](shared-minor-profile-relaxation.md): all shared-line and rectangular-subalphabet prefix constraints admit explicitly certified negative relaxed pressure at four and five symbols; these are construction targets, not realizing tensors or factoring improvements.

- [Growing quartic words and removable endpoint poles](reflection-charged-quartic-tiles.md): ordered creation/annihilation gates leave the Gaussian class; an explicit charge budget proves exact pole cancellation and crossing endpoints after paired type restriction, with complete degree costs and original-tensor checks.

- [Complete one-support multiplier automata](boolean-multiplier-automata.md): arbitrary degree, characteristic-two coefficient extensions, and ramification reduce exactly to finitely many continuous cutoff functions and periodic binary tails; global code supports retain every actual physical mode energy.

- [Global code multiplier searches](global-code-multiplier-searches.md): 13,248 explicit supports on 64, 128, and 256 modes use the complete multiplier envelope; fifteen completed full-ratio searches select the baseline, with independent support and root replays.

- [Complete-entry larger-alphabet synthesis](full-alphabet-two-cut-affine-synthesis.md): all degree-bounded F4 coefficients satisfy exact endpoint-preserving two-cut equations; compressed fibers support two completed 20,000-walk four-symbol searches with 1,046 full-profile classes and eighty complete witness replays.

- [Nonresonant filters for arbitrary Boolean supports](nonresonant-boolean-filter-floor.md): the complete pressure floor extends to every Boolean pattern in a two-power family outside a finite resonant window; an exact counterexample preserves the excluded resonance.

- [Affine Boolean filter tiles](affine-boolean-filter-tiles.md): literal coefficient relabeling and exact monomial normalization leave the symmetric family; all 124,080 specified maps and witnesses regenerate 65,475 retained full-state profiles, with eight complete root replays.

- [Coherent physical split projections](physical-split-crossings.md): constant left/right splittings preserve strict endpoints and permit unequal alphabets; full profiles and actual new-symbol heights are recomputed, with 5,480 completed finite projections and independent F4 restriction-of-scalars controls.

- [Reflection coordinate-code alphabets](reflection-code-alphabets.md): non-Cartesian physical state restrictions leave the Gaussian class while keeping the endpoint; a 236-alphabet natural-gauge screen and exact original-tensor checks retain the full physical entropy cost.

- [Growing Gaussian physical frames](growing-gaussian-physical-frames.md): endpoint-trivial creation and annihilation words may change a growing fraction of modes, with arbitrary independent mode heights, complete degree costs, and full original-tensor controls.

- [Coupled Gaussian gates](coupled-gaussian-gates.md): cross-letter quadratic operations have an exact rank-two Clifford update with a determinant-mean correction; all rational denominator costs and original-Fock profiles are retained.

- [Symmetric filters in every prime characteristic](symmetric-prime-multilevel-filters.md): Vandermonde compression retains complete digit multiplicities, endpoint signs, exact polynomial feasibility, and tensor costs at actual physical energies.

- [Completed symmetric tile catalogues](symmetric-tiled-filter-searches.md): 3,830 exact polynomial profiles and 32 full positive-root searches across binary and qutrit intervals, with every saved witness reconstructed and all best points selecting the baseline.

- [Universal binary Smith-pressure floor](universal-binary-smith-floor.md): every Laurent-polynomial binary crossing seed, every coefficient field, arbitrary degree, all regular frames, and all independent physical-letter heights have strictly positive critical pressure; shared minors and an exact 421-leaf compact certificate close this entire two-place Smith route.

- [Symmetric multilevel filter synthesis](symmetric-multilevel-filter-synthesis.md): occupation-weight compression and exact tensorization preserve full degree costs and actual physical energies; a literal four-mode filter gives strict natural-charge relative gains without a better full crossing root.

- [Exact Gaussian fixed-type duals](reflection-gaussian-fixed-type-duals.md): rational active-minor certificates prove global minima over all independent gauges and scales for explicit physical marginals at intervals 4, 8, 16, and 32, without claiming an all-type obstruction.

- [Multipoint Hermite crossing transfer](multipoint-hermite-crossings.md): local jets at several points share the degree budget, reconstruct exact crossing pieces, and admit independent two-cut scales with complete factor-construction costs; an exact ternary example validates the criterion below the existing performance level.

- [Global homogeneous one-level filter floor](global-homogeneous-boolean-filter-floor.md): all overlapping homogeneous Boolean occupation filters have positive complete pressure at every charge and scale, with unrestricted admissible degree and an exact compact cover for the remaining nine lengths.

- [Complete Gaussian profiles from Clifford matrices](clifford-gaussian-smith-profiles.md): an elementary characteristic-two Cartan argument recovers all Fock Smith orders from a matrix of twice the mode count, with closed reflection matrices and no internal regularity restriction on physical mode gauges.

- [Nonconstant regular binary flags](nonconstant-binary-row-flags.md): arbitrary regular polynomial row frames keep complete two-cut synthesis affine, finite jets suffice for fixed targets, and the full first-jet F4 catalogue retains only the reference line.

- [Terminal operators above a power of two](above-power-identity-operator-reduction.md): an existing unfiltered low identity direction closes a small square Weyl core while keeping its physical spectator and degree accounting explicit.

- [Completed above-power operator catalogue](above-power-operator-searches.md): all 42,560 specified interval-33 projector kernels yield 1,111 exact profile classes, with complete numerical grid and refinement records.

- [Prime-field affine operator synthesis](prime-deficit-affine-operator-synthesis.md): all occupation-preserving qutrit coefficient entries become exact affine unknowns, with strict identity endpoints, explicit source-precision guards, constant infinity flags, and four completed finite path searches.

- [Operators on arbitrary one-deficit selections](arbitrary-one-deficit-operator-search.md): a compatible selected-index API with full original controls and four complete parent-exponent qutrit catalogues.

- [Energy-conjugated one-deficit operators](deficit-energy-conjugated-operators.md): selected-index monomial conjugation preserves the strict endpoint and charges its complete Laurent support; four further qutrit catalogues completed.

- [Arbitrary selections one below a prime power](arbitrary-one-deficit-identity-cores.md): weighted principal minors retain only 2r+1 actual modes for every selection, with one explicit local eigenchannel and its endpoint restriction.

- [Arbitrary selections above a prime power](surplus-identity-cores.md): adding unchanged actual identity directions produces a square core with at most 2(r+2d) modes in every prime characteristic; terminal selections use 2(r+d).

- [Fixed deficits below prime powers](fixed-deficit-identity-cores.md): mixed initial and terminal directions have a uniformly bounded rectangular core, exact Gaussian-digit division, and a rank-at-most-deficit contraction correction in every prime characteristic.

- [Complete-entry binary crossing synthesis](binary-two-cut-affine-synthesis.md): exact affine systems at both cuts cover free Laurent entries, a support lemma bounds the required degree for specified targets, and complete constant F4 flag searches decide those targets.

- [Joint residue-compressed polynomial synthesis](joint-binary-filter-synthesis.md): at most twice the number of digit states plus one residue patterns give exact degree/order feasibility in every prime characteristic, with a constructive fixed-field endpoint lift.

- [Completed joint-filter catalogues](joint-filter-search-records.md): every retained literal witness is reconstructed, and finite all-index grids keep mode-assignment coverage and interruption limits explicit.

- [Polynomial filter flags and zero-unit removal](weyl-filter-flags-and-unit-reduction.md): constant extension-field flags and exact polynomial unit division preserve valid row-lattice bounds and reduce large non-prime-power sources.

- [Identity operators below a power of two](near-power-identity-operator-reduction.md): weighted interpolation retains one extra auxiliary mode, exact Gaussian multiplicity division enables a small Weyl evaluator, and the contraction has an explicit scalar-plus-rank-one formula.

- [Completed all-order bulk searches](reflection-hot-all-bulk-searches.md): all sixteen ordering cones, fixed and interval-scaled height boxes, complete numerical records, and replay of the saved best points.

- [Completed separate-leg bulk searches](reflection-hot-leg-bulk-searches.md): both 25-candidate searches at intervals 64 and 128 completed with full cold and physical costs; the resumed boundary is recorded.

- [All physical-leg bulk orderings](reflection-hot-all-bulk-gauges.md): two explicit regular cross pivots cover every nonnegative bulk ordering with the same bounded support, complete coefficient sources, and independent original-tensor controls.

- [Separate physical-leg bulk heights](reflection-hot-leg-bulk-gauges.md): explicit ordering and hot-cold inequalities preserve a bounded reflection core, with distinct physical type costs and exact original-tensor controls.

- [Completed cold-compensated searches](reflection-hot-compensated-searches.md): full coefficient grids on negative-height rays, scoped continuous searches through interval 64, and nonlinear OR-state penalties include every compensating cold cost.

- [Completed reflection hot-state screens](reflection-hot-state-searches.md): complete positive-height coefficient grids, scoped unrestricted heights, and nonlinear two-bit penalties retain full type costs and report no improvement within their recorded scopes.

- [Four-mode resonant occupation transfers](binary-four-mode-resonance.md): complementary equal-energy states reduce to five scalar factors and a two-by-two block, with exact profiles and an all-scale obstruction for polynomial blocks of degree at most twice the interval scale.

- [Complete Laurent resonance catalogues](resonant-laurent-catalogues.md): exact multiplicities and profile signatures cover every invertible-endpoint F4 matrix in two specified three-power windows; numerical searches over all resonant interval-32 quadruples found no gain.

- [Operator lattices over fixed extension fields](identity-operator-field-flags.md): nonbinary row bases and F16 endpoint lifts combine with a small scalar-Pascal source formula for every selected subset.

- [One-hole minimum-order tradeoff](identity-operator-minimum-tradeoff.md): with a monomial top coefficient, any positive infinity-minimum gain forces a quantified zero-minimum loss over every fixed coefficient extension and arbitrary selected indices.

- [All identity directions at prime-power lengths](prime-power-joint-identity-filters.md): a single regular reflected-Pascal change separates arbitrary selected indices, gives exact joint finite-difference profiles, and reuses one coefficient catalogue across every selected subset.

- [Negative hot-state heights with cold-bit compensation](reflection-hot-compensated-gauges.md): retaining only the internal hot-core factors permits negative selected-state heights, with every compensating cold height and both independent boundary bulks included in the exact cost.

- [Exact operator profiles over F16](f16-identity-operator-profiles.md): restriction of scalars duplicates each Smith order exactly twice, enabling complete extension-field profiles with literal original-tensor and all-minor controls.

- [Fixed-field endpoint lifts for polynomial lattices](extension-field-lattice-lifts.md): exact packed `F4` constraints admit constructive `F16` lifts whenever none of the eight required endpoint coordinates is forced to zero, with a checked equality-case fallback.

- [Two-cut synthesis for occupation-transfer kernels](identity-operator-lattice-synthesis.md): all operator coefficients are solved together under independent zero and infinity row constraints, with explicit Laurent support and an exact repaired-endpoint test.

- [Physical state gauges on reflection hot cores](reflection-hot-state-gauges.md): a regular reduction supports arbitrary state heights and endpoint-preserving masks on at most eleven modes per hot bit; monotone exterior-defect heights require only four.

- [Polynomial synthesis through reduced Smith lattices](weyl-lattice-filter-synthesis.md): linear coefficient constraints certify every occupation-block Smith prefix, with complete infinity costs and an exact eight-state endpoint-cover test over `F4`.

- [Mixed initial and terminal identity modes](split-identity-mode-reduction.md): a weighted interpolation proof concentrates any split of selected endpoint modes into a bounded core, including three-mode operator kernels and exact original-tensor controls.

- [Uniform eight-parameter hot-core comparison](reflection-hot-core-finite-grids.md): exact enumeration of all 65,536 F4 coefficient rows and a bivariate minor-degree argument prove no natural-gauge gain at any scale or any power-of-two interval.

- [Nonlinear three-direction filters over F4](f4-weyl-filter-screen.md): complete polynomial compilation and a portable catalogue cover 369,603 colored monomial filters, with every distinct exact profile retained and independently replayable.

- [General extreme-core transition masks](raw-core-transition-masks.md): arbitrary binary or `F4` masks preserving every diagonal transition keep the typed endpoint and admit exact charge-block profiles with arbitrary physical core gauges.

- [State filters with nonzero endpoint repair](repaired-endpoint-state-filters.md): a free diagonal physical repair enlarges scalar and joint polynomial synthesis, including a strict feasibility example, signed prime-field filters, and colored binary filters over `F4`.

- [Weighted Weyl reduction of coupled operators](identity-operator-weyl-reduction.md): an exact differential-operator identity reduces the active matrix order from `p^(2r)` to `p^r`, with explicit auxiliary shifts, weighted regularity, and complete polynomial coefficient arithmetic.

- [Bounded hot cores after reflection compression](reflection-hot-core-defects.md): arbitrary even polynomial defects on the hottest logical modes reduce to three active modes per selected mode, with complete Laurent coefficients and all Gaussian spectators retained.

- [Coupled occupation-transfer operators](coupled-identity-operator-filters.md): exact bounded active matrices extend the joint diagonal filters to occupation-preserving kernels, with polynomial degree costs and a larger energy-filtered endpoint class.

- [Coupled polynomial filter synthesis](coupled-filter-linear-synthesis.md): exact linear constraints certify a complete zero-partition upper bound while charging all infinity degrees, with brute-force feasibility controls and scoped finite searches.

- [Reflection-pair Gaussian compression](reflection-dual-rail-compression.md): an explicit characteristic-two homology calculation reduces the physical alphabet from `2^M` to `2^(M/2)`, with a polynomial endpoint, complete two-cut profiles, and the identical natural rate at twice the scale.

- [Complete finite frozen-row searches](frozen-core-complete-search.md): exact coefficient flags, affine cell counts, and rational-scale profile comparisons cover all 4,194,304 rows in each of 86 stated gauge/scale/orbit cases, with no same-gauge Gaussian improvement.

- [Strict fixed-core charge projections](strict-core-charge-projection.md): regular support on at most four times the physical core size, exact physical-core gauges, and complete fourteen-class screens at intervals eight and sixteen.
- [Raw extreme-core charge projections](raw-extreme-charge-projections.md): endpoint-preserving masks on both extreme modes retain all physical states, with ten active fermion modes and exact original weighted-cut controls.

- [Exact prime-state filter searches](prime-state-filter-searches.md): portable monomial, dual, cancellation-branch, and polynomial synthesis searches with exact profiles and complete degree costs; the recorded finite grids give no threshold improvement.

- [Coupled identity-filter reduction](coupled-identity-filter-reduction.md): regular paired shears concentrate arbitrary correlated digit weights into twice as many active zero modes, with terminal and opposite-end extensions at prime-power interval lengths.

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
- [Accelerated exact interval certificates](accelerated-interval-certificates.md): a positive Euler-series remainder removes the scan through modes and certifies an interval of 2^1636-1 modes, with a rational upper bound less than 2e-980 above the unchanged threshold.
- [Greedy Pascal complements](schur-pascal-complement.md): an exact binary-linear-algebra formula for the infinity spectrum of every finite mode set, with a traversal that skips large index gaps and an all-length interval corollary.
- [Constant-flag crossings](constant-flag-crossings.md): polynomial Gaussian crossings from two arbitrary constant flags and a unitriangular endpoint; an exact four-mode, fixed-gauge screen covers all 4,096 lower flags over the four-element field.
- [Strict Gaussian endpoints and deformations](strict-gaussian-relabeling-crossings.md): a profile-preserving Pascal repair enables arbitrary affine relabelings, binary charge projections with exact blocks, and nonseparable entrywise Laurent shifts.
- [Compact polynomial frozen-row profiles](frozen-core-row-profiles.md): twelve compiled source matrices give exact order-64 profiles at every interval length, with two parameter orbits and recorded finite cancellation searches.
- [Slow corrections after ramification](slow-ramification-crossings.md): an explicit separated-degree threshold proves exact affine profiles, and a rational qutrit example exhibits a favorable finite correction at one scale.
- [Prime-state filter diagonalization](prime-identity-filter-diagonalization.md): explicit polynomial Pascal factors give all active orders from finite differences of arbitrary Laurent digit weights, uniformly in every prime characteristic and interval length.
- [Qutrit identity-mode reduction](qutrit-identity-zero-reduction.md): regular concentration reduces an all-interval perturbation to two active modes, with exact charge-block factorizations.
- [Qutrit state-filter comparison](truncated-symmetric-identity-filters.md): polynomial endpoint-preserving filters have complete profiles, and every nonidentity monomial filter increases the Gaussian product in the stated charge range.
- [Exact polynomial filter synthesis](polynomial-state-filter-synthesis.md): degree bounds, zero orders, and the crossing endpoint reduce to at most `2p+1` residue systems; explicit witnesses retain every degree and denominator cost.
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

- [Nine-letter quadratic-inverse probe](nine-letter-quadratic-inverse-probe.md): a general first-two-coefficient tradeoff and complete 31-path schedule, with every exact profile, 248 independent Toeplitz controls and all-scale exclusions.
- [Inverse-degree compatibility and noncommuting paths](inverse-degree-crossing-compatibility.md): a new cut constraint rules out the formal seven-letter optimum at all integer ramifications; a quadratic-inverse constructor supplies exact endpoints and a complete 5,005-path first-jet exclusion.
- [Pieces-to-factoring transfer](pieces-to-factoring.md): conditional crossing-parameter transfer, with full bit-cost accounting.
- [Scale-one modular transfer](scale-one-modular-transfer.md): fixed-field transfer and machine costs.
- [Sparse divisor nonzeroness](sparse-divisor-nonzeroness.md): detection of divisors in weighted polynomial counts.
- [Modular crossing](modular-crossing.md): normalized finite-field pieces and quotient construction.
- [Modular parameter descent](modular-parameter-descent.md): the state ideal also annihilates the operator parameter, giving a direct reduced-dual compiler and a simplicity restriction on further quotients; the pieces and factoring score are unchanged.
- [Homogeneous duals](homogeneous-duals.md): alternative algebraic pieces and their checked identities.
- [Near-ones label families](near-ones-labels.md): label construction and rank experiments.
- [Binary-prefix top channel](binary-prefix-top-channel.md): structural analysis of channel widths.
- [Primary central Gram matrix](primary-central-gram.md): source-linked rank calculations.
- [Half-piece ranks](half-piece-rank.md): finite crossing-piece rank evidence.

- [Qutrit operators below a power of three](qutrit-deficit-operator-search.md): a complete prime-field Weyl evaluator with exact dummy-factor removal, original-tensor controls, and two complete 1,129-kernel projector searches at interval 26.

## Obstructions and work beyond 1/2

- [Regular cyclic-algebra transforms](cyclic-algebra-crossings.md): polynomial changes preserve ideal column orders, with exact dense-operator controls; 1,024 eight-letter transforms are excluded for the bare two-place criterion.

- [Locally joined shift registers](shift-register-crossings.md): full-cycle seeds from maximal feedback rules, exact linear cut-rank formulas, and complete bare-profile exclusions through thirty-two physical letters.

- [Copied-wire register barrier](shift-register-copy-wire-barrier.md): a colour-correlation argument forces limiting bare parameter at least `3/4` for every prime-power digit full-cycle register, including arbitrary nonlinear feedback; a full-cycle counterexample separates two-power and three-power rank inequalities.

- [Complete six-bit nonlinear catalogue](quadratic-register-crossings.md): all 65,536 quadratic-plus-zero-indicator feedback rules, their 1,952 full cycles, and all 19,520 physical partitions fail the bare criterion by exact certificates.

- [Quantum-Pascal trace crossings](quantum-pascal-trace-crossings.md): divided-power polynomial endpoints through fifteen letters, an exact degree formula, and a unitriangular trace-basis obstruction at `rho>=3/4` for bare and unweighted multipoint profiles; physical word types remain separate.

- [Convex fixed-type Smith capacity](fixed-type-smith-capacity.md): global entropy-weighted minor-flag supports, a rational dual-certificate criterion, exact convexity controls, and a bounded numerical paired-type solver with separate primal and LP estimates.

- [Matrix-vacuum boundary barrier](matrix-vacuum-boundary-barrier.md): arbitrary polynomial matrix vacua retain a square-basis determinant; an all-field boundary difference bounds the zero-cut rank, and the complete 7,350-case binary rank-two affine catalogue forces bare parameter above `7/8`.

- [Steiner cubic crossings](steiner-cubic-crossings.md): three controlled involutions give determinant-one crossing paths for every Steiner quasigroup; exact design-independent profiles exclude the bare criterion at every size, with a six-dimensional orbit calculation and full-cut controls.

- [Steiner coefficient barrier](steiner-coefficient-barrier.md): shared coefficients force actual ranks at every word length; the complete multivariate span has `rho>=3/4` for `q>=7`, and the three-letter univariate span has asymptotic floor `0.8065306485`, also with a sublinear number of independent position classes.

- [Ordered physical-gauge crossings](ordered-physical-gauge-crossings.md): noncommuting matrix gauges preserve every endpoint and the polynomial coefficient count; exact rank-sum dual potentials exclude 216 stated binary controls, with a distinct clock/shift phase construction left for further work.

- [Ordered Weyl boundaries](weyl-ordered-boundaries.md): explicit quadratic-phase word kernels and an all-length rank-one cold-boundary classification; exact four-jet certificates reject all 4,096 normalized F4 boundaries at two positions, and complete middle-coefficient tests reject actual F4 spans at two selected boundaries; longer ordered words remain open.

- [Subgroup-return crossings](subgroup-return-crossings.md): group-word membership preserves the crossing and a rank-one starting cut at every length, including an effective free-group filter with infinite Hankel rank; complete binary coefficient spans at lengths one through three reject two examples without a uniform no-go claim.

- [Ordered charge-area crossings](ordered-charge-area-crossings.md): growing clocks have an exact integer charge-walk normal form; polynomial-degree field descent preserves the exponent, and a binary area-zero filter retains both strict endpoints with linearly many coefficients; complete binary spans fail at two selected lengths.

- [Joint energy conjugations](joint-energy-conjugations.md): joint-pair changes preserve one entire weighted Smith profile; a 16-dimensional bimodule enumerates the complete 1,296-member F4 orbit fixing the five-letter reference's leading coefficient, and exact prefix dominance excludes any scale improvement at its natural heights; a tangent identity classifies all formal leading-rank-one motions near the identity.

- [Polynomial joint-energy dressings](joint-energy-polynomial-dressings.md): a four-dimensional square-zero sector preserves one full profile; 180 quadratic identities and two explicit minors reduce every degree over every characteristic-two field to a constant conjugate under the original pole/rank constraints; exact quadratic/Frobenius elimination also closes every square-zero scalar-polynomial energy generator of arbitrary rank.

- [Group-shear crossings](group-shear-crossings.md): multiplication/division shears give cubic swap paths over every finite `p`-group in characteristic `p`; exact early-jet certificates exclude six displayed binary seeds, without an all-group claim.

- [All-degree homogeneous Boolean charge bands](homogeneous-boolean-charge-bands.md): a degree-independent loss cap and stronger Gaussian floors exclude all exterior charges and a linearly growing interior band, including overlapping monomials.

- [All-digit charged Gaussian floor](all-digit-charged-gaussian-floor.md): strict monotonicity in every real digit parameter and positive full critical pressure for all `q>1`, charges, scales, and interval lengths; exact compact certificates plus analytic tails.

- [Read-once Boolean filter floor](read-once-boolean-filter-floor.md): a degree-controlled log-odds theorem excludes every homogeneous read-once filter at interior charges and every parity filter at all charges, without a selected-mode bound.

- [Multilevel canonical-packet collisions](multilevel-thue-morse-packet.md): explicit one-bit Thue--Morse inputs collide with zero on every retained valuation level below a chosen cutoff; exact local-polynomial controls supplement the all-size proof.

- [Global charged Gaussian and one-mode Laurent floor](global-charged-gaussian-floor.md): the binary Gaussian has positive critical pressure at every integer length, real aggregate charge, and scale; the complete admissible one-mode Laurent family is also excluded at all charges.

- [Continuous charge folding and the wider tail floor](folded-charged-gaussian-floor.md): exact all-digit folding and a convex-tail prefix comparison extend the binary continuous floor to absolute charge at most 2n-2-2sqrt(2).

- [Small charged Gaussian bases](small-charged-gaussian-bases.md): analytic tails and 1,627 complete integer-outward rectangles prove all-scale residual bases at lengths one, two, and three.

- [Exterior Gaussian charge floor](exterior-charged-gaussian-floor.md): exact ground cancellation and convexity prove positive full pressure for every digit count at every real charge with absolute value at least twice the interval length.

- [Charged one-mode Laurent filters](charged-one-mode-binary-laurent-filters.md): a coefficient barrier proves Gaussian domination throughout the full interior charge interval, while a sharp universal loss cap and 70-interval certificate exclude every exterior critical witness.

- [Charged Gaussian pressure floor](charged-gaussian-pressure-floor.md): an integer interval certificate and analytic tails prove all-scale positivity on a linearly growing real-charge band; an exact identity treats every even charge in the physical range for all digit counts.

- [All one-mode binary Laurent filters](one-mode-binary-laurent-filter-floor.md): a complete degree/order dichotomy proves Gaussian domination at zero charge for arbitrary degree and every characteristic-two coefficient field, with exact equality classification.

- [Homogeneous Boolean occupation filters](homogeneous-boolean-filter-bounds.md): an exact two-probability cost formula proves all-scale natural-charge domination, a uniform large-scale comparison at every charge, and a bounded-charge small-scale critical-pressure bound.

- [Prime-power vacuum locality](prime-power-vacuum-locality.md): Artin--Hasse transcendence excludes the stated scalar finite-locality replacement, with independent finite coefficient controls and no factoring lower-bound claim.
- [One-sided charge pruning](one-sided-charge-pruning.md): common-type restriction removes every strictly one-directional charge perturbation; both-cut Smith majorization gives a matching untyped obstruction.
- [Typed Gaussian generic blocks](typed-gaussian-generic-blocks.md): arbitrary coefficient recombinations have cost `Omega(D^2)` at fixed nonvacuum occupations, with exact finite controls and a separate growing-density bound.
- [Exact typed charge-block compression](typed-gaussian-charge-blocks.md): polynomially many orbit blocks compute every coefficient-combination rank, including complete binary searches through 1,260 physical labels and independent dense controls.
- [Sparse bulk species-gauge floor](sparse-bulk-species-gauge-floor.md): an exact entropy dual proves `rho >= 3/4` in the stated common-density limit model; no exclusion of general sparse families or unbounded energy slopes is claimed.
- [Controlled monomial Smith floor](controlled-permutation-smith-floor.md): a fiber-span argument extends the nilpotent-path obstruction to dependent, noncommuting, and weighted monomial controls.
- [Qutrit Weyl crossing](qutrit-weyl-crossing.md): a complete auxiliary crossing with an exact all-length constant-coefficient rank formula and a scoped raw-piece obstruction.
- [Inverse-degree aspect bounds](inverse-degree-aspect-bounds.md): entropy slack gives degree/inverse-degree limits, stronger quadratic caps 14,20,35,91, and an abstract two-case exclusion for all 2^31 binary level weights on the fixed nine-letter flag.
- [Leading-rank entropy floor](leading-rank-entropy-floor.md): rational entropy duals exclude a repeated minimum Smith order on 7 through 10 letters, and close every positive kernel-level weighting of the displayed seven-letter transpose channel.
- [Shared-minor floor through six letters](shared-minor-smith-floor.md): common q-by-q minors couple both cuts and give an all-degree, all-field bare-profile exclusion for q<=6, including fixed auxiliary states; seven letters are the first case left open by this relaxation.
- [All-degree small-alphabet Smith floor](all-degree-small-alphabet-smith-floor.md): over every field and arbitrary fixed auxiliaries, bare ungauged local profiles on 2, 3 or 4 physical labels cannot improve the half-rate threshold; a sharper scalar binary parity floor holds at every degree.
- [Complete binary quadratic catalogue](binary-quadratic-catalogue.md): all 65,536 fixed-endpoint scalar binary quadratics, 162 paired profiles, and an exact all-scale pressure floor of 256/27 for the stated ungauged cuts.
- [Extended binary polynomial catalogues](extended-binary-polynomial-catalogues.md): all nonzero rank-one column constants in quadratics, complete parity cubics and quartics, and a separately marked million-path unrestricted cubic sample; no improved pressure.

- [Directions beyond half](post-half-directions.md): experiments, scope, and unresolved questions.
- [Network tradeoff](half-top-structure.md): restriction for the current local-enumeration and wire-channel constructor.
- [Restricted-label counterexample](half-restricted-labels.md): failure of a proposed low-mode dimension reduction.
- [Restricted intermediate spans](half-restricted-span.md): quantitative obstruction, with precise selected-label conditions.
- [Half-rate recurrence barrier](half-recurrence-barrier.md): limits of the existing operators and recurrence.
- [Graded-space barrier](graded-space-barrier.md): lower bound under the stated graded-space assumptions.
- [Rectangular-grid frontier](rectangular-grid-frontier.md): limits of unequal splits under the declared-profile recurrence.
- [Top-channel rank](top-channel-rank.md): exact factor rank on the full partition domain.
- [Coefficient basis](coefficient-basis.md): scoped minimax calculation for the chosen coefficient span.

## Binary actual-rank tools and physical codes

- [Coherent bounded-spin crossing](bounded-spin-crossing.md): uniform height-three path encoders over F5, Fibonacci growth, strict endpoint, exact binary widths and a finite-state entry constructor; a post-half actual-rank bound remains open.
- [Bounded-spin primary rank algebra](bounded-spin-rank-algebra.md): exact six-dimensional cut algebras at length eight, a rank formula over every characteristic-five extension, and a sharp minmax product 2304 certified by 759 common kernels and a separating hyperplane.
- [Second coefficient and the transposition class sum](bounded-spin-second-coefficient.md): a uniform compressed-operator formula identifies the eight-site row algebra with the replica swap and diagonal transposition sum, including its square-zero component; larger-length ranks remain open.
- [Three uniform spin-zero jets](bounded-spin-three-jets.md): integral scalar one-site transition matrices force three coefficient relations and an `n-2` span bound at every even length at least six; a four-qubit invariant model separates this from the open rank-growth problem.
- [Four-qubit generators in odd characteristic](four-qubit-orthogonal-generators.md): a direct `SO_4 x SO_4` proof of the free degree `2,4,4,6` invariant ring, complete F5 pencil relations and the exact universal `n-2` scalar Hilbert function; ordinary cut ranks remain open.
- [Enlarged coherent six-vertex crossing](bounded-spin-six-vertex.md): a uniform degree-`O(n^2)` substitution exposes the full monomial span and adds a seventh direction at eight sites; each cut has an exact invariant commutator module of dimension 28, with uniform favorable ranks still open.
- [Spin-family modular interface](bounded-spin-modular-interface.md): the existing fixed-characteristic transfer applies without rational rank lifts; dense factorization has the permitted cost, while automatic coefficient selection is proved for the original linear-size span and remains an obligation for arbitrary enlarged-space recombinations.
- [Constant-state spin projectors and cut reciprocity](bounded-spin-projector.md): exact thirty-state full and sixty-state binary-selected projectors preserve both ambient cut ranks; swapping the two balanced-product exponents identifies the normalized cut operators at every length.
- [Bounded physical support and uniform rank factors](bounded-spin-support-ranks.md): the shared extreme tensor has explicit channels growing as `(2+sqrt(2))^(n/2)`; near-corner coefficients have uniform combinatorial factors, while sublinear corner bands alone cannot span the crossing.
- [Binary site-symmetry actual-rank barrier](binary-site-symmetry-barrier.md): a universal degree-four certificate forces `rho>=3/4` over any fixed field, including all scalar polynomial degrees, coefficient recombinations, and common partitions into sublinearly many position classes.
- [Full auxiliary coefficient-span optimization](auxiliary-coefficient-span.md): all 2,588 binary functionals in 15 finite cases, strict recombination gains, exact minmax certificates, and rational rank-sum optima certified by shortest-path dual potentials.
- [Nonlinear downset-code crossings](downset-code-crossings.md): exact global characteristic-two Smith factors and a constant basis normal form for actual coefficient ranks on every downset and every paired-downset encoding, with complete small-domain controls and nonlinear regular-code searches.
- [Matching-value curve controls](matching-curve-controls.md): degree-five Nullstellensatz certificate rules out common projected cofactors over all extensions; a degree-eleven inverse certificate proves disjoint cofactor spans. Both remain subject to the broader site-symmetry barrier.
- [Common projected-cofactor actual-rank barrier](common-cofactor-barrier.md): exact certificates for every one of 32,265 normalized ternary rank-product-two quadratics cover all 12,852,225 coefficient seeds and force `rho>=log(3)/(2log(2))` at every tensor length.

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

## Retained PR #23 time-space and crossing checkpoint

These ordinary proofs and finite controls were integrated from [exact PR #23 source](https://github.com/1337isnot1337/integer-factoring/tree/147652574e7ce1ff68b2c378712846cbf438da34). Integration is not independent mathematical acceptance. The complete time exponent remains `1/2`.

- [Crossings from partitions into blocks of total weight zero](balanced-partition-crossings.md)
- [The original constant-weight three-label restriction](balanced-three-label-barrier.md)
- [Exact ranks of the elementary binary coefficient family](binary-coefficient-entropy.md)
- [Recover witnesses with bit tags and deterministic isolation](bit-tagged-witness-recovery.md)
- [A half-rate family with nearly linear-time factor queries](blocked-half-family.md)
- [Conserved convolution crossings and a two-cut obstruction](conserved-convolution-crossings.md)
- [A rank obstruction for coordinate matching support](coordinate-support-rank.md)
- [Odd cyclic orbit sums of a finite crossing seed](cyclic-orbit-crossing.md)
- [Direct exact divisor counts in sublinear workspace](direct-sublinear-divisor-counts.md)
- [A sharp rank-sum budget with one distinguished piece](distinguished-piece-rank-budget.md)
- [Dual gauges and exact mode filters](dual-gauge-filters.md)
- [Faster baselines for the recursive clique evaluator](fast-matrix-clique.md)
- [Reflected monomial duals and a four-piece crossing certificate](fermionic-piece-polynomial.md)
- [Exact crossings from finite-order algebraic approximations](finite-order-crossing-transfer.md)
- [Amplifying a finite crossing certificate](finite-seed-amplification.md)
- [A rank obstruction for permutation-invariant defect families](fixed-defect-invariant-barrier.md)
- [Fourier queries for every original multiplication gate](fourier-multiplication-gates.md)
- [Compact arithmetic-gate Fourier queries in odd characteristic](fourier-odd-characteristic.md)
- [The ideal Hankel row profile cannot have a rank-one opposite endpoint](ideal-hankel-endpoint.md)
- [Exact divisor counts in polynomial space](integral-counting-transfer.md)
- [Integral Schur crossing families and prime-power block composition](integral-schur-crossing.md)
- [Matching crossings with monomial duals](intersection-matching-crossing.md)
- [A constant-composition obstruction after an isotropic basis change](isotropic-three-label-barrier.md)
- [Interpolation, addition, and carry blocks close under Fourier queries](linear-carry-blocks.md)
- [Half-rate deterministic factoring in linear space](linear-space-half-rate.md)
- [A multiplication network with logarithmic wire widths](logarithmic-wire-network.md)
- [A multigraded crossing on uniform set partitions](matching-code-lift.md)
- [Homogeneous matching duals from Jack integrality](matching-gram-unipotence.md)
- [A minimum-cut rank budget and a twenty-one-label cutoff](minimum-cut-rank-budget.md)
- [Streamed modular factors in every fixed characteristic](modular-factor-space-odd.md)
- [Polynomial-space evaluation of the ambient modular factors](modular-factor-space.md)
- [Linear space for the regular modular factors](near-linear-factor-space.md)
- [Self-dual nilpotent cycle labels](nilpotent-cycle-labels.md)
- [Research progress: 30 October 2026](oct30-progress.md)
- [Complementary-pair labels give exact subset-intersection blocks](paired-mode-johnson.md)
- [Disjointness minors from paths and cycles](partial-matching-minors.md)
- [Coordinate restrictions of the four-permutation seed](permutation-coordinate-restrictions.md)
- [Large homogeneous permutation label spaces](permutation-homogeneous-labels.md)
- [Single-cycle blocks obstruct the full permutation coefficient span](permutation-span-cycle-obstruction.md)
- [An explicit top clique with polynomial-time table queries](polynomial-top-queries.md)
- [Prime-power Gram involutions and streamed Schur factors](prime-power-gram-involution.md)
- [Rational crossings need no separate polynomial coefficient-height bound](rational-height-transfer.md)
- [Simultaneous time and workspace from rational crossing families](rational-simultaneous-transfer.md)
- [Far-reflected labels reduce to coordinate matching](reflected-support-reduction.md)
- [The regular homogeneous central piece has the half-rate exponent](regular-fusion-rank.md)
- [Small-factor preprocessing improves the modular transfer](rough-input-transfer.md)
- [The homogeneous Schur central piece as a fusion Gram matrix](schur-fusion-gram.md)
- [Asymptotically full Schur fusion rank](schur-fusion-rank.md)
- [Exact rank-sum searches and dual label restrictions](seed-subspace-recombination.md)
- [A simultaneous time and workspace transfer for crossing families](simultaneous-time-space-transfer.md)
- [Six labels cannot pass the rank-sum amplification criterion](six-label-rank-sum-obstruction.md)
- [Sum-of-squares identities and the two crossing ranks](sos-partial-transpose-profile.md)
- [A six-label permutation crossing with rank product 150](split-permutation-seed.md)
- [Streamed finite seeds and exact polynomial certificates](streamed-seed-transfer.md)
- [Sublinear workspace for the original top-table queries](streamed-top-space.md)
- [A transfer theorem with more than polynomially many crossing pieces](subexponential-piece-transfer.md)
- [Half-rate factoring with sublinear working space](sublinear-space-half-rate.md)
- [A swap/product rank bound and the finite rank-sum cutoff](swap-product-rank-sum.md)
- [A sharp rate barrier for symmetric three-matching pieces](symmetric-three-matching-barrier.md)
- [Two pieces cannot pass the finite rank-sum criterion](two-piece-rank-sum-obstruction.md)

## Retained PR #24 effective-family and finite-rank checkpoint

Integrated from [exact PR #24 source](https://github.com/1337isnot1337/integer-factoring/tree/05d9c1b8aed4f8494ce230f0869f9048ef879474), preserving its attribution and scope. Conditional transfers and finite barriers do not establish an exponent above `1/2`; integration is not independent mathematical acceptance.

- [Binary code restrictions: an affine-code obstruction](affine-code-barrier.md)
- [The full three-pairing type route is obstructed over every field](all-field-pairing-barrier.md)
- [Current crossing frontier and resolved shortcuts](frontier-search-ledger.md)
- [Exact finite cutoff 25 for the joint antipodal-sector span](hamming-joint-optimum.md)
- [Polynomial crossing profiles with bounded bit weights](hamming-pfaffian.md)
- [Exact cutoff four in the relabeled width-two span](hamming-relabeling-optimum.md)
- [Joint label relabelings and a product-sixteen seed](hamming-relabeling.md)
- [Crossing sectors on several antipodal pairs](hamming-sectors.md)
- [Local tensor powers and weight symmetry of the crossing deformation](hamming-tensor-power.md)
- [Complete central ranks from the small kernels of the two factors](modular-central-rank.md)
- [The exact half endpoint for the normalized modular family](modular-half-endpoint.md)
- [A reflected particle-hole identity for the complete modular profile](modular-particle-hole.md)
- [Effective recombinations of the complete modular pieces](modular-recombination.md)
- [Closing the small characteristic-two pairing seeds](pairing-small-seeds.md)
- [A characteristic-two barrier for the full pairing type span](pairing-word-barrier.md)
- [From rational-function crossings to one fixed prime field](rational-family-transfer.md)
- [Polynomial-size rank blocks and an exactly analyzed reciprocal pencil](sl2-coefficient-ranks.md)
- [Crossing families from a sparse sequence of widths](sparse-width-interpolation.md)
- [A rank barrier for position-invariant word families](word-symmetry-barrier.md)

## Retained PR #26 effective-transfer and crossing checkpoint

Integrated from [exact PR #26 source](https://github.com/1337isnot1337/integer-factoring/tree/dd28d375ea65863356c67b051f5d4f217562475d), with original credits and scope retained. Conditional transfers, universal-search observations and finite certificates do not prove a new numerical factoring bound. Integration is not independent mathematical acceptance.

- [A uniform all-field barrier for binary word and balanced-layer tensors](all-field-binary-word-barrier.md)
- [A uniform all-field barrier for symmetric binary-energy shells](binary-energy-symmetry-barrier.md)
- [An all-field binary two-piece inequality and an F4 boundary pencil](binary-two-piece-rank-sum-and-f4-boundary.md)
- [An exact bivariate trace-span control at eight labels](bivariate-two-patch-span-control.md)
- [Constructive extraction from border and nonunit crossing certificates](border-and-nonunit-seed-transfer.md)
- [Catalytic crossing seeds with a common boundary vector](boundary-stable-catalytic-seeds.md)
- [Exact flow certificates for channel-walk bounds](channel-flow-dual-certificates.md)
- [Matrix-valued crossing seeds with commuting factor halves](commuting-halves-matrix-seed-transfer.md)
- [Component-path refinement of channel-walk bounds](component-channel-walk-transfer.md)
- [Diagonal-pair rank-one auxiliary transfers](diagonal-rank-one-word-transfer.md)
- [Eight local parity features in a constant-bond crossing](eight-feature-parity-crossings.md)
- [Exact exponential growth of matching channel types](exact-channel-type-growth.md)
- [The exact finite divisor-count parameter](exact-divisor-count-parameter.md)
- [Exterior monodromy and occupation-graph rank certificates](exterior-monodromy-and-occupation-graphs.md)
- [Exact channel-support feasibility by bounded binary carries](finite-carry-channel-feasibility.md)
- [Finite crossing rank profiles can be realized over finite fields](finite-field-rank-specialization.md)
- [Finite-field functional certificates for observable transfers](finite-functional-rank-detector-certificates.md)
- [Finite-state polynomial crossing families](finite-state-pairing-crossings.md)
- [Four-feature matrix-product crossing identities](four-feature-parity-crossings.md)
- [A four-piece F4 boundary pencil with two free ket functions](four-piece-f4-boundary-pencil.md)
- [A uniform score transfer from arbitrarily sparse finite seeds](gap-free-seed-score-transfer.md)
- [Generic full factor pairs need not admit a common scalarization](generic-full-factor-nonscalarization.md)
- [The logarithmic correction in the existing half-rate channels](half-rate-second-order-audit.md)
- [Rank-sensitive Fourier contraction with digit restrictions](hankel-digit-box-contraction.md)
- [What a lacunary characteristic-two kernel would still need](lacunary-kernel-boundary.md)
- [Linear rank detectors for observable transition tables](linear-rank-detectors-for-observable-tables.md)
- [Local alphabet bounds for observable word-rank transfers](local-alphabet-and-pairing-spectral-barriers.md)
- [A matrix-unit monodromy crossing and its product-boundary restriction](matrix-unit-monodromy.md)
- [A constructive commutator projector for weights 1 and 3](mismatched-commutator-projector.md)
- [Mismatched mode weights: a noncommutative Weyl family and a spectral obstruction](mismatched-mode-weyl-crossings.md)
- [A characteristic-two barrier for all binary position-invariant word tensors](modular-binary-word-barrier.md)
- [Realizing an observable matrix table as a crossing transfer](observable-transition-realization.md)
- [Odd-characteristic Weyl spectra and a finite component certificate](odd-characteristic-weyl-spectra.md)
- [An integer divisor-count bound with leading coefficient one](optimized-divisor-count-bound.md)
- [Exact certificates for overlapping three-label subalphabets](overlapping-subalphabet-spectral-certificates.md)
- [Base-prime digits, exact carries, and the cost of lifting Fourier queries](padic-digit-carry-transfer.md)
- [Partial-mode Weyl curves and an exact ordered endpoint construction](partial-mode-weyl-sweeps.md)
- [Minimum-channel jet factorizations with polylogarithmic workspace](polylog-space-jet-factorization.md)
- [Matrix rank and affine recovery with polylogarithmic workspace](polylog-space-rank-recovery.md)
- [Polynomial shear controls and a partial long-profile checkpoint](polynomial-shear-diagnostics.md)
- [A reduction-map Fierz crossing and its exact word bound](reduction-map-fierz-crossing.md)
- [Product-subspace automata and joint scalarization](scalarized-channel-automata.md)
- [Further-bound research session: 30 October 2026](session-2026-10-30-further-bounds.md)
- [Bivariate site-weighted crossing identities](site-weighted-pairing-crossings.md)
- [Quasipolynomial-time matrix queries with linear auxiliary space](small-space-berkowitz-rank.md)
- [Complete small controls for a tagged trace coordinate](tagged-two-patch-span-controls.md)
- [An exact binary control for three rank-one boundary pieces](three-piece-binary-boundary-control.md)
- [Binary trace kernels with small crossing ranks and large ordinary rank](trace-bilinear-conservation-kernels.md)
- [Exact syndrome ranks for restricted bilinear trace kernels](trace-kernel-restriction-and-composition.md)
- [A smaller binary trace-kernel counterexample with exact support](trace-kernel-small-certificate.md)
- [Scalar twists and finite channel walks for matrix-valued seeds](twisted-channel-walk-transfer.md)
- [Five-label twisted-pairing preflight](twisted-pairing-observable-preflight.md)
- [Exact recombination controls for two small trace spans](two-patch-span-recombination.md)
- [A complete two-patch trace crossing and its finite rank control](two-patch-trace-crossing.md)
- [Exact rank compression for two-state path crossings](two-state-run-rank-compression.md)
- [One factoring machine can select its finite field with its seed](universal-field-seed-selection.md)
- [Universal verified search already attains the supremum factoring score](universal-verified-factoring-search.md)
- [A score-preserving transfer for arbitrary vanishing width error](vanishing-slack-score-transfer.md)
- [Weighted binary-energy shells do not rescue the reciprocal pencil](weighted-energy-reciprocal-pencil.md)
- [The exact multiplication law of the weighted Weyl coefficients](weighted-weyl-additive-row-law.md)
- [Weighted Weyl crossings and multivariate coefficient extraction](weighted-weyl-and-multivariate-extraction.md)
- [A protected minor in the weighted Weyl symmetry span](weighted-weyl-protected-minor.md)
- [A finite square-zero ideal in the full parameterized Weyl algebra](weighted-weyl-rank-one-ideal.md)
- [An all-energy row filtration for the binary weighted Weyl curve](weighted-weyl-row-filtration.md)
- [The constant alternate-order Weyl algebra](weyl-alternate-order-algebra.md)
- [A p-typical Witt state space with perfect modular duality](witt-duality-and-vertex-jets.md)

## Retained PR #27 conservation and crossing checkpoint

Integrated from [exact PR #27 source](https://github.com/1337isnot1337/integer-factoring/tree/b49b9e19480c65d139ef15beecacb92064c1a960), retaining original credits and scope. Ordinary transfers and restricted rank barriers do not establish a stronger complete factoring exponent.

- [An exact exterior-algebra model for both crossing cuts](bi-exterior-crossing-algebra.md)
- [Degree and quadratic-rank budgets for the full exterior model](bi-exterior-degree-budget.md)
- [Verified divisor recovery in a varying characteristic](characteristic-independent-recovery.md)
- [Dual code sections: exact normal forms and complete finite screens](code-section-normal-forms.md)
- [A rank-sum barrier for all unsymmetrized diagonal words](diagonal-word-rank-sum-barrier.md)
- [A rank quotient for arbitrary fixed-target Feistel relabelings](fixed-target-feistel-barrier.md)
- [Flag quotients exclude finite rank-sum seeds through eighty-two labels](flagged-cut-rank-budget.md)
- [A common rank quotient for all affine-selector frames](frame-filtration-barrier.md)
- [A rank-sum obstruction for every binary conservation kernel](full-conservation-rank-sum-barrier.md)
- [Formal crossing limits are excluded through sixty-four labels](geometric-border-cutoff.md)
- [Geometric compression excludes exact rank-sum seeds through ninety-five labels](geometric-cut-rank-budget.md)
- [The geometric compression bounds still admit the dyadic scalar profile](geometric-cut-relaxation-limit.md)
- [Universal subspace compression from integer Jacobian certificates](geometric-subspace-compression.md)
- [One global block length permits growing coefficient fields](global-block-field-transfer.md)
- [A support barrier for Boolean-incidence crossing kernels](incidence-kernel-rank-sum-barrier.md)
- [Nonlinear label maps give a genuinely new layer family, but no gain yet](nonlinear-layer-probe.md)
- [Pauli/Walsh atoms: an exact sixty-four-label dictionary exclusion](pauli-walsh-dictionary.md)
- [Polynomial-time coefficient queries on quadratic Boolean codes](quadratic-code-oracle.md)
- [Factor kernels exclude formal crossing limits through seventy-one labels](rank-one-formal-cutoff.md)
- [A quantitative global transfer from raw rank sums](rank-sum-global-transfer.md)
- [A two-product approximation leaves a half-dimensional complement](two-product-complement.md)
- [The two-product theorem excludes rank-sum seeds through 54 labels](two-product-seed-cutoff.md)
- [Weighted cut budgets survive formal one-parameter degenerations](weighted-cut-border-stability.md)
- [Positive cut duals exclude seeds through fifty-two labels](weighted-cut-dual-certificate.md)
- [Weighted cut budgets and a forty-six-label cutoff](weighted-cut-rank-budget.md)
- [An exact limitation of the scalar cut budgets](weighted-cut-relaxation-limit.md)
- [Residual cut budgets exclude seeds through fifty labels](weighted-cut-residual-budget.md)
- [Arbitrary binary-word restrictions of the 23-term profile](word-restriction-barrier.md)

## Retained PR #29 Ward and crossing checkpoint

Integrated from [exact PR #29 source](https://github.com/1337isnot1337/integer-factoring/tree/7154cec09b982b54a0fcc0bb7a303318fc1f90da), retaining its original credits and scope. Entry algorithms, conditional interfaces and construction-specific finite or all-size barriers do not establish an exponent above `1/2`.

- [Complete control permutations and bipartite Pfaffian codes](additional-crossing-screens.md)
- [Quadratic charges, residue syndromes, and label entropy](binary-quadratic-charge-shells.md)
- [Moment identities and a cumulative bosonic crossing profile](bosonic-hermite-profile.md)
- [Crossing families from coded central quotients](coded-central-crossings.md)
- [Common-functional rank compressions without coefficient algebras](common-functional-rank-compressions.md)
- [Explicit complement shells and their finite local tensor](complement-shell-crossings.md)
- [Exact rank blocks and protected minors for complementary shells](complement-shell-rank-blocks.md)
- [Cubic crossings with arbitrary dual output maps](cubic-arbitrary-output-identity.md)
- [A bounded non-diagonal cubic probe over F8](cubic-complement-F8-probe.md)
- [Allowing a nonzero complementary output image](cubic-complementary-output-maps.md)
- [Coordinate-image output mixtures: an all-power barrier](cubic-coordinate-image-mixtures.md)
- [Cyclic complementary maps: an all-field, all-power barrier](cubic-cyclic-all-power-barrier.md)
- [Fixed-image rank-one outputs: a universal all-power barrier](cubic-fixed-image-rank-one-barrier.md)
- [Jordan complementary maps: a triangular all-power barrier](cubic-jordan-all-power-barrier.md)
- [Full mixed cubic coefficient spaces over F8](cubic-mixed-F8-spaces.md)
- [Mixing the recorded cubic output spaces](cubic-mixed-output-pool.md)
- [A gate for changing a few cubic output rows](cubic-output-perturbations.md)
- [Cubic radical codes with syndrome-preserving output maps](cubic-radical-output-maps.md)
- [Varying rank-one images: an affordable-span escape](cubic-varying-image-affordable-escape.md)
- [Structure and bounded optimization of the varying-image space](cubic-varying-image-normal-form-and-peels.md)
- [Exact weighted aggregation of complete folded-Gram channels](folded-gram-aggregation.md)
- [Complete folded-Gram channels and exact Johnson ranks](folded-gram-channels.md)
- [The first characteristic-two residual jet at even state energy](homogeneous-central-eighth-jet.md)
- [Quotient symbols for the modular partition crossing](homogeneous-central-modular-symbols.md)
- [Complete homogeneous profiles and their first nonzero jets](homogeneous-central-profiles.md)
- [Homogeneous modular duals at every label width](homogeneous-central-selection.md)
- [Mixed products of the certified coefficient spaces](mixed-common-functional-barrier.md)
- [A nonbinary output section for the 24-label code](nonbinary-two-even-site-sections.md)
- [Abelian orbit kernels: midpoint blocks and an all-cyclic exclusion](noncommutative-crossing-abelian-midpoints.md)
- [Finite rank-sum gates for crossing constructions](noncommutative-crossing-finite-gates.md)
- [Graded Hopf crossing families and scalar rank factorizations](noncommutative-crossing-graded-hopf.md)
- [Class-two group-word crossings with arbitrary central controls](noncommutative-crossing-group-degeneration.md)
- [Involution-free minors of group orbit tensors](noncommutative-crossing-involution-support.md)
- [Optimal common retracts that cancel two Pfaffian vertices](noncommutative-crossing-optimal-retracts.md)
- [Crossing pieces from noncommuting group words](noncommutative-crossing-orbits.md)
- [A PBW degeneration closes the graded Hopf crossing lane](noncommutative-crossing-pbw-degeneration.md)
- [Pfaffian boundary sections](noncommutative-crossing-pfaffian-boundaries.md)
- [Sparse isotropic Pfaffian sections](noncommutative-crossing-pfaffian-isotropic.md)
- [Ward factors from commuting symplectic forms](noncommutative-crossing-radical-ward.md)
- [Sparse symplectic radical codes with two vanishing endpoints](noncommutative-crossing-symplectic-radicals.md)
- [Charge filtration for the Ward rank factors](noncommutative-crossing-ward-charge-filtration.md)
- [Euler refinement and a small Ward interface](noncommutative-crossing-ward-euler.md)
- [Graded rank descent after the Ward factor](noncommutative-crossing-ward-rank-descent.md)
- [A uniform binary barrier for cancelling P6 boundary sections](pfaffian-p6-boundary-orbit-barrier.md)
- [Exact one-particle reduction of the four-state Pfaffian](pfaffian-particle-hole.md)
- [Exact ranks of Pfaffian retract families](pfaffian-retract-rank-profiles.md)
- [A three-parameter Pfaffian family on fixed-length codes](pfaffian-ternary-edge-family.md)
- [Exact high-mode residue folding and two-cut channel dimensions](residue-folded-channels.md)
- [Simultaneous time and workspace at the finite-seed endpoint](simultaneous-seed-endpoint.md)
- [A certified limitation of the stabilized centered-bound computation](stable-boundary-envelope.md)
- [Complete coordinate changes in the stable Pfaffian evaluator](stable-pfaffian-coordinate-search.md)
- [A stable two-factor decomposition of the complete Pfaffian](stable-pfaffian-modules.md)
- [Exact rank profiles for symplectic boundary codes](symplectic-boundary-code-profiles.md)
- [Mixed exterior grades: endpoint contacts and exact repairs](symplectic-mixed-grade-codes.md)
- [The coefficient algebra of the 26-label radical seed](symplectic-radical-seed-algebra.md)
- [A common-character barrier for every power of the 24- and 26-label seeds](symplectic-seed-character-barrier.md)
- [Exact field screens for shared-parameter seed powers](symplectic-seed-tensor-field-screens.md)
- [Exact small modules for the 26-label radical seed](symplectic-seed-tensor-modules.md)
- [Weighted radical-seed powers and exact small modules](symplectic-seed-weighted-powers.md)
- [A sumset lower bound from the two syndrome block decompositions](syndrome-sumset-rank-cost.md)
- [Bounded Frobenius modules at ten sites](two-even-site-P10-frobenius-probe.md)
- [An all-field, all-power barrier for the 240-label seed](two-even-site-P6-character-barrier.md)
- [A 672-label all-field certificate from small Jordan blocks](two-even-site-P7-jordan-barrier.md)
- [Eight-site barrier from bounded Jordan modules](two-even-site-P8-character-barrier.md)
- [All-size assembly from four-site cores](two-even-site-core-factorization.md)
- [Small row corners close the nine-site coefficient span](two-even-site-row-corner-barrier.md)
- [Two-even-site codes with different output sections](two-even-site-section-search.md)
- [Independent operator and state grades in the complete Gram construction](unequal-grade-folded-gram.md)
- [Computed Euler and charge-filtered Ward bounds](ward-charge-bound-computations.md)
- [Label entropy for binary Ward codes](ward-code-charge-entropy.md)
- [Computable Ward bounds and actual overlapping codes](ward-code-rank-bounds.md)
- [Complete coefficient ranks with the Ward operator retained](ward-combined-operator-ranks.md)
- [Sixteen scalar channels for disjoint Ward forms](ward-disjoint-operator-transfer.md)
- [The complete lowest Ward boundary vanishes for the 84-label fan](ward-fan-boundary.md)
- [The 84-label fan has exact Ward vanishing order four](ward-fan-exact-order.md)
- [An overlapping Ward family with interval width two](ward-fan-family.md)
- [An exact vanishing block improves the first fan upper profile](ward-fan-leading-block.md)
- [Two-state formal Ward queries for the fan](ward-fan-polynomial-queries.md)
- [Several Ward forms in the same eight-dimensional operator algebra](ward-form-operator-transfer.md)
- [Ward entry evaluation controlled by support overlap](ward-interval-operator-transfer.md)
- [Polynomial-time physical Ward entry queries](ward-pfaffian-entry-queries.md)
- [An asymptotic floor for the computed Ward rectangular bound](ward-rectangular-bound-envelope.md)
- [An eight-dimensional interface for the actual Ward operators](ward-stable-operator-interface.md)
