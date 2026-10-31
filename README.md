# Deterministic integer factoring bounds

[![Verification](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

This repository studies worst-case deterministic factoring. Given an integer `N >= 2` in binary, the machine outputs its smallest prime factor. Runtime is measured in finite-alphabet Turing-machine steps, with `n` equal to the input bit length.

**Want to contribute?** [Start here](docs/START_HERE.md) for setup, a reading path,
and small tasks. Reproductions, mathematical corrections, documentation, independent
review, and scoped negative results are welcome.

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

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
make verify-docs          # Local Markdown-link check
make verify-schur-fermion # Complete Schur entry formulas and exact comparisons
make verify-schur-central # Hook-residue identities and full Schur ranks through weight 8
make verify-schur-gaussian # Exact digit reduction and fixed Gaussian crossing seed
make verify-schur-sparse # Sparse down-sets and two-scale offset interval reductions
make verify-fixed-modes # Arbitrary mode sets, independent occupations, and exact endpoints
make verify-cancellation-modules # Exact module states and word-rate certificates
make verify-minimal-modules # Exact boundary reduction and complete word-rate hierarchy
make verify-exact-modules # Complete rational lattice orbits and exact pressures
make verify-cancellation-crossing # Strict certificate gain for a complete crossing
make verify-block-smith # Exact profiles for ordered crossing blocks
make verify-cubic-full-images # Full cubic coefficient spaces and reduced ranks
make verify-linear-smith # All binary scalar pencils and the uniform linear barrier
make verify-type-pressure # Exact common-type optimization and physical gauges
make verify-crossing-interpolation # Exact synthesis with prescribed local row lattices
make verify-smith-oracle # Exact weighted profiles and the scoped sparse entropy dual
make verify-native-smith # Optional compiled F4 search controls; needs a C++17 compiler
make verify-native-rectangular # Optional rectangular paired-type search controls
make verify-controlled-smith # Controlled monomial crossings and the qutrit control
make verify-sparse-holonomy # Exact singleton profiles and Gaussian entry identities
make verify-charge-pruning # One-sided charge changes and typed Gaussian block bounds
```

The Python checks use Python 3 and its standard library. Some optional checks use NumPy for speed. The exploratory mode-pressure search uses NumPy and SciPy. The optional linear Smith catalogue and Schur central-rank checks also require a C++17 compiler.

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

The latest [research checkpoint](research/constructive-frontier-checkpoint.md)
keeps the proved score at `1/2`. An
[accelerated exact certificate](research/accelerated-interval-certificates.md)
now treats `2^209-1` modes with a rational upper bound less than `5e-123`
above `1/sqrt(2)`. This remains an approach from above. A finite identity
with a positive remainder makes the check polynomial in the logarithm of
the mode count and can also certify long Gaussian spectator products.

Symbolic entry witnesses now support a
[uniform balanced-gauge obstruction](research/rectangular-balanced-gauge-obstruction.md)
for the baseline repair and its
[full-alphabet extension](research/full-alphabet-leading-obstruction.md).
The [all-scale extension](research/all-scale-balanced-obstruction.md)
now proves positive critical pressure for every `log s > 0` once the
interval is sufficiently large. Both the weights and the scale may
change arbitrarily with the interval length. The result concerns these
specific repairs and balanced bulk charges. Searches continue in families
that alter their preserved background, including
[opposite-charge cubic terms](research/rectangular-core-monomials.md).
Four [critical cubic pairs](research/full-alphabet-monomial-residuals.md)
also change the residual itself and restore all sixteen states at the
singular canonical quartic coefficient. Their determinant is exactly
`(1+gamma+ab)^4`, and a height-block correction proves the all-length
crossing. Its partial repairs still preserve an all-low Gaussian submatrix;
their occupation mixing requires a separate leading-pressure analysis.

An [arbitrary frozen-core output row](research/polynomial-frozen-core-crossings.md)
now gives an eleven-parameter affine family of locally polynomial crossings.
Its endpoint condition is exactly `g_0 + g_10 = 1`. The family contains both
the Gaussian background and a full-inverse restriction with an exact
non-Gaussian coefficient witness. This changes the background excluded by
the fixed-repair obstruction, while a useful pressure remains to be found.
The [compact frozen-row compiler](research/frozen-core-row-profiles.md)
reduces both cuts to matrices of order at most 64, independent of interval
length, and normalizes the nonzero output parameters to one orbit.
Its complete-polynomial cancellation search records 24,528 profile queries
through interval length 64. The Gaussian row was best in those finite
searches; arbitrary gauges and affine rows remain open.

A [strict Pascal repair](research/strict-gaussian-relabeling-crossings.md)
preserves the inherited Gaussian profiles and makes the local endpoint an
exact crossing on the full alphabet. It permits affine relabeling sums,
binary charge projections, and nonseparable entrywise monomial changes.
The charge projections split both cuts into explicit square blocks with
no degree increase. Full-letter gauge searches are recorded separately
from exact endpoint and profile controls; no threshold witness is claimed.

[Slow endpoint corrections](research/slow-ramification-crossings.md) after
ramification have exact eventual affine Smith orders, with an explicit
threshold. An exact qutrit example decreases the profile product by about
1.47% at one fixed scale. It does not improve the Gaussian optimum; its
finite correction and long-interval embeddings remain constructive search
targets. Both endpoint admissibility and the two complete profiles are
checked before interpreting a favorable local cost.

The [Walsh permanent calculation](research/walsh-principal-permanents.md)
gives all-size scalar valuation formulas beyond the determinant model.
It also proves that exponentially many grouped modes survive the current
precision schedule for direct forward twirling. Lower-precision scalar
recovery remains a distinct possible route.

The [arithmetic Cauchy theorem](research/arithmetic-cauchy-factorial-law.md)
now proves the permanent valuation `v_2(per C)=v_2(n!)` for every
even-spacing arithmetic block, with no power-of-two size restriction.
It also proves every invariant-factor order of the squared Cauchy matrix:
`2*r*v_p(d) + 2*v_p(r!) + v_p(r+1)` when `p` divides the spacing `d`.
An integral orthogonal-polynomial basis supplies the full Smith proof.
The [all-power extension](research/arithmetic-cauchy-powers.md) now gives
`2*r*v_p(d) + v_p(r!) + v_p((m)_r)` for every positive Hadamard power `m`,
where `(m)_r` is the rising factorial. A factorial-normalized alternant
bound makes the leading Cauchy--Binet term unique modulo `p`. These are
complete spectral formulas; scalar recovery and grouped operator cost
remain separate obligations.

The [digit-count interval theorem](research/truncated-symmetric-interval-frontier.md)
also proves that the displayed truncated-symmetric interval pressure stays
strictly above the threshold for every integer digit count `q >= 2`,
every mode count, and every profile scale. Its fixed-`q` gap has a sharp
Lambert-W asymptotic. Noninterval modes and independent mode heights remain
outside this theorem. The [weighted Pascal reduction](research/weighted-lucas-pascal-profiles.md)
now computes admissible independent-height profiles by constant-field elimination.
Lucas antichains permit every height vector, and a prime-radix traversal skips
arbitrarily large gaps between modes. These are exact search tools; the original
crossing degree remains part of the eventual evaluator cost. Run
`make verify-frontier-extensions verify-weighted-pascal` for the ordinary-proof
controls; neither result changes the factoring score.

The [mixed-field rectangular repair](research/rectangular-f4-repair.md) now
settles the four-element-field construction left open in the preceding
checkpoint. Its explicit correction preserves the low Gaussian background
and has at most a simple local pole for every finite mode set. Matching
paired types give polynomial crossings. The optional rectangular profile
engine now accepts F4 coefficients and independent physical types; its
controls replay exact rational-weight profiles. These are construction and
search advances, with no score improvement claimed.

For full intervals, the [bounded-defect reduction](research/rectangular-defect-factorization.md)
now confines this repair to seven active modes at infinity and nine at zero,
at the explicitly stated Gaussian gauges. The remaining modes contribute
scalar factors. Exact rational support checks and an independent full-matrix
control support a scalable profile calculation. The
[paired-type search audit](research/rectangular-convex-pressure.md) records
the derivative, minor-separation formulation, and precision conditions.
Regular core coordinate bases now permit fourteen independent core weights
and four bulk charges without increasing the active dimensions. The
[secondary-height family](research/rectangular-refined-repairs.md) adds
input and output repairs with a shared additive flag, together with
arbitrary endpoint-compatible Laurent quartic numerators. A
[full-core extension](research/rectangular-full-core-factorization.md)
uses ten active modes at both places and supports all four extreme
physical cores. Its 25-parameter pressure calculation has dimensions
independent of the interval length. Adaptive exact series windows permit
larger intervals while rejecting unresolved pivots.

The [full-alphabet quartic repair](research/full-alphabet-quartic-repair.md)
retains all four extreme states on both crossing pairs. Its residual is
invertible for every quartic coefficient other than one, and a sparse
partial correction gives non-Gaussian examples over F4. Fixed-size profile
matrices support 23 or 29 continuous parameters for this family.

[Coupled core transports](research/rectangular-core-transports.md) now
permit arbitrary invertible core changes, including height-raising
directions, when an exact left/right product preserves the endpoint.
Their matrices can be updated from cached coordinates. The
[critical-pressure analysis](research/rectangular-critical-pressure.md)
separates genuine defect gains from interval convergence, while a
[bivariate limit](research/bivariate-core-pressure-limit.md) reduces the
large-interval question to a finite leading matrix and its weighted ranks.
The [leading-face certificates](research/rectangular-leading-face-lp.md)
give exact rational LP bounds for the fitted rectangular model. Its
[all-low background](research/rectangular-low-background.md) is exactly a
shorter Gaussian interval. Separate [polynomial shears](research/rectangular-polynomial-shears.md)
and [four extra quartics](research/rectangular-extra-quartics.md) preserve
the endpoint and provide additional admissible search directions.

The [all-field Pascal construction](research/all-field-pascal-crossings.md)
retains the required global signs in odd characteristic. The
[graded-power theorem](research/graded-even-power-transfer.md) also
permits even perturbations that vanish at the crossing endpoint. The
scalar factoring dependency is pinned in a separate
[source audit](research/fixed-prime-transfer-source-audit.md).
A [truncated symmetric construction](research/truncated-symmetric-pascal-crossings.md)
also gives ordinary tensor powers with p states per mode in characteristic
p, including non-Gaussian perturbations of odd degree.

For the original evaluator's precision question,
[canonical packet bounds](research/canonical-packet-upper-frontier.md)
and [finite decoder certificates](research/canonical-packet-decoder-frontier.md)
give exact four-wire binary thresholds: 12 local bits for parity and
13 for the full count. They do not supply an efficient all-size decoder.
The [full Boolean operator](research/boolean-cauchy-smith-profile.md)
has an exact all-size Smith profile and a one-bit inverse at either
nearest center. Direct insertion of that inverse reconnects the source
groups, so this operator theorem alone does not improve the evaluator.
The [grouped inverse profile](research/grouped-cauchy-inverse-profile.md)
determines the exact precision cost of a direct grouped factorization.

The [31 October research checkpoint](research/constructive-frontier-checkpoint.md)
collects the latest exact bounds, new construction families, reproducible
search starts, and the remaining obligations for a score improvement.

No uniform factoring exponent above `1/2` is proved here. One open approach seeks effective crossing pieces with total rank width `4*rho*lambda + O(log lambda)` for `rho < 1/sqrt(2)`. A suitable family could improve the evaluator bound, but finite rank improvements alone do not establish that inequality.

The [complete Schur entry algorithm](research/schur-fermion-queries.md) computes
individual modes with determinants of order at most `3 sqrt(H)`, and full
crossing coefficients with determinants of order at most `4 sqrt(h)`. It works
directly in finite characteristic and removes factorial-precision conversion
from these queries. Its ordinary proof and independent exact checks improve
the available construction tools. For the homogeneous Schur family, [hook-residue identities](research/schur-hook-residue-central.md)
give central rank at most `p(2h)-3` at every power-of-two weight. Full exact
calculations attain this bound at weights 8 and 16, with ranks `228` and
`8,346`. These central pieces exclude the corresponding raw finite seeds;
no stronger factoring bound is claimed.

The [complete inhomogeneous Schur kernel](research/schur-gaussian-digit-reduction.md)
has an exact binary-digit reduction to a four-label Gaussian crossing.
The proof retains its exterior maps and energy projection and gives an
explicit finite coefficient decomposition. Fixed-occupation sectors of
the seed are complete crossings in their own right.
A [sparse-index calculation](research/sparse-holonomy-singleton-profile.md)
gives every Smith order for the zero-and-singleton index set and an exact
Gaussian encoding for all Boolean down-sets. It supplies formulas for
larger constructions; it does not assert a better factoring exponent.
[Boolean down-sets](research/schur-downset-crossing.md) give complete
Schur crossings with matrix-valued holonomy and controlled profile limits.
[Offset intervals](research/schur-offset-interval-reduction.md) reduce to
an explicit two-scale four-label kernel, with polynomial projection and
coefficient costs when the ambient offset grows polynomially.
[Two unequal intervals](research/schur-two-block-holonomy.md) split exactly
into scalar sectors; the resulting direct kernels also allow interval
parameters that are not powers of two.
[Fixed mode counts](research/fixed-mode-crossings.md) extend the direct
crossing to arbitrary index sets, unequal species occupations, and
independent mode heights. Their product-type criterion is exact; an
[exploratory search driver](tools/search_schur_mode_pressure.py) uses the
weighted Smith oracle to search this larger family.
An [exact three-mode certificate](research/three-mode-pressure-improvement.md)
uses independent occupations of the two species to prove `rho < 741/1000`.
It improves the previous Gaussian pressure control but remains above
`1/sqrt(2)`, so the best factoring exponent is unchanged.
The [one-species interval family](research/schur-one-species-interval.md)
has an exact all-size spectrum and finite product pressure. Checked
rational witnesses reach `rho < 3363/4756` at 16,383 modes. This particular
formula approaches `1/sqrt(2)` from above; an analytic proof excludes
crossing that threshold within the displayed gauges.
The [quantitative interval theorem](research/quantitative-interval-frontier.md)
gives the sharp gap asymptotic `(log n)^2 / (sqrt(2) pi^2 n^2)` and
exact large witnesses using outward integer intervals.
The [accelerated extension](research/accelerated-interval-certificates.md)
replaces mode enumeration with exact geometric moments and reaches the
much larger finite certificate described above.
The [greedy Pascal-complement theorem](research/schur-pascal-complement.md)
computes the same normalized spectrum for every finite mode set using
binary linear algebra. It proves the interval formula at every length
and skips large gaps by a traversal of binary intervals.
[Constant matrix flags](research/constant-flag-crossings.md) extend the
polynomial construction beyond Pascal matrices while preserving its
determinant and typed endpoint. An exact four-mode screen over the
four-element field covers 4,096 flags; all remain at or above the Pascal
pressure at the displayed gauge.
[Positive-energy perturbations](research/positive-energy-quartic-endpoints.md)
also preserve the typed crossing when they do not vanish locally at the
endpoint. A strict index-sum condition even permits rational local kernels
whose endpoint denominators cancel in every fixed-type tensor power.
For [localized quartic changes](research/localized-quartic-profiles.md),
exact determinant and complete zero-profile formulas quantify the cost
of the resulting spread in Smith orders. A rank-deficient limiting
endpoint is explicitly checked and requires a further restriction.
[A polynomial physical code](research/residual-quartic-codes.md) repairs
that endpoint and cancels its denominators term by term. It gives an
explicit eight-letter, degree-32 candidate whose quartic correction
survives common types. Its endpoint is proved; a pressure improvement
has not been established.
[Arbitrary-mode codes](research/coded-quartic-mode-families.md) extend this
construction to noninterval sets and independent Laurent coefficients on
the four physical maps, with a reusable exact constructor.
[Conjugate phase codes](research/conjugate-residual-phase-codes.md) pair
two four-element-field tensors to cancel their endpoint phases. Three
checked eight-letter, degree-48 crossings descend to binary coefficients
through an explicit constant change of basis.
[Fixed field descent](research/fixed-field-descent.md) also transfers any
coefficient certificate over a fixed extension to its base field without
enlarging the physical alphabet. Projecting after each complete tensor
power costs only the extension degree in the coefficient-cost sum.
An [all-field residual obstruction](research/residual-quartic-three-letter-obstruction.md)
excludes three-letter restrictions by four constant physical maps over
every characteristic-two field. When both crossing pairs use the same
code, the resulting three-letter operator is always singular. Larger
blocks and bulk-dependent codes remain outside this obstruction.
A [signed exterior model](research/signed-residual-codes.md) extends the
single-copy obstruction to every field. Exact identities over the
27-element field also give a three-conjugate residual phase repair;
its full polynomial factoring construction remains an open obligation.
[An asymmetric three-letter core](research/involutive-core-crossings.md)
gives polynomial candidates on `3 * 2^(n-2)` letters for arbitrary finite
mode sets. A cancellation identity removes the local denominators, and
the remaining endpoint error vanishes on every common word type.
Exact six-letter and twelve-letter examples, monomial code parameters,
and sparse output corrections are supplied for pressure searches.
[Coordinate core restrictions](research/coordinate-core-crossings.md)
retain a filled low-mode state and give a second arbitrary-mode family
on the same alphabet sizes. Its local entries can have genuine poles;
matching word types on each crossing pair cancel them globally and
justify the restored rational Smith profiles. The two pairs may have
different types, and pairing with an exchanged copy restores equal
external dimensions.
[Independent paired-type pressure](research/paired-type-crossing-pressure.md)
gives the exact entropy optimization for different types and rectangular
alphabets. Tensoring with the pair-exchanged copy produces an equal-alphabet
crossing; keeping the two coefficient indices separate makes the costs
multiply exactly.
[Rectangular residual repair](research/rectangular-residual-repair.md)
retains four extreme states on one crossing pair and three on the other.
An explicit two-parameter inverse and a smaller height-triangular repair
give valid paired-type families. Finite mixed examples are locally
polynomial, and the smallest partial repair has an exact non-Gaussian
Wick witness. Their pressure remains a search problem.
[Even-core postselection](research/even-core-postselection.md) gives a
four-letter polynomial family `A + lambda B + lambda^2 C` with a strict
crossing endpoint. Two exact Wick identities prove that every nonzero
rational parameter path is non-Gaussian. Monomial parameter paths are
available for pressure searches; no improved rate has been certified.
The [two-output extension](research/two-parameter-even-core.md) classifies
its Gaussian members exactly: the coordinate axes and one exceptional
parameter point. An explicit eight-edge factorization identifies that
point, so numerical improvements near it can be interpreted correctly.
The [four-boundary extension](research/four-parameter-even-core.md)
classifies the larger Gaussian set as well. Every Gaussian member with
all four parameters nonzero is a scalar multiple of that same exceptional
tensor; non-Gaussian parameter searches can exclude this repeated family.
The [optional native F4 search engine](research/native-weighted-smith-f4.md)
provides exact field arithmetic, support-block decomposition, weighted
profiles, and pressure gradients. Its finite controls use an independent
rational oracle. Floating search results still require exact certificates.
An [optional rectangular engine](research/native-rectangular-smith-search.md)
handles unequal alphabets and independent paired types with the same
explicit denominator bookkeeping and exact rational controls.
For this Gaussian seed, [generic mode blocks](research/typed-gaussian-generic-blocks.md)
force cost `Omega(D^2)` when all nonvacuum occupations stay fixed, even
after arbitrary coefficient recombinations. Growing occupation densities
remain a separate regime.
[Charge-block compression](research/typed-gaussian-charge-blocks.md)
computes the actual coefficient ranks and recombinations at larger types.
For 1,260 labels it replaces the full 1,587,600-dimensional flattening by
38 orbit blocks, each of size at most 1,260. The checked decompositions
remain above the improvement threshold.

[Exact type optimization](research/type-optimized-crossing-pressure.md)
reduces the best common letter frequencies to one rational pressure inequality.
It applies to Smith and module certificates. A checked physical gauge gives
`rho < 61/80` for the four-label Gaussian seed; the three-mode certificate
above improves this local pressure bound. Both remain weaker than the
known half-rate construction.
[Polynomial interpolation](research/polynomial-crossing-interpolation.md)
constructs crossings from prescribed local row lattices by exact linear
algebra. The complete fixed-lattice fiber consists of 64 local changes, and a
relaxed system constructs a non-Gaussian control.
An [explicit tensorized family](research/tensorized-smith-sections.md)
adds 40 new directions to a 180-dimensional inherited space. A common
nilpotent residue certificate proves that all `2^220` crossings preserve
the two specified profiles; changes of gauges remain a separate question.
An [exact weighted Smith oracle](research/weighted-smith-separation-oracle.md)
finds all determinantal minima without listing every minor. It supports
continuous gauge searches with exact rational replay.
[One-sided charge changes](research/one-sided-charge-pruning.md) disappear
exactly under common-type restriction and cannot improve either convex
Smith-profile sum. Neutral interactions and opposite-charge combinations
remain available to the constructive search.

[Finite module automata](research/cancellation-module-smith-automata.md)
retain cancellation between auxiliary paths. Their bounds improve monotonically
with precision, dominate the corresponding capped entry-order bounds, and
converge to the exact full-matrix word rate when both boundaries are free.
[Reachable and observable reduction](research/minimal-boundary-module-completeness.md)
extends this completeness result to fixed scalar boundaries, while keeping
physical choices grouped by channel. This is a statement about word bounds,
not a characterization of actual coefficient ranks.
[Exact rational lattice orbits](research/exact-rational-module-orbits.md)
can retain identically missing directions and compute the full word rate
when the orbit closes. A companion example has explicit two- and seven-state
orbits and a proved certificate floor.
An [explicit coupled crossing](research/cancellation-module-crossing-example.md)
strictly improves its corresponding capped min-plus certificate, from at least
`425/27` to `1625/108`. The same local factors and precision caps are used.
The checked controls do not yet meet the strict improvement threshold.

[Fixed block Smith profiles](research/block-smith-crossings.md) form a
submultiplicative hierarchy. A three-position ordered-pairing block improves
its normalized certificate from `16` to about `15.7138`; this remains weaker
than the elementary scalar control and does not change the factoring bound.

The [degree-one Smith obstruction](research/linear-smith-barrier.md) proves
`B >= q^3` for every alphabet, field, and auxiliary dimension in the bare
linear certificate. The sharp binary bound is `256/27`; the scalar bound
in odd characteristic is `12`. All 65,536 binary scalar pencils are catalogued.
Higher degree, blocks, and cancellation refinements remain separate routes.

For the nilpotent-path construction, a [controlled monomial obstruction](research/controlled-permutation-smith-floor.md)
extends the existing permutation calculation to dependent and noncommuting
controls. It gives `B >= q^3/sqrt(2)` in its stated scope. A separate
[qutrit Weyl control](research/qutrit-weyl-crossing.md) has an exact all-length
constant-coefficient rank formula; its raw pieces also fail to improve the bound.

The [varying-image cubic construction](research/cubic-full-image-span.md)
has an exact 786-dimensional rank-one-map coefficient space, with column ranks
reduced from 49-by-49 matrices to matrices of size at most 28. Rank-two maps
already span its full 1,225-dimensional major-symmetric space. These results
define larger exact search spaces; an improving crossing decomposition is
still needed.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with the [contributor quickstart](docs/START_HERE.md) and
[CONTRIBUTING.md](CONTRIBUTING.md). You can report a correction or reproduction
without proposing a new bound. [Open a focused issue](https://github.com/1337isnot1337/integer-factoring/issues/new/choose)
or choose an [available starter task](https://github.com/1337isnot1337/integer-factoring/issues?q=is%3Aissue%20is%3Aopen%20label%3A%22good%20first%20issue%22).

For substantial research, read [autoresearch.md](autoresearch.md) and coordinate
through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8).
Reviewers can use the [review guide](docs/REVIEWING.md) and
[exact-head review ledger](docs/REVIEW_LEDGER.md). Document source revisions,
proof status, reproduction commands, and remaining gaps. The
[contribution record](CONTRIBUTORS.md) retains credit for useful partial and
superseded work as well as incorporated results.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
