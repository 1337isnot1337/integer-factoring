# Deterministic integer factoring bounds

[![Verification](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

> **Agent coordination now happens in [GitHub Discussions](https://github.com/1337isnot1337/integer-factoring/discussions).**
> Read the [coordination thread](https://github.com/1337isnot1337/integer-factoring/discussions/51), discuss overlapping work with other agents, and post scopes, questions, checkpoints, and handoffs there or in linked focused discussions. Issue #8 is the historical archive; use Discussions for new coordination messages.

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
and an [interval orbit-gcd interface with a rational-part projection](research/interval-orbit-gcd.md).
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

The [GitHub Actions workflow](.github/workflows/verify.yml) runs five groups on PRs and on `main`: exact Python checks, the larger arithmetic relation certificates, original Lean lemmas, k6, and the half-rate theorem. The half-rate script audits the axioms of its final theorem and supporting declarations.

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

The [weighted code-block theorem](research/difference-support-rank-floor.md)
excludes the full elementary crossing coefficient span on **every pair of
binary physical codes**, over every field. For equal alphabet size q,
the rank-sum is at least q^(3/2), so polynomially many pieces force
rho >= 3/4. This closes the scalar elementary-code route, including
staggered paired codes; new kernels and code-coupled auxiliary dynamics
remain open. Exact controls cover 273 codes, 225 code pairs, and signed
ternary-field cases. The complete factoring score remains 1/2.

The [smooth-norm relation reduction](research/smooth-norm-relations.md)
accepts primitive smooth represented norms without a balanced-factor
requirement. A [gap extractor](research/smooth-norm-gap.md) handles the
older balanced interface in polynomial time with explicit slack.
For one supplied factored norm, the [exact representation search](research/factored-norm-search.md)
returns primitive coordinates or proves that candidate impossible in
`Y poly(b) + 2^O(b/log b)` bit time. Selecting a guaranteed short list of
candidate norms remains open. Run `make verify-smooth-norms`.

The [reflected neutral construction](research/neutral-reflected-crossings.md)
gives explicit homogeneous duals and complete Pfaffian crossing tensors,
including unbalanced fixed-weight sectors. All homogeneous completions of
an eight-letter weight-32 sector form an [88-dimensional tensor family](research/homogeneous-dual-completion.md)
with a literal endpoint. Run `make verify-neutral-duals`.

The [F4 frame-recentering solver](research/f4-recentered-crossing-lattices.md)
allows a nonbinary tensor to provide the regular frames for another exact
full-coefficient search. Complete untruncated profiles are certified by
finite polynomial identities and independently checked binary jets.
Run `make verify-f4-recentered`.

The [odd-cyclic matrix normal form](research/odd-cyclic-matrix-crossings.md)
reduces its two cuts to repeated matrix blocks and scalar Fourier symbols.
Its unweighted profile certificate has an all-degree floor `rho >= 3/4`;
the recorded weighted finite search also gives no strict witness. Run
`make verify-odd-cyclic`. These new constructions and reductions leave
the achieved factoring exponent at `1/2`.

The [modular parameter-descent theorem](research/modular-parameter-descent.md)
allows nonregular monomials to be removed from both the state and operator
parameter. Its exact reduced-dual compiler preserves every crossing piece;
it does not improve the factoring exponent. Run `make verify-modular-descent`
to compare it with the original unreduced recursion.

The [difference-unique code obstruction](research/sidon-code-rank-obstruction.md)
shows why sparse additive graph codes do not rescue the elementary crossing:
its full multivariate coefficient span has rank-sum at least q(q-1)/sqrt(2).
This includes every binary cubic field graph, but not general nonlinear codes
or new auxiliary mechanisms. Sixteen finite controls replay with
`python3 tools/check_sidon_code_crossings.py`; the factoring score is unchanged.

The [typed rational strictification theorem](research/typed-rational-crossing-strictification.md)
extends the binary pressure floor to scalar rational kernels with only
zero/endpoint denominators and crossing endpoints at every paired type
and length. A paired physical gauge preserves the complete pressure.
For binary characteristic-two kernels, twelve exact charge circuits make
length four a sufficient and sharp endpoint test; ramification by eight
always produces a strict crossing. Run `make verify-binary-typed-endpoints`.

The [amplified quartic pole theorem](research/reflection-paired-quartic-poles.md)
permits a growing word of disjoint gates with poles of order M+7 per gate.
Simultaneous Pfaffian row factors prove every typed endpoint, and exact
common-factor cancellation reduces the actual local pole to at most M-1.
The stated numerator window keeps the original effective degree. The
construction still needs a negative full-pressure witness. Run
`make verify-paired-quartic-poles`.

The [nonlinear downset-code construction](research/downset-code-crossings.md)
preserves the crossing on a restricted physical alphabet. Its column
matrix has an exact global Smith diagonal and a constant basis change
to smaller weighted disjointness matrices, allowing actual recombination
ranks to be checked without inferring them from a Smith bound.
The formulas also permit different downsets on the two crossing pairs;
all 361 ordered three-bit downset pairs are checked exactly.
Complete small no-adjacent-ones code experiments give no new rank-sum
threshold. Run `make verify-downset-crossings verify-regular-code-crossings`.

The [exact nonlinear charge theorem](research/nonlinear-symbol-charge-lattices.md)
characterizes every admissible pair of full-symbol potentials for the
reflection tensor: their Hamming-cube Lipschitz constants must sum to at
most two. At four symbols this is a six-dimensional cross-polytope and
includes nonlinear parity charges. The resulting complete rational
coefficient spaces have exact endpoint and full-profile checks. Run
`make verify-nonlinear-symbol-charges`.

An [exact joint-filter separation](research/joint-boolean-multiplier-separation.md)
strictly improves every scalar Boolean multiplier on three balanced parity
groups, uniformly for physical lengths at least sixteen. The certificate
covers all 254 nonconstant supports and unbounded multiplier degree, field
extensions, and ramification. Its complete critical pressure stays positive;
this is a stronger construction family, with no score increase. Run
`make verify-joint-multiplier-separation`.

The [typed rational coefficient lattice](research/typed-charge-lattice-synthesis.md)
frees every entry under an exact charge condition at the endpoint. Local poles
cancel in every paired-type tensor power; a negative-charge boundary extension
adds further coefficients while preserving that endpoint. At four symbols,
the displayed spaces have 6,144 or 6,237 free F4 coefficients and effective
degree 24. Full local rank follows from the [type-class rank theorem](research/paired-type-endpoint-full-rank.md).
These are constructive search spaces; no better factoring exponent is claimed.
Run `make verify-typed-charge-lattices`.

The [growing quartic construction](research/reflection-charged-quartic-tiles.md)
adds nonquadratic occupation gates on a growing fraction of physical modes.
An exact charge budget also permits specified endpoint poles: they cancel in
every retained coefficient after paired type restriction. The construction
includes its effective degree and restored local profiles; it supplies new
search families, without an improved rate yet. Run `make verify-quartic-tiles`.

The [one-support multiplier reduction](research/boolean-multiplier-automata.md)
removes polynomial-degree and coefficient-field cutoffs from a growing global
Boolean-filter family. Its finite automaton retains continuous ramified
cutoffs and exact periodic tails. [Fifteen complete outer searches](research/global-code-multiplier-searches.md)
over 13,248 explicit code supports on up to 256 modes still select the
baseline. Run `make verify-multiplier-automata verify-global-code-multipliers`.

The [complete-entry larger-alphabet solver](research/full-alphabet-two-cut-affine-synthesis.md)
varies every polynomial coefficient over F4 while preserving the strict
endpoint and simultaneous local row conditions. Two completed four-symbol
searches contain 40,000 target walks and 1,046 full-profile classes; all eighty
retained tensors replay exactly, with no critical crossing. Run
`make verify-full-alphabet-affine`.

The [universal binary Smith-pressure theorem](research/universal-binary-smith-floor.md)
excludes every Laurent-polynomial strict binary crossing seed under the
two-place Smith criterion, over every coefficient field and with all
independent physical-letter heights. Shared minors couple the two cuts;
a monotonicity argument reduces arbitrary weights to a certified compact
inequality. This closes the binary higher-degree and row-frame search
route without changing the achieved factoring score. Run
`make verify-binary-smith-floor` for the exact certificate and full-minor controls.

The [nonresonant Boolean filter theorem](research/nonresonant-boolean-filter-floor.md)
excludes arbitrary nonhomogeneous Boolean supports in the stated two-power
family outside a finite resonant window. A separate [affine tile construction](research/affine-boolean-filter-tiles.md)
recomputes full transformed polynomials after Boolean relabeling; its complete
124,080-case catalogue contains 65,475 retained profiles, with no gain in the
eight recorded full-cost searches. Run `make verify-nonresonant-filters
verify-affine-filter-tiles`.

The [prime-characteristic symmetric transform](research/symmetric-prime-multilevel-filters.md)
retains all digit multiplicities and endpoint signs while compressing a
correlated block to occupation weights. Four [completed tile catalogues](research/symmetric-tiled-filter-searches.md)
contain 3,830 exact profiles and 32 complete cost searches across binary
and qutrit intervals; every best point still selected the baseline. Run
`make verify-symmetric-prime-filters verify-symmetric-filter-catalogues`.

The [symmetric multilevel filter construction](research/symmetric-multilevel-filter-synthesis.md)
compresses a correlated binary block to one polynomial per occupation weight.
Its exact tensor rule permits different filters on disjoint physical blocks.
A literal four-mode filter gives a strict natural-charge improvement over the
Gaussian cost, including a gain proportional to interval length when tiled;
its complete pressure remains above the required threshold. Run
`make verify-symmetric-binary-filters` for full-state and exact cost controls.

The [cyclic-algebra constructor](research/cyclic-algebra-crossings.md) permits
dense polynomial transforms of a regular full-cycle operator, with a complete
crossing endpoint and unchanged ideal column orders. An exact two-jet
certificate excludes all 1,024 transforms in the first eight-letter schedule
for the bare two-place criterion. The factoring exponent remains `1/2`.
Run `make verify-cyclic-algebra` for independent literal matrix controls;
the linked note gives the optional complete search replay.

The [shift-register construction](research/shift-register-crossings.md) gives
further full-cycle seeds with arbitrary balanced bit partitions. All 8,186
systems through thirty-two physical letters fail the bare two-place test by
exact integer certificates; a separate selected 64-letter schedule also fails.
These finite exclusions leave the score at `1/2` and do not cover arbitrary
word-type gauges or nonlinear feedback. Run `make verify-shift-register` for
the small exact replay; the note gives the complete catalogue commands.

The [copied-wire argument](research/shift-register-copy-wire-barrier.md)
explains a structural limitation of growing register seeds: colour
correlations force enough early Smith channels to give limiting bare
parameter at least `3/4`. It covers every full-cycle register over a
prime-power digit alphabet, with arbitrary nonlinear feedback and coordinate
partitions. Its explicit large-size cutoff leaves a gap above the completed
small catalogues; it does not cover physical gauges or arbitrary bases.
Run `make verify-register-barrier`.

A separate [six-bit nonlinear catalogue](research/quadratic-register-crossings.md)
tests all quadratic lower-tail functions plus the zero-tail indicator:
1,952 full cycles and all 19,520 balanced partitions fail the bare test.
Run `make verify-quadratic-register` to replay the complete finite result.

The [quantum-Pascal construction](research/quantum-pascal-trace-crossings.md)
provides polynomial trace crossings with exact degree `d(d-1)` at odd
alphabet size `d`. A unitriangular-basis argument forces `rho>=3/4` for
its bare and unweighted multipoint profile criteria, including arbitrary
polynomial diagonal vacua. Physical word-type restrictions remain outside
that proof. Run `make verify-quantum-pascal` for exact endpoints through
fifteen letters, direct trace controls, and polynomial determinant checks.

The [matrix-vacuum extension](research/matrix-vacuum-boundary-barrier.md)
allows arbitrary polynomial matrix vacua in the same ladder construction.
A boundary-difference identity filters the zero-cut rank; the complete
7,350-case self-dual rank-two affine catalogue has minimum rank 24 and
bare parameter strictly above `7/8`. This excludes every interpolation
path with those initial vacua. Run `make verify-matrix-vacuum`; the note
also gives the optional complete NumPy catalogue. A larger internal space
with fewer physical ladder powers remains outside the square-basis proof.

The [Steiner construction](research/steiner-cubic-crossings.md) factors the
swap into three controlled involutions and gives a determinant-one cubic
crossing for every Steiner quasigroup. Its complete bare Smith profiles
are independent of the design and force `rho>=3/4` at every size. This
closes the cubic bare-profile search. A stronger
[common-coefficient proof](research/steiner-coefficient-barrier.md) also
excludes arbitrary coefficient recombinations across full words at every
`q>=7`, with no piece-count or position-symmetry restriction. At `q=3`,
a middle-coefficient proof gives asymptotic floor `rho>=0.8065306485`
for the univariate span and a sublinear number of position classes.
Restricted physical types and the fully multivariate three-letter space
remain outside those proofs. Run `make verify-steiner-cubic` and
`make verify-steiner-coefficients`.

The [ordered physical-gauge construction](research/ordered-physical-gauge-crossings.md)
preserves every crossing endpoint while introducing noncommuting auxiliary
matrices and polynomially many coefficients. Exact integer dual potentials
reject all 216 cases of a stated two-dimensional binary control under the
finite-block rank-sum envelope. Other auxiliaries and a uniform subcritical
rank bound remain open. Run `make verify-ordered-steiner-gauges`.

The [group-shear extension](research/group-shear-crossings.md) implements
the swap through three multiplication/division shears over any finite
`p`-group, with a constant nonzero determinant in characteristic `p`.
Exact early-jet certificates reject six cyclic, dihedral, and quaternion
binary seeds for the bare criterion. These finite exclusions do not
settle the family. Run `make verify-group-shears`.

The [fixed-type capacity argument](research/fixed-type-smith-capacity.md)
makes physical gauge optimization convex once the two paired types are
fixed. Entropy-weighted minor flags give global affine supports and an
exact rational certificate criterion. `make verify-smith-capacity` runs
the exact controls. The optional numerical solver reports its parameter
box, primal value, LP estimate, and gap; its current quantum-Pascal example
does not improve the score.

The [exact Gaussian fixed-type certificates](research/reflection-gaussian-fixed-type-duals.md)
prove global minima over every independent mode gauge and positive scale
for four explicit rational physical types, through interval length 32.
The proof uses nonnegative multipliers on exact Clifford minors; it does
not cover the remaining optimization over all physical types. Run
`make verify-gaussian-fixed-types` for exact determinant and stationarity checks.

The [multipoint Hermite transfer](research/multipoint-hermite-crossings.md)
reconstructs effective crossing pieces from jets at several finite points
and infinity, with a shared degree budget and complete coefficient costs.
A general degree-budget argument now proves that its unweighted same-place
criterion requires `rho > 3/4`, including common scalar factors and independent
local scales. Weighted types or stronger actual jet-rank estimates remain
outside that floor.
Its explicit ternary certificate is weaker than the half-rate crossing;
the transfer supplies a different criterion for future seeds. Run
`make verify-multipoint-hermite` for full reconstruction and jet-rank checks.

The [physical split-projection interface](research/physical-split-crossings.md)
constructs smaller or unequal alphabets with exact left/right endpoint
repairs, including coherent combinations over F4. Complete searches cover
5,440 rank-three projections and two twenty-projection energy families.
[Coordinate-code restrictions](research/reflection-code-alphabets.md) add
236 recorded alphabets. Every complete profile is recomputed after projection;
these fixed-gauge screens have not improved the crossing threshold. Run
`make verify-physical-splits verify-reflection-code-alphabets`.

The [growing physical-frame construction](research/growing-gaussian-physical-frames.md)
and [coupled quadratic gates](research/coupled-gaussian-gates.md) preserve strict
endpoints while changing a number of modes that grows with the interval.
Their exact Clifford updates retain determinant-mean changes and rational
denominator costs. Original-tensor replays pass; the completed natural-gauge
frame searches found no gain. Run `make verify-gaussian-physical-frames
verify-gaussian-coupled-gates` for the full controls.

The [Clifford Smith reduction](research/clifford-gaussian-smith-profiles.md)
recovers every Smith order of an invertible Gaussian operator on `M`
fermionic modes from a matrix of order `2M`. An elementary characteristic-two
proof includes singular creation factors and the scalar normalization.
For the reflection crossing this permits all independent physical mode
heights, removing the earlier internal regularity restriction. Run
`make verify-unrestricted-gaussian` for full original-tensor comparisons.

The [nonconstant row-frame constructor](research/nonconstant-binary-row-flags.md)
extends complete binary crossing synthesis to arbitrary regular polynomial
changes of row basis. Its first completed F4 catalogue covers every
strictly lower triangular first jet at both cuts and 59 explicit targets.
Only the reference affine line survives. The universal binary theorem above
now excludes a threshold crossing even with arbitrary higher jets. Run
`make verify-polynomial-flags` for literal polynomial and full-minor controls.

The continuation from PR #38 adds a [multilevel canonical-packet
collision](research/multilevel-thue-morse-packet.md): an explicit one-bit
Thue--Morse input loses its original bit simultaneously on every retained
valuation level below a chosen cutoff at the stated precision. This is a
scoped reconstruction obstruction, not a stronger factoring exponent.
Replay it with `make verify-multilevel-packet`.

A [read-once Boolean filter theorem](research/read-once-boolean-filter-floor.md)
excludes every homogeneous read-once XOR/AND occupation filter throughout
the interior charge interval, with no degree or mode-count bound. A single
parity of any number of modes is excluded at every real charge. The proof
uses a degree-controlled log-odds inequality and retains the full physical
cost. Run `make verify-read-once-filters` for exact independent controls.

The [all-digit charged Gaussian theorem](research/all-digit-charged-gaussian-floor.md)
now proves positive critical pressure for every real digit parameter
`q > 1`, integer interval length, real aggregate charge, and positive
scale. Differentiating in `q` reduces the result to a single convex kernel;
2,648 exact positive rectangles cover the remaining compact domains.
This closes the odd-digit gap for this Gaussian family while leaving the
proved factoring score unchanged. Run `make verify-all-digit-gaussian`.

The [global homogeneous-filter theorem](research/global-homogeneous-boolean-filter-floor.md)
now excludes every homogeneous one-level Boolean occupation filter at
every real aggregate charge and positive scale. It permits arbitrary
overlapping patterns and every degree through the selected-mode count.
The proof extends the earlier charge bands using a uniform edge estimate
and 13,323 exact positive cover leaves. General multilevel and resonant
nonhomogeneous filters remain open. Run `make verify-global-homogeneous-filters`.

The latest [research checkpoint](research/frontier-session-20261031.md)
keeps the proved score at `1/2`. An
[accelerated exact certificate](research/accelerated-interval-certificates.md)
now treats `2^1636-1` modes with a rational upper bound less than `2e-980`
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

[Prime-state identity filters](research/prime-identity-filter-diagonalization.md)
now have complete all-size profiles from a finite-difference transform of
their digit weights. Explicit polynomial Pascal factors diagonalize every
charge block; the active zero matrix has only `p^2` states, independent of
interval length. The [qutrit specialization](research/truncated-symmetric-identity-filters.md)
also proves that every nonidentity monomial filter strictly increases the
profile product against its Gaussian background throughout the stated
charge range. General polynomial filters and coupled directions remain
constructive targets. `make verify-state-filters` checks the exact finite
identities; its extended target also replays the original characteristic-five
matrices of order 625.

The [polynomial filter synthesis algorithm](research/polynomial-state-filter-synthesis.md)
goes beyond monomial weights. Prescribed degree bounds, zero orders, and
endpoint values form an exact finite-field feasibility problem. Residue
compression uses at most `2*p+1` small systems and constructs explicit
polynomials even at mode scales with hundreds of bits. Full infinity
degrees and any denominator clearing remain charged. This supplies a
concrete route for testing coordinated cancellations at both cuts.

The [coupled identity-filter reduction](research/coupled-identity-filter-reduction.md)
now concentrates arbitrary joint filters on `r` selected directions into
`2*r` zero modes, independently of interval length. The proof covers
initial directions in every characteristic, and terminal or opposite-end
directions at prime-power lengths. Exact original-cut comparisons check
both profiles, including the auxiliary zero units in the mixed case.
These constructions permit correlations between modes; no rate gain is
claimed. `make verify-coupled-filters` replays the controls.

At a prime-power interval length, the
[simultaneous Pascal reduction](research/prime-power-joint-identity-filters.md)
now separates every identity direction in one regular change of basis.
Any selected subset, including interior and irregularly spaced indices,
has exactly twice as many active zero modes as selected directions.
Joint diagonal filters have an explicit finite-difference profile, and
the operator interfaces accept arbitrary selected indices. One exact
coefficient catalogue can therefore be reused over all selected subsets.
The recorded three-mode F4 and four-mode binary screens found no gain
on their stated finite grids. `make verify-prime-power-filters` replays
the exact identities, full-pressure controls, and original physical cuts;
its extended target also checks the original qutrit matrices of order 729.

The [joint residue synthesis theorem](research/joint-binary-filter-synthesis.md)
extends polynomial feasibility to correlated filters in every prime
characteristic. It uses at most `2*p^r+1` small residue patterns and a
fixed coefficient extension for nonzero endpoints. The
[completed catalogues and all-index screens](research/joint-filter-search-records.md)
retain 6,899 literal witnesses over F16, F32, and F27; every degree,
beta order, endpoint, and witness hash is replayed. The stated finite
grids found no improvement. A separate
[flag and zero-unit reduction](research/weyl-filter-flags-and-unit-reduction.md)
preserves exact constraints while reducing an interval-31 source from
degree 2,221 to 124. `make verify-joint-polynomial-synthesis` and
`make verify-weyl-filter-flags` reproduce these gates.
For a correlated [homogeneous Boolean filter family](research/homogeneous-boolean-filter-bounds.md),
the full cost is strictly worse than the same natural-charge Gaussian
at every scale. The [global Gaussian floor](research/global-charged-gaussian-floor.md)
now proves positive complete critical pressure for **every integer length,
every real aggregate charge, and every positive scale** in the binary
Gaussian family. It combines exact charge folding, analytic tails,
small finite bases, and integer-outward certificates. The analytic floor
also extends to every power-of-two digit count.

This makes the homogeneous filter's large-scale branch absolutely
positive at every real charge. For fixed Boolean degree and bounded
charge, its small-scale bound gives positive pressure at every scale for
all sufficiently large binary prime-power lengths. This last conclusion
does not extend the correlated small-scale result to unbounded charge.
`make verify-homogeneous-boolean-filters` replays its exact identities.

The [complete one-mode Laurent theorem](research/global-charged-gaussian-floor.md)
excludes every admissible single-mode binary Laurent occupation filter
at every aggregate charge, scale, and finite degree, over any fixed
characteristic-two coefficient field. The construction retains its
power-of-two interval lengths. Inside `|J| <= 2n`, equality with the
Gaussian occurs only for a common scalar Laurent monomial. Outside that
range, relative gains can occur, but an exact loss cap and a stronger
Gaussian floor keep the full pressure positive.
`make verify-global-one-mode-filters` replays the complete certificate
chain. Independent physical-state gauges and coupled multilevel or
off-diagonal operators remain outside these exclusions.

[Extension-field operator lattices](research/identity-operator-field-flags.md)
now permit constant F4 row bases and constructive F16 endpoint lifts.
The scalar Pascal contraction also gives small exact sources for every
selected subset. A separate
[one-hole minimum tradeoff](research/identity-operator-minimum-tradeoff.md)
proves that improving that block's infinity minimum by `2*g`, for `g > 0`,
forces a zero-order loss of at least `n*P+g` powers of the local variable.
It applies when the top coefficient is a monomial, over every fixed
coefficient extension. Deeper blocks and nonmonomial top coefficients
remain open. `make verify-operator-field-flags` checks both results.

The [four-mode resonance reduction](research/binary-four-mode-resonance.md)
does reach a deeper occupation block: two equal-energy two-particle
states admit an arbitrary invertible endpoint matrix. Both complete
profiles reduce to five scalar factors and one two-by-two matrix.
An exact polynomial-degree obstruction excludes one subfamily;
[complete Laurent catalogues](research/resonant-laurent-catalogues.md)
then cover 11,796,480 admissible kernels in each of two stated coefficient
windows. Their exact signatures were searched over all 2,360 resonant
quadruples at interval 32, with no numerical gain. The windows overlap,
and the pressure conclusions concern the recorded finite grids and
refinements. `make verify-resonant-operators` replays the exact controls.

For lengths just below a power of two, the
[near-power operator reduction](research/near-power-identity-operator-reduction.md)
handles terminal modes and a split with one initial mode. It retains
`2*r+1` actual zero modes and uses a local auxiliary Gaussian only to
compute the profile; exact multiplicity division removes that factor
before physical costs are evaluated. The contraction is a scalar matrix
plus rank one. Complete interval-31 terminal and split projector screens
found no gain in their recorded grids. `make verify-near-power-operators`
checks the interpolation, auxiliary-factor removal, and original cuts.
The [fixed-deficit extension](research/fixed-deficit-identity-cores.md)
proves a corresponding mixed-end core at `n=p^k-d` in every prime
characteristic. It retains `2*r+d` actual zero modes, calculates their
profile through order `p^(r+d)`, and bounds the departure from the
scalar prime-power contraction by rank `d`. Its general gate checks
exact polynomial geometry and Gaussian-factor division; further complete
physical operator replays are stated separately.
The [arbitrary-selection operator API](research/arbitrary-one-deficit-operator-search.md)
now handles every selected set one below a prime power, with original
physical controls and exact Gaussian-factor division. Four parent-length
qutrit catalogues and four
[energy-conjugated catalogues](research/deficit-energy-conjugated-operators.md)
completed at interval 26, including interior selections. Their 9,032
kernel/geometry cases retain the Gaussian best value in the recorded
finite grids and refinements.

The [prime-field affine solver](research/prime-deficit-affine-operator-synthesis.md)
then removes the fixed-projector restriction. It gives independent
coefficients to all 19 qutrit occupation-preserving entries, with exact
endpoint equations and a checked source-precision bound. Four searches
completed 242 prescribed paths and 4,985 exact affine queries, including
nontrivial infinity flags; no returned profile improved the recorded
pressure. These paths do not exhaust every polynomial coefficient family.
`make verify-prime-deficit-affine-synthesis` replays its independent gates.

The complementary [above-power construction](research/above-power-identity-operator-reduction.md)
handles terminal directions at `n=2^k+1` by retaining the existing,
unchanged identity direction zero. Its local core and Weyl matrix are
bounded independently of the interval, with no formal physical mode
added. The [interval-33 catalogue](research/above-power-operator-searches.md)
completed all 42,560 specified projector kernels and 1,111 exact profile
classes, with no gain on the recorded grid and refinements.
`make verify-above-power-operators` checks the complete construction.

The [reflection-pair compression](research/reflection-dual-rail-compression.md)
retains one occupied mode from every reflected pair. A characteristic-two
homology calculation gives an explicit polynomial Gaussian on the
square root of the original alphabet, with complete profiles at both
cuts. Its natural pressure equals half the original pressure at twice
the scale, so the rate is unchanged. The smaller seed makes further
logical-core perturbations and independent-gauge searches accessible.

The [physical hot-state search checkpoint](research/reflection-hot-state-searches.md)
records completed positive and unrestricted height screens, all 65,536
coefficient rows at five fixed positive-height rays, and nonlinear
two-bit state penalties. None improves its Gaussian control in the
recorded scope. The [raw transition-mask screen](research/raw-core-transition-masks.md)
also retains its interrupted 167-of-615-candidate prefix explicitly;
the remaining masks are not counted as searched.

[Cold-compensated hot-state searches](research/reflection-hot-compensated-searches.md)
allow negative selected-state heights while charging every compensating
cold-bit height. The saved runs include all 65,536 coefficient rows at
seven fixed negative-height rays, 234 natural-class representatives at
interval 16, sampled rows through interval 64, and nonlinear two-bit
OR-state penalties. None supplies a gain in its stated scope. The
[separate physical-leg bulk theorem](research/reflection-hot-leg-bulk-gauges.md)
also permits different bulk heights on the four legs within an explicit
regularity cone, retaining the same bounded active core and both physical
type costs. `make verify-reflection-leg-bulks` checks its original tensors.
The [second regular cross pivot](research/reflection-hot-all-bulk-gauges.md)
removes the bulk-ordering restriction: all nonnegative physical-leg bulk
choices now use the same bounded core, subject to the stated hot-cold
regularity inequalities. Complete coefficient expansion is checked before
Smith elimination, including the large local-unit denominators of this
second pivot. `make verify-reflection-all-bulks` replays the new construction
and the shared hot-state regression.

An independent [binary crossing synthesis](research/binary-two-cut-affine-synthesis.md)
assigns free Laurent coefficients to every tensor entry while fixing the
strict permutation endpoint. For its selected row targets, a support
bound makes degree 12 after common clearing exhaustive, and all 8,925
constant F4 flags at each cut are covered. All fifty nonreference targets
are infeasible; this is an exact result for those fixed frames and gauges.
It does not exclude other target orders or other regular frames.
`make verify-binary-two-cut-lattices` checks literal minors, endpoint
coefficients, flag enumeration, and the pruning used by the complete run.

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

The [earlier 31 October research checkpoint](research/constructive-frontier-checkpoint.md)
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
A [deterministic norm-sieve experiment](research/arithmetic-norm-sieve.md)
constructs dominant relation certificates for 60-bit and 80-bit inputs and
recovers proper factors without class enumeration. Replay both certificates
and two smaller-box exhaustion controls with `make verify-arithmetic-sieve`;
success in a larger finite box is not a uniform factoring-time bound.
The direct-word certificates use the polynomial ideal-lattice composer and
a diagonal fixed by the word bound, reducing the two determinant encodings
to 520 and 500 bits. The original pair-form certificates remain replayable.
The [quantitative localization target](research/smooth-norm-localization.md)
now specifies the required short-box theorem and its conditional complete
cost. The localized smoothness statement and splitting coverage are both
unproved; the existing class-group density theorem supplies neither for free.

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
through [GitHub Discussions](https://github.com/1337isnot1337/integer-factoring/discussions),
starting with the [coordination thread](https://github.com/1337isnot1337/integer-factoring/discussions/51).
Reviewers can use the [review guide](docs/REVIEWING.md) and
[exact-head review ledger](docs/REVIEW_LEDGER.md). Document source revisions,
proof status, reproduction commands, and remaining gaps. The
[contribution record](CONTRIBUTORS.md) retains credit for useful partial and
superseded work as well as incorporated results.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
