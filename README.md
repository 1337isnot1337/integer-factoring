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
| Polynomial-time physical Ward kernel queries with arbitrary overlapping forms | Ordinary exact algorithm; numerical and formal coefficient entries, no rank bound | [One streamed Pfaffian with exact zero handling](research/ward-pfaffian-entry-queries.md) |
| Uniform eighth-order residual identity at even state energy | Ordinary exact theorem; first contact term identified | [Determinant-minor pairing and residue charge](research/homogeneous-central-eighth-jet.md) |
| Homogeneous modular duals at every label width | Ordinary uniform construction; retains exponent `1/2` | [Rectangular vacuum minor selection](research/homogeneous-central-selection.md) |
| Noncommutative 64-label constant-cut rank reduced from 31 to 28 | Exact finite construction and exhaustive search; complete profile does not improve the score | [Group-word crossings and all central shears](research/noncommutative-crossing-orbits.md) |
| Seven-label cubic cost reduced from `7+12 sqrt(7)` to 36 | Exact optimum in the stated five-dimensional `F4` space; target is about 15.673 | [Complementary output maps](research/cubic-complementary-output-maps.md) |
| Cumulative bosonic crossing profile | Ordinary exact identities and finite construction; no stronger exponent | [Moment pairing, exact rank, and both channel cuts](research/bosonic-hermite-profile.md) |
| Crossing families from binary codes | Ordinary all-size construction and one-cut bounds; complete column profile open | [Coded central quotients](research/coded-central-crossings.md) |
| Explicit complement-shell crossings | Ordinary all-size construction and exact four-state tensor reduction; tested seeds do not improve the score | [Identity Gram, sparse support, and local tensor](research/complement-shell-crossings.md) |
| Graded noncommutative Hopf crossings | Ordinary construction; full PBW-domain seeds are excluded by a rank-sum theorem | [Construction](research/noncommutative-crossing-graded-hopf.md), [four-leg degeneration](research/noncommutative-crossing-pbw-degeneration.md) |
| Boundary-section Pfaffian crossings | Ordinary all-size formula; every four-pair section is excluded over all characteristic-two fields | [Boundary Gram and exact rank-cost certificate](research/noncommutative-crossing-pfaffian-boundaries.md) |
| Class-two group-word crossings with arbitrary central controls | Ordinary all-size obstruction, including mixed lookup families, on the full label domain | [A common weighted four-leg degeneration](research/noncommutative-crossing-group-degeneration.md) |

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

The Python checks use Python 3 and its standard library. Some optional checks use NumPy for speed. Native crossing checks also require a C++20 compiler.

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

[Exact syndrome blocks](research/complement-shell-rank-blocks.md) make both full shell cut ranks accessible through 496 labels. Protected central minors rule out full central shells at any fixed length at least three when the cutoff grows. Separately, a [PBW degeneration](research/noncommutative-crossing-pbw-degeneration.md) transfers the full conservation-kernel rank-sum obstruction to the complete PBW-domain Hopf families, over every characteristic-two extension field. Selected label codes, non-graded retracts, and modified local tensors require separate analysis.

The [three-parameter Pfaffian family](research/pfaffian-ternary-edge-family.md) preserves crossing on arbitrary fixed-length occupation codes and permits linear label retracts. Six edge choices collapse to three products, giving a bounded coefficient space for complete searches. Binary and four-element-field screens exclude the recorded small spaces. [Exact residue folding](research/residue-folded-channels.md) supplies explicit incoming and outgoing factors with finite channel dimensions while retaining every high-mode contribution and the exact energy boundary.

[Boundary sections](research/noncommutative-crossing-pfaffian-boundaries.md) change the single-occupation lifts while preserving crossing. Their complete correction has an explicit formula at every size and is divisible by `p^2(p+q+r)^2`. On four pairs, every section reduces to one extra tensor direction; an integer charge certificate proves rank cost at least 120 over every characteristic-two field. Larger coupled sections remain a separate question.

[Sparse isotropic sections](research/noncommutative-crossing-pfaffian-isotropic.md) cancel the entire `p^P` endpoint coefficient for every `P >= 6` while retaining `2^P` labels. For `P >= 7`, each lift contains only four exterior monomials. At six pairs, the full family reduces to at most nine tensor directions and an explicit six-piece curve. This is an all-size exact cancellation; a rank-sum improvement has not been established.

[Optimal common retracts](research/noncommutative-crossing-optimal-retracts.md) cancel both noncrossing vertices and retain asymptotically half the full middle exterior grade. The exact maximum is `m` labels when the complement involution has a fixed monomial, and `m-1` otherwise, where `m` is its number of two-cycles. Sparse and triangular attaining maps are explicit. [Complete finite rank profiles](research/pfaffian-retract-rank-profiles.md), including extension-field searches and all 512 orientations in the ten-label case, show that the recorded choices still miss the seed threshold.

[Symplectic radical codes](research/noncommutative-crossing-symplectic-radicals.md) combine two-vertex cancellation with a constant fraction of the full label dimension at every fixed dense grade ratio. Their explicit basis vectors have polynomial support. A [Ward identity](research/noncommutative-crossing-radical-ward.md) gives an all-size factor `(q+r)^k` for codes from any product of `k` commuting forms, including overlapping forms, and an exact homogeneous dimension formula. [Complete finite profiles](research/symplectic-boundary-code-profiles.md) exclude the nine recorded binary spaces and six four-element-field spaces. The 26-label example reduces to six scalar coefficient polynomials and remains a source for tensor-power analysis.

The [26-label coefficient algebras](research/symplectic-radical-seed-algebra.md) have dimensions 45 and 11, with explicit radical filtrations and constant adapted bases. Their character minors prove that literal coefficient expansions of every normalized polynomial parameter curve have asymptotic cost at least `106.4396…^n / poly(n)`, above the required `26^(sqrt(2)n)`. The stronger common-character theorem below also covers arbitrary coefficient recombinations and rational curves.

[Exact constant module decompositions](research/symplectic-seed-tensor-modules.md) reduce the two cuts of `R(t,0)` to blocks of size at most 22 and 6. The tensor-square cuts, originally of dimension 456,976, can therefore be checked through blocks of size at most 484 and 36. Complete binary recombination has certified minimum cost in `[250528,250530]`, above the required `10049.61…`; all basis changes and rank tables are included.

[Weighted direct-summand field screens](research/symplectic-seed-tensor-field-screens.md) exclude the complete shared-parameter coefficient spaces of `R(t,0)` at powers 2, 3, and 4 over both `F2` and `F4`. The fourth-power screen covers all 22,369,621 projective vectors over `F4`. Larger omitted blocks are replaced by zero only to obtain rigorous rank lower bounds; the resulting affordable supersets already miss the crossing. A span-preservation lemma also covers scalar weighting and rational reparameterization of the same curve.

[The full weighted six-direction family](research/symplectic-seed-weighted-powers.md) also has exact constant modules of size at most 28 and 8. Complete binary searches exclude its shared-parameter square and cube, including all 268,435,455 nonzero cube coordinates. Every eligible combination was independently replayed. All 28 pairs of the eight recorded local cubic curves also miss the threshold.

A [common-character theorem](research/symplectic-seed-character-barrier.md) now excludes the entire coefficient spans of both the 24- and 26-label seeds over every characteristic-two field and at every tensor power. It also covers mixed products, independent local parameter curves, and arbitrary recombinations. Explicit row triangular blocks and a shared 43-dimensional column compression force cost greater than `d^(3/2)`. Larger codes and changed label maps remain outside this theorem.

A [direct common-functional compression](research/common-functional-rank-compressions.md) extends that exclusion to the eight-direction, 80-label five-site seed and every mixed product of the 24-, 26-, and 80-label spaces. Its reusable subspace criterion avoids constructing a coefficient algebra. The exact 80-label witnesses have compression dimensions 2,072 and 286; their product exceeds `80^3`. The checker regenerates the literal tensor and verifies both constant maps entry by entry.

[Changed four-site output sections](research/two-even-site-section-search.md) were also searched: all 728 choices of omitted, lower, or upper sections on the six pair locations fail over `F2`. Every one of their 42,536 nonzero coefficient vectors has independently checked cut ranks. The same package contains complete `F2` and `F4` profiles for the 80-label five-site code.

A [uniform nonbinary section](research/nonbinary-two-even-site-sections.md) changes all 24 output rows over `F4`, outside the fixed binary coefficient span. Its complete six-dimensional `F4` space also fails: all 4,095 nonzero vectors were tested, and the affordable vectors span four dimensions and miss the crossing. This statement concerns the specified section and field.

An [all-size small-core factorization](research/two-even-site-core-factorization.md) expresses the two-even-site family through cores on at most four sites and elementary factors elsewhere. Its complete polynomial span has at most `2P-2` directions. The exact assembler reaches 240 labels at six sites, avoiding 113 million literal contractions; all 1,023 binary combinations in its ten-dimensional span fail the seed gate. Exhaustive five-site matrix comparison and independent six-site polynomial and rank controls accompany the computation.

A [constant-compression theorem for the 240-label seed](research/two-even-site-P6-character-barrier.md) now excludes its full coefficient span over every characteristic-two field and every tensor power. The row and column compression dimensions are 15,688 and 718. Their product gives `S > 240^(37*n/25)`, with `37/25 > sqrt(2)`. Mixed tensor products of the fixed 24-, 26-, 80-, and 240-label spaces obey the same exponent obstruction. Altered output maps and larger constructions require separate analysis.

[Small Jordan modules and principal restrictions](research/two-even-site-P7-jordan-barrier.md) extend this all-field, all-power exclusion to the 672-label seven-site seed. A triangular row certificate of size 108,315 and a column compression of size 1,150 give `S > 672^(10*n/7)`. Mixed products of all five certified spaces obey exponent `10/7 > sqrt(2)`. The largest audited matrix has dimension 2,560; no full seven-site matrix is assembled.

The [eight-site extension](research/two-even-site-P8-character-barrier.md) reaches 1,792 labels using matrices of size at most 5,120. Its row and column certificates have dimensions 769,060 and 2,106. Exact integer comparisons prove `S > 1792^(99*n/70)`, where `99/70 > sqrt(2)`, over every characteristic-two field and every tensor power. The same exponent applies to all mixed products of the six certified spaces.

[Jointly nilpotent five-site row corners](research/two-even-site-row-corner-barrier.md) strengthen the eight-site result and reach **4,608 labels** at nine sites. The nine-site dimensions are 7,522,480 and 3,332, giving `S > 4608^(17*n/12)`. The seven fixed seed sizes from 24 through 4,608 can be mixed with the same exponent `17/12 > sqrt(2)`. Explicit certificate formulas hold at every larger size, but these particular witnesses fall short at ten sites; they do not exclude all larger constructions.

A [bounded Frobenius-module calculation](research/two-even-site-P10-frobenius-probe.md) raises the ten-site column certificate from 5,052 to 5,738, paired with a row certificate of 44,638,912. Its cost exponent is about 1.40448, below `sqrt(2)`, so it does not exclude the 11,520-label space. Explicit constant direct sums and coefficient-rank ceilings record the remaining compression gap without constructing a full ten-site matrix.

The [graded Ward factorization](research/noncommutative-crossing-ward-rank-descent.md) bounds restricted coefficients through smaller exterior-grade rectangles. Its [exact rectangular evaluator](research/ward-code-rank-bounds.md) handles full codes with more than `4e14` labels using explicit blocks of size at most 80 in the recorded example. The resulting upper bound still misses the target. Two actual overlapping-form codes on 6 and 34 labels also have complete binary exclusions, with all centered jets checked against the Ward bounds.

An [Euler refinement](research/noncommutative-crossing-ward-euler.md) reduces the global-form bound to two rectangles for either cut. For several forms, a [dual-incidence identity](research/noncommutative-crossing-ward-charge-filtration.md) removes the channel multiplicity from both cuts, even for arbitrary form weights in characteristic two. It uses at most `binom(k+2,2)` grade rectangles; an independent binary charge argument gives each constant factor entry as one `k` by `k` minor. The [recorded computations](research/ward-charge-bound-computations.md) lower the large global-code estimate from `rho = 0.85443` to `0.83464`; the needed value is below `1/sqrt(2)`. The binary several-form theorem also gives a joint vanishing order controlled by the form Gram matrix when outputs preserve the admissible syndromes.

An [exact label-entropy theorem](research/ward-code-charge-entropy.md) quantifies the cost of adding binary Ward forms. Across all grades, the full code dimension is at most `2^(2P-2k-q)`, where `q >= k` is the quadratic-charge rank. Equality `q=k` characterizes disjoint form supports after a row change. A positive-parameter bound also controls every individual grade, so overlap must offset its label loss with a sharper restricted-rank calculation.

A [protected-block theorem](research/ward-rectangular-bound-envelope.md) shows that the current dual-incidence **complete-rectangle estimate** has limiting `rho >= 29/40` whenever `k/P -> kappa`, `L/P -> lambda`, and `0 <= 2*kappa < lambda <= 1`. It covers every minimum implemented by that evaluator, including physical and logical dimension caps. It does not bound actual restricted ranks, additional factorizations, the separate Euler parity refinement, or regimes outside those hypotheses. The next rank calculation must retain the constant Ward maps or the actual label correlations.

An [exact stable operator interface](research/ward-stable-operator-interface.md) now retains both incidence and full-pair creation through the stable splitting. Its ancillary algebra has dimension eight; one global form has a matrix-product representation of bond at most nine. [Several forms](research/ward-form-operator-transfer.md) use the same eight ancillary channels over a square-free form-variable ring, whose scalar expansion can still require `2^k` coordinates. These are constructive identities over every characteristic-two field, not coefficient-rank complexity bounds. [Finite combined-operator calculations](research/ward-combined-operator-ranks.md) improve four recorded upper costs, with independent full-tensor checks, but remain above the seed target.

For [disjoint form supports](research/ward-disjoint-operator-transfer.md), processing each support in order reduces the scalar operator width to at most 16, independent of the number of forms. Exact grade conservation makes every fixed operator entry a single possible parameter monomial. Its coefficient can be evaluated in `O(P)` field operations with constant field-element workspace and logarithmic counters, given an ordered read-only input stream. This computes entries; the large coefficient-matrix ranks remain open.

For [overlapping supports](research/ward-interval-operator-transfer.md), the operator width is at most `8*2^w`, where `w` counts active form intervals. An explicit row reduction attains the optimal fixed-order width `max_i(rank A_{<=i} + rank A_{>=i} - k)`. An [all-size fan family](research/ward-fan-family.md) has width two for the form intervals, operator width 32, and an exact label-dimension recurrence. Its central-grade label entropy is preserved up to a polynomial factor; the first complete-rectangle cost estimate still misses the target.

A [streamed Pfaffian query algorithm](research/ward-pfaffian-entry-queries.md) evaluates physical occupation-basis entries of the combined Ward kernel for arbitrary overlapping forms. It uses `O(P*k^2+k^3)` field operations and `O(k^2)` field elements of workspace, including exact treatment of zero local factors. Finite-field interpolation then gives polynomial-time formal coefficient queries for a fixed number of parameters. The proof and implementation cover every characteristic-two field at the algebraic level; they leave the large coefficient-matrix ranks open.

For the [duplicated fan](research/ward-fan-polynomial-queries.md), the physical entry is a weighted spanning-tree polynomial with a two-state, division-free recurrence. It evaluates the entire formal polynomial in `O(P)` ring operations, including zero and nilpotent local entries. The physical operator has bond dimension two between duplicated-edge groups and at most four inside a group. An explicit sparse interface for the 84-label ten-site code reduces exact coefficient-rank calculations to 31 syndrome blocks, with largest dimension 592.

For odd-dimensional cubic radical codes, an [all-field output-map theorem](research/cubic-radical-output-maps.md) gives `S >= d(d-1)/sqrt(2)` for every syndrome-preserving left inverse and even for mixtures of such output maps. This closes that output family while leaving arbitrary left inverses open. Exact cyclic and circulant profiles, including cyclic examples through 1,023 labels, accompany the proof.

An [untouched-output gate](research/cubic-output-perturbations.md) also excludes every single-output-row perturbation of the seven-label cubic construction over all characteristic-two fields. The corresponding 1,960 triple-monomial corrections were checked literally. [Mixed-grade radical codes](research/symplectic-mixed-grade-codes.md) require a separate endpoint repair; constant row normalization and grade balancing both give exact 28-label crossings, but their recorded finite coefficient spaces do not pass the seed threshold.

An [all-size cubic output identity](research/cubic-arbitrary-output-identity.md) allows simultaneous changes to every output row while fixing the dual action on the two complementary cubic input spaces. It reduces the complete coefficient span to two fixed and two variable tensors, and gives a universal two-piece upper bound. Among 680 exact multirow maps, one new tensor has cut ranks `(16,15)` and individual cost below `7^sqrt(2)`, but its complete coefficient space still misses the seed criterion. All 10,200 nonzero coefficient rank pairs were independently replayed.

[Allowing a nonzero complementary image](research/cubic-complementary-output-maps.md) gives a five-direction formula and a three-piece cost bound `(d+1)(d+2)/2` for the stated odd-size scalar-image construction. At seven labels its `F4` optimum is exactly **36**, down from `7+12 sqrt(7)`. Complete recorded searches cover 1,195 binary maps and 56 `F4` maps. A [sumset rank-cost theorem](research/syndrome-sumset-rank-cost.md) separately excludes all tensor powers and mixtures obeying a common syndrome support law; it gives `S >= 22^n` for the seven-label diagonal family. Non-diagonal complementary images and off-syndrome corrections remain outside that theorem.

A [one-entry certificate for mixed cubic outputs](research/cubic-mixed-output-pool.md) excludes combinations of all 659 distinct affordable binary tensors from the earlier 680-map pool, even with arbitrary characteristic-two scalars on those fixed tensors. The union of all coefficient directions does contain the crossing, so newly formed low-rank combinations across spaces remain a separate possibility.

A [non-diagonal `F8` probe](research/cubic-complement-F8-probe.md) tests 428 further complementary maps, including scaled order-seven spectral projectors, cyclic and Jordan shifts, and seeded low-rank images. All 2,003,468 projective coefficient vectors were covered; every recorded affordable span misses the crossing. Independent field elimination checks all affordable vectors and rejected controls, and literal physical expansions verify six maps. This finite result does not exclude other maps, larger fields, or cross-space recombinations.

[Full mixed output spaces](research/cubic-mixed-F8-spaces.md) allow new linear combinations across maps before taking ranks. Five structured pairs, representing four distinct spaces of dimension seven or eight, have complete `F8` exclusions across 5,692,269 projective vectors. All 10,452 affordable vectors and 60 rejected controls have independently checked ranks. The scalar-family interpolation argument identifies the larger parameter spaces covered by these finite records without inferring an extension-field rank bound.

The [full cyclic scalar-output space](research/cubic-cyclic-all-power-barrier.md) is now excluded over every characteristic-two field and every tensor power: a constant row compression and a triangular column evaluation on a seven-point grid give `S >= 336^(n/2)`. The [full Jordan scalar-output space](research/cubic-jordan-all-power-barrier.md) has a separate constant restriction and triangular filtration giving `S >= 392^(n/2)`. Both proofs allow arbitrary recombinations within their respective tensor spans; directions from outside those spans require separate analysis.

[Fixed-image rank-one output mixtures](research/cubic-fixed-image-rank-one-barrier.md) obey `S >= 252^(n/2)` for every characteristic-two field and every tensor power, with arbitrary covectors and a fixed image vector. [Coordinate-image mixtures](research/cubic-coordinate-image-mixtures.md), even across all seven changed rows, obey `S >= 392^(n/2)`. These theorems close the remaining structured configurations from the mixed `F8` search.

[Varying non-coordinate image vectors](research/cubic-varying-image-affordable-escape.md) passes the necessary affordable-span test: an explicit family of 154 binary rank-one maps generates 22 individually affordable pieces that sum to the crossing. Their total cost is 49; a separate three-piece decomposition in the same space has cost **36 over `F2`**. Neither reaches the strict budget `7^sqrt(2)`, and optimal total cost in this enlarged mixed space remains open.

The [mixed-product corollary](research/mixed-common-functional-barrier.md) combines all seven certified symplectic spaces with the cyclic, Jordan, fixed-image, and coordinate-image cubic spaces. Arbitrary coefficient recombinations in any tensor product of these fixed local spaces have `S > D^(17/12)`, where `D` is the product label count. The theorem does not include arbitrary unions of distinct local spaces; the varying-image escape is outside its hypotheses.

The [uniform six-pair certificate](research/pfaffian-p6-boundary-orbit-barrier.md) now excludes every boundary section that cancels `p^6`, with binary coefficient recombinations: its rank sum is at least 834. This covers all active isometries and radical additions. Separate [finite-size gates](research/noncommutative-crossing-finite-gates.md) exclude every crossing seed on at most five labels, over every field, and every orbit-kernel seed on a non-elementary group of order at most eight.

[Midpoint blocks](research/noncommutative-crossing-abelian-midpoints.md) give an exact column-rank multiplicity for every abelian orbit kernel and exclude all cyclic groups over every field. [Involution-free diagonal minors](research/noncommutative-crossing-involution-support.md) extend the obstruction to all odd-order groups, groups with an odd-order subgroup of index two, and groups with a unique nonidentity involution. These results cover arbitrary coefficient kernels on the full group label domain.

[Complete folded Gram channels](research/folded-gram-channels.md) combine all source angles before taking rank. An explicit eight-dimensional change of basis and Frobenius recursion compute the resulting large Johnson-matrix ranks using polynomial-size blocks. [Weighted aggregation](research/folded-gram-aggregation.md) retains exact shell energy and makes dense, noncentral scans practical. The recorded 542 energy cases do not yield a sufficient rank-sum upper bound; they do not exclude smaller actual ranks.

The [particle-hole identity](research/pfaffian-particle-hole.md) expresses the complete local Pfaffian as an exterior representation of a four-dimensional product matrix, preserving the exact grade and energy conditions. A [constant 24-dimensional intertwiner](research/stable-pfaffian-modules.md) splits its remaining extension after adjoining one fixed four-dimensional factor. This gives polynomial-size blocks for rigorous complete coefficient-rank upper bounds, including bivariate jets. The current bounds miss the seed threshold; they do not give rank lower bounds. [Unequal-grade Gram channels](research/unequal-grade-folded-gram.md) also retain the exact shell energies and odd Frobenius weights; the recorded 217 cases do not improve the equal-grade bound.

[Complete coordinate searches](research/stable-pfaffian-coordinate-search.md) include every endpoint-preserving binary pencil and affine chart in the stated small spaces, using all monomial ancillary entries. A separate [boundary-block theorem](research/stable-boundary-envelope.md) proves that the six-ancilla centered full-grade upper-bound computation has asymptotic `rho >= 29/40` at every fixed dense grade ratio. This is a limitation of that computed bound; it does not bound actual coefficient ranks or restricted codes.

A larger centered-chart record contains 29 completed grade cases, reaching an upper-bound value `rho = 0.77380906` at `P=128,L=61`. The exact cost uses all 1,953 coefficient bounds and explicit matrices of size at most 240. Two proposed higher-grade cases were not completed and are explicitly marked in the archive. The needed value remains below `1/sqrt(2)`.

The [homogeneous selection theorem](research/homogeneous-central-selection.md) constructs duals of one weight at every binary label width. A large nonsingular minor exists in a homogeneous block of the proved source pairing, even when the same-weight Gram matrix is singular. The selected operator and state weights differ by `O(sqrt(h) log h)`, so the construction retains the uniform half-rate bound. Its raw-parameter Pfaffian makes complete finite profiles accessible without listing intermediate partition spaces.

[Noncommuting group words](research/noncommutative-crossing-orbits.md) give a different exact crossing family. An exhaustive screen of all `33554432` normalized functions on `F_4^2` finds the optimal constant-cut rank 28 in that central-shear family. The other pieces remain too expensive, including after all binary recombinations of the displayed winning profile. This is a finite rank improvement and a reproducible search result; the factoring score remains `1/2`.

The [cumulative bosonic profile](research/bosonic-hermite-profile.md) uses a vanishing first moment and a scalar second moment to build new exact pieces. Its scalar parameter pairing has rank `p(h)-p(h-2)`, proved using an explicit Virasoro operator. Prefix and suffix formulas give both channel factorizations without dropping contact terms. These pieces retain the leading half-rate channel scale; their complete ranks are measured separately in the finite examples.

[Coded central quotients](research/coded-central-crossings.md) let a binary code control the noncommuting group relations. A column-matroid formula bounds the constant cut, and a weighted filtration bounds every row coefficient. A [heavier central filtration](research/noncommutative-crossing-group-degeneration.md) now rules out full-domain ordered word families with arbitrary central lookup corrections, even when coefficient directions from several such families are mixed. Controls that mix central and horizontal coordinates, and arbitrary label retracts, require separate analysis.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
