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
| Uniform eighth-order residual identity at even state energy | Ordinary exact theorem; first contact term identified | [Determinant-minor pairing and residue charge](research/homogeneous-central-eighth-jet.md) |
| Homogeneous modular duals at every label width | Ordinary uniform construction; retains exponent `1/2` | [Rectangular vacuum minor selection](research/homogeneous-central-selection.md) |
| Noncommutative 64-label constant-cut rank reduced from 31 to 28 | Exact finite construction and exhaustive search; complete profile does not improve the score | [Group-word crossings and all central shears](research/noncommutative-crossing-orbits.md) |
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

The [graded Ward factorization](research/noncommutative-crossing-ward-rank-descent.md) bounds restricted coefficients through smaller exterior-grade rectangles. Its [exact rectangular evaluator](research/ward-code-rank-bounds.md) handles full codes with more than `4e14` labels using explicit blocks of size at most 80 in the recorded example. The resulting upper bound still misses the target. Two actual overlapping-form codes on 6 and 34 labels also have complete binary exclusions, with all centered jets checked against the Ward bounds.

An [Euler refinement](research/noncommutative-crossing-ward-euler.md) reduces the global-form bound to two rectangles for either cut. For several forms, a [dual-incidence identity](research/noncommutative-crossing-ward-charge-filtration.md) removes the channel multiplicity from both cuts, even for arbitrary form weights in characteristic two. It uses at most `binom(k+2,2)` grade rectangles; an independent binary charge argument gives each constant factor entry as one `k` by `k` minor. The [recorded computations](research/ward-charge-bound-computations.md) lower the large global-code estimate from `rho = 0.85443` to `0.83464`; the needed value is below `1/sqrt(2)`. The binary several-form theorem also gives a joint vanishing order controlled by the form Gram matrix when outputs preserve the admissible syndromes.

An [exact label-entropy theorem](research/ward-code-charge-entropy.md) quantifies the cost of adding binary Ward forms. Across all grades, the full code dimension is at most `2^(2P-2k-q)`, where `q >= k` is the quadratic-charge rank. Equality `q=k` characterizes disjoint form supports after a row change. A positive-parameter bound also controls every individual grade, so overlap must offset its label loss with a sharper restricted-rank calculation.

A [protected-block theorem](research/ward-rectangular-bound-envelope.md) shows that the current dual-incidence **complete-rectangle estimate** has limiting `rho >= 29/40` whenever `k/P -> kappa`, `L/P -> lambda`, and `0 <= 2*kappa < lambda <= 1`. It covers every minimum implemented by that evaluator, including physical and logical dimension caps. It does not bound actual restricted ranks, additional factorizations, the separate Euler parity refinement, or regimes outside those hypotheses. The next rank calculation must retain the constant Ward maps or the actual label correlations.

For odd-dimensional cubic radical codes, an [all-field output-map theorem](research/cubic-radical-output-maps.md) gives `S >= d(d-1)/sqrt(2)` for every syndrome-preserving left inverse and even for mixtures of such output maps. This closes that output family while leaving arbitrary left inverses open. Exact cyclic and circulant profiles, including cyclic examples through 1,023 labels, accompany the proof.

An [untouched-output gate](research/cubic-output-perturbations.md) also excludes every single-output-row perturbation of the seven-label cubic construction over all characteristic-two fields. The corresponding 1,960 triple-monomial corrections were checked literally. [Mixed-grade radical codes](research/symplectic-mixed-grade-codes.md) require a separate endpoint repair; constant row normalization and grade balancing both give exact 28-label crossings, but their recorded finite coefficient spaces do not pass the seed threshold.

The [uniform six-pair certificate](research/pfaffian-p6-boundary-orbit-barrier.md) now excludes every boundary section that cancels `p^6`, with binary coefficient recombinations: its rank sum is at least 834. This covers all active isometries and radical additions. Separate [finite-size gates](research/noncommutative-crossing-finite-gates.md) exclude every crossing seed on at most five labels, over every field, and every orbit-kernel seed on a non-elementary group of order at most eight.

[Midpoint blocks](research/noncommutative-crossing-abelian-midpoints.md) give an exact column-rank multiplicity for every abelian orbit kernel and exclude all cyclic groups over every field. [Involution-free diagonal minors](research/noncommutative-crossing-involution-support.md) extend the obstruction to all odd-order groups, groups with an odd-order subgroup of index two, and groups with a unique nonidentity involution. These results cover arbitrary coefficient kernels on the full group label domain.

[Complete folded Gram channels](research/folded-gram-channels.md) combine all source angles before taking rank. An explicit eight-dimensional change of basis and Frobenius recursion compute the resulting large Johnson-matrix ranks using polynomial-size blocks. [Weighted aggregation](research/folded-gram-aggregation.md) retains exact shell energy and makes dense, noncentral scans practical. The recorded 542 energy cases do not yield a sufficient rank-sum upper bound; they do not exclude smaller actual ranks.

The [particle-hole identity](research/pfaffian-particle-hole.md) expresses the complete local Pfaffian as an exterior representation of a four-dimensional product matrix, preserving the exact grade and energy conditions. A [constant 24-dimensional intertwiner](research/stable-pfaffian-modules.md) splits its remaining extension after adjoining one fixed four-dimensional factor. This gives polynomial-size blocks for rigorous complete coefficient-rank upper bounds, including bivariate jets. The current bounds miss the seed threshold; they do not give rank lower bounds. [Unequal-grade Gram channels](research/unequal-grade-folded-gram.md) also retain the exact shell energies and odd Frobenius weights; the recorded 217 cases do not improve the equal-grade bound.

[Complete coordinate searches](research/stable-pfaffian-coordinate-search.md) include every endpoint-preserving binary pencil and affine chart in the stated small spaces, using all monomial ancillary entries. A separate [boundary-block theorem](research/stable-boundary-envelope.md) proves that the six-ancilla centered full-grade upper-bound computation has asymptotic `rho >= 29/40` at every fixed dense grade ratio. This is a limitation of that computed bound; it does not bound actual coefficient ranks or restricted codes.

The [homogeneous selection theorem](research/homogeneous-central-selection.md) constructs duals of one weight at every binary label width. A large nonsingular minor exists in a homogeneous block of the proved source pairing, even when the same-weight Gram matrix is singular. The selected operator and state weights differ by `O(sqrt(h) log h)`, so the construction retains the uniform half-rate bound. Its raw-parameter Pfaffian makes complete finite profiles accessible without listing intermediate partition spaces.

[Noncommuting group words](research/noncommutative-crossing-orbits.md) give a different exact crossing family. An exhaustive screen of all `33554432` normalized functions on `F_4^2` finds the optimal constant-cut rank 28 in that central-shear family. The other pieces remain too expensive, including after all binary recombinations of the displayed winning profile. This is a finite rank improvement and a reproducible search result; the factoring score remains `1/2`.

The [cumulative bosonic profile](research/bosonic-hermite-profile.md) uses a vanishing first moment and a scalar second moment to build new exact pieces. Its scalar parameter pairing has rank `p(h)-p(h-2)`, proved using an explicit Virasoro operator. Prefix and suffix formulas give both channel factorizations without dropping contact terms. These pieces retain the leading half-rate channel scale; their complete ranks are measured separately in the finite examples.

[Coded central quotients](research/coded-central-crossings.md) let a binary code control the noncommuting group relations. A column-matroid formula bounds the constant cut, and a weighted filtration bounds every row coefficient. A [heavier central filtration](research/noncommutative-crossing-group-degeneration.md) now rules out full-domain ordered word families with arbitrary central lookup corrections, even when coefficient directions from several such families are mixed. Controls that mix central and horizontal coordinates, and arbitrary label retracts, require separate analysis.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
