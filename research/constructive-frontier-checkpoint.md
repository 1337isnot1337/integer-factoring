# Constructive frontier checkpoint: 31 October 2026

## Result of this research session

The complete deterministic factoring bound remains

\[
T(n)=\exp\!\left(O\!\left(\frac{n}{\sqrt{\log n}}\right)\right),
\qquad \alpha=\frac12.
\]

Here \(n\) is the input bit length. No strict improvement to the score is
claimed. This session ran for over three hours, beginning on 30 October 2026 at
22:53 UTC, with three parallel research workers alongside the coordinating
worker. The work
produced ordinary proofs, exact finite certificates, and reusable search
tools. No new Lean formalization was undertaken.

The required next certificate is still \(\rho<1/\sqrt2\). In the
independent-pair formulation, the equivalent target is a strict rational
instance of \(B(s)<[Z_X(s)Z_Y(s)]^R\) with \(R<\sqrt2\).

## Follow-through beginning at 04:02 UTC on 31 October

Further work continues on the same research branch, with three parallel
workers and frequent published checkpoints. The score remains \(1/2\).
The following advances are ordinary proofs and exact finite controls;
none is being presented as a strict pressure certificate.

- The formerly open mixed-field repair is now constructed explicitly,
  with a uniform simple-pole bound and independent paired-type endpoints.
  Its full-interval profile computation reduces to seven active modes at
  infinity and nine at zero. Arbitrary weights on all four extreme cores
  use ten active modes at each place.
- [Coupled constant transports](rectangular-core-transports.md),
  [polynomial shears](rectangular-polynomial-shears.md), and
  [additional quartics](rectangular-extra-quartics.md) give exact families
  beyond the first repair. Cached coordinate updates make their screens
  practical without rebuilding each full tensor.
- [Full-alphabet quartic repairs](full-alphabet-quartic-repair.md) retain
  all four extreme states on both crossing pairs. A rank-four projector
  gives the determinant \((1+\gamma)^4\), and a five-entry partial
  correction admits non-Gaussian F4 seeds. Their profiles also reduce to
  fixed-size matrices, with 23 or 29 continuous parameters.
- [Critical pressure](rectangular-critical-pressure.md) and the
  [bivariate limit](bivariate-core-pressure-limit.md) separate an actual
  defect gain from ordinary interval convergence. Exact rational LP
  certificates for the [fitted leading models](rectangular-leading-face-lp.md)
  isolate an all-low face; continuing the required entry bounds uniformly
  is a separate symbolic task. The preserved
  [all-low background](rectangular-low-background.md) is proved to be
  precisely the shorter interval of length \(M-2\).
- [All-field Gaussian endpoints](all-field-pascal-crossings.md),
  [graded powers](graded-even-power-transfer.md), and
  [truncated symmetric powers](truncated-symmetric-pascal-crossings.md)
  give further characteristic-dependent search families. The truncated
  symmetric algebra has p states per mode and permits odd-degree
  perturbations. Its interval bulk still has critical rate \(1/\sqrt2\).
- The [Boolean Cauchy operator](boolean-cauchy-smith-profile.md) has an
  all-size binomial Smith profile and an implicit one-bit inverse at a
  nearest center. The [grouped inverse theorem](grouped-cauchy-inverse-profile.md)
  accounts for the normalization precision and shows why directly
  inserting that inverse does not provide a cheaper grouped factorization.

Three [frozen rectangular search states](../evidence/rectangular-search-snapshots.json)
replayed with zero discrepancy in both complete profiles and the pressure
value: \(0.7201918203\) at eight modes, \(0.7111834960\) at sixteen,
and \(0.7084273909\) at thirty-two. These are finite search states,
not global optima or new certified bounds. They remain above the target
and do not improve the much larger interval certificate below.

The completed four-mode constant-output screen covered all 264 retained
shears at five frozen starts, with 1,320 successful profile queries and
no improvement over its control. Larger-interval and newly enlarged
construction screens are recorded separately as they complete. A failed
finite screen excludes only its stated candidates and starts.

## Quantitative progress with exact certificates

The [interval frontier](quantitative-interval-frontier.md) now has the
following exact rational upper bounds. They all remain above the target.

| Number of modes | Certified upper bound for \(\rho\) |
|---:|---:|
| 63 | \(17/24\) |
| 511 | \(99/140\) |
| 4,095 | \(577/816\) |
| 16,383 | \(3363/4756\) |

The bound at the largest mode count is approximately \(0.7071068124474\), about
\(3.13\cdot10^{-8}\) above \(1/\sqrt2\). Outward integer intervals
certify the pressure comparison; it is not inferred from floating point.
The sharp asymptotic for this displayed family is

\[
\rho_n=\frac1{\sqrt2}
 +\frac{\sqrt2}{8\pi^2n^2}
   \left(W(12n^2)^2+2W(12n^2)\right)+o(n^{-2}).
\]

This explains the positive leading gap
\((\log n)^2/(\sqrt2\pi^2n^2)\) for the stated gauges. It is a scoped
result, not a lower bound for every Gaussian construction or every crossing.
The [Pascal-complement theorem](schur-pascal-complement.md) also gives an
exact spectrum algorithm for every finite mode set, including large gaps.

## New admissible construction families

The principal constructive advance is a substantially larger set of seeds
whose endpoint and polynomiality obligations are settled.

| Family | Established feature | Remaining task |
|---|---|---|
| [Polynomial residual codes](residual-quartic-codes.md) and [arbitrary-mode extension](coded-quartic-mode-families.md) | Explicit polynomial typed crossings, with Laurent code parameters and surviving quartic terms | Improve both weighted profiles together |
| [Conjugate phase codes](conjugate-residual-phase-codes.md) | Three eight-letter polynomial crossings over the four-element field, with explicit binary descent | Find a useful physical basis and types |
| [Asymmetric involutive cores](involutive-core-crossings.md) | General polynomial family on \(3\cdot2^{n-2}\) letters; six- and twelve-letter exact examples | Make the third core symbol useful in pressure |
| [Coordinate cores](coordinate-core-crossings.md) | Genuine local poles cancel on independent crossing-pair types | Exploit the quartic interaction without suppressing its core |
| [Rectangular residual repair](rectangular-residual-repair.md) | Four-by-three extreme cores, an unrestricted two-parameter inverse, polynomial mixed examples, and a smaller height-triangular repair | Optimize the repair and two type distributions jointly |
| [Even-core postselection](even-core-postselection.md) | A strict four-letter polynomial crossing family with symbolic non-Gaussian witnesses | Improve pressure off the exactly classified Gaussian set |

The constant three-letter obstruction and these constructions are
compatible. The obstruction concerns four constant physical maps on one
residual copy. The new constructions use polynomial maps, physical height
scalings, type restrictions, or a joint output repair.

The [paired-type theorem](paired-type-crossing-pressure.md) allows the two
crossing pairs to use different alphabets and different distributions.
Tensoring with the exchanged copy restores equal physical dimensions.
Keeping both coefficient indices separate makes the costs multiply exactly.
The [field-descent theorem](fixed-field-descent.md) shows that fixed extension
fields cost only a constant factor after the full typed tensor power.

The [two-output](two-parameter-even-core.md) and
[four-boundary](four-parameter-even-core.md) classifications prevent a
particular false positive: on the nonzero Gaussian surface, the larger
boundary-parameter family is only a scalar copy of one explicit Gaussian
tensor. An improvement over a restricted parameter slice need not represent
a useful non-Gaussian interaction.

## What the numerical searches did and did not show

The numerical values below are observed search values, not certified global
minima and not exclusion theorems. Nonsmooth optimizers can stop at valuation
walls. A favorable determinant alone does not control the full profile.

- A fixed four-mode positive-energy quartic search evaluated 205,796
  candidate profiles, including joint perturbations. It did not improve its
  displayed Gaussian control, approximately \(0.731080936\).
- The three-letter no-bulk monomial screen covered 125 parameter triples,
  two gate choices, and two cut orientations: 500 optimized candidates.
  Its best observed value was approximately \(0.9244445\).
- A central-window postselection screen covered 252 cases. Its best control
  was Gaussian, approximately \(0.8112785\); the nonzero parameter cases
  tested did not improve that control.
- All 243 constant phase masks for the chosen one-copy four-element-field
  family were screened at the stated natural gauges. They did not improve
  the filled-low background, approximately \(0.7602278\).
- The coordinate twelve-letter search reached approximately \(0.74743\)
  with a common type and \(0.73784\) with independent paired types.
  The common-type run suppressed the third core symbol. The paired run
  approached a boundary on one crossing pair; its value is also only an
  exploratory floating result.
- Full rectangular inversion initially damaged the useful background.
  The smaller two-entry repair recovered an observed value near
  \(0.75583\), near a boundary whose exact restriction is Gaussian. The
  no-bulk partial tensor for \(S=\{0,1\}\), \(M=2\) is already
  non-Gaussian: an accessible Wick residual is exactly \((1+w)^4\).
- Twenty-five polynomial mixed rectangular monomial seeds were tested at
  seven inherited/core gauges and both cuts. None improved the undeformed
  control near \(0.81705\).
- The even-core screen covered 159 polynomial endpoint perturbations and
  refined 21 of them. None improved the exceptional Gaussian control.
  The separate input-bias and four-element-field constant screens also
  gave no improvement; their scopes are recorded in the
  [search summary](../evidence/even-core-search-summary.json).

An obsolete experimental subclass briefly applied the scalar denominator
restoration twice. Those outputs were discarded and the affected screens
were rerun. They are not used by any theorem or certificate. The published
native tools apply the correction once and test it independently.

Saved [square-family starts](../evidence/constructive-frontier-search-starts.json)
include the physical gauges, both types, constructor parameters, and seed
hashes. They can be reconstructed and replayed with the separate rectangular
engine. The [rectangular starts](../evidence/rectangular-search-starts.json)
likewise retain seven independently replayed states and the scopes of 1,540
bounded profile-query attempts. These are floating consistency checks, not
exact bounds:

```sh
python tools/replay_constructive_search_starts.py
```

## Reproduction and next mathematical obligations

The ordinary `make verify` suite passed during the session. Every later
construction checker was also replayed independently. The new search tools
were compared with exact rational-weight profile calculations, and their
failure behavior and denominator normalization were checked. No new Lean
build is required for these additions.

```sh
make verify
make verify-fixed-modes
make verify-crossing-interpolation
make verify-native-smith verify-native-rectangular
```

The native targets are optional and require a C++17 compiler. Their Python
wrappers and finite controls use the standard library. Numerical optimization
programs may additionally use NumPy and SciPy.

The original checkpoint proposed the mixed coefficients
\(\alpha=1,\beta=\omega\) as an untested extension. That construction,
its denominator proof, and the four-element-field rectangular engine are
now complete. The current targets are a strict finite witness in the
expanded repair families, a negative finite limiting pressure in a family
that changes the preserved Gaussian background, and a scalar recovery
method that can use the improved operator precision without recreating
the original grouping cost. Fixed-field descent handles the eventual
binary transfer once a strict crossing certificate exists.

For any numerical candidate, the next gate is an exact weighted
profile or coefficient-rank bound and a strict rational pressure comparison.
Only after that comparison passes the threshold does the effective
crossing-to-factoring transfer change the score.
