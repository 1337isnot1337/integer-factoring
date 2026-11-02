# Color parity obstructs large one-hot cube profiles

## Result and scope

Let a polynomial cube tensor over any field have one vacuum and `d=q-1`
nonzero colors. Assume its entries have homogeneous degree three times
occupied vertex count, its `(1,1,1)` endpoint is antipodal equality, and
an entry vanishes whenever any individual color occurs an odd number of
times. For every `q>=34`, its common-scale, unweighted three-axis
Smith-profile method with a common physical type has pressure

```math
\boxed{\rho_{\mathrm{prof}}\ge\frac34>\frac1{\sqrt3}.}
```

Gaussian factorization, color symmetry and axis symmetry are unnecessary.
This includes the proposed one-hot color-class Gaussian at `q=256` and
all larger alphabets. Large constructions retaining separate color parity
are stopped under this criterion. Smaller alphabets, types differing
between vertices, weighted gauges, actual coefficient ranks and tensors
breaking a displayed hypothesis remain outside the theorem. The complete
factoring exponent remains `1/2`.

The polynomial inequality is the existing
[scalar color-clone certificate](cube-scalar-color-clone-floor.md), source
`0c78cd26c5a7eefad1e107077b6385cc2bbdd9fb`. The contribution here is the
extension from scalar clones to arbitrary parity-preserving tensors by
using **odd-color sectors**, rather than exact color multiplicities. The
worker-4 [construction and specialized controls](color-class-cube-profile.md)
are retained from `1772970a7dbaab83eeb0607014a1111e46718d86`, with their
correction history and original source credit.

## 1. Closed parity sectors in every cut

For a four-vertex side word, define `C` as the set of nonzero colors
occurring an odd number of times. Individual color parity forces its
column word to have exactly the same `C`. Consequently every cut matrix
is a direct sum over `C`, where `|C|<=4`. Its endpoint is a permutation
matrix on each sector, so all sector Smith orders are finite.

Write `m_(a,C)` for the smallest order in sector `C` at axis `a`. Choose
`r=|C|` antipodal pairs, assign the colors of `C` once each to those pairs,
and set the remaining pairs to vacuum. The same physical tensor entry
`f_C(x,y,z)` lies in sector `C` at every cut. It is nonzero, has endpoint
one and degree `6r`. Hence

```math
\sum_a m_{a,C}\le\sum_a\operatorname{ord}_a f_C\le6r. \tag{1}
```

The first inequality uses the minimum entry valuation characterization of
the first Smith order. The second holds because the three coordinate
powers divide a nonzero polynomial of degree `6r`. For `C` empty the
all-vacuum entry is the constant one and all three minima are zero.

Retain one channel per sector in each axis partition function `Z_a(t)`.
Convexity of log-sum-exp, followed by (1) and `0<t<1`, gives

```math
\frac13\sum_a\log Z_a(t)
 \ge\log\left(\sum_C t^{\frac13\sum_a m_{a,C}}\right)
 \ge\log Z_d(t^2),\qquad
 Z_d(u)=\sum_{r=0}^4\binom dr u^r.                    \tag{2}
```

Use `binom(d,r)=0` for `r>d`. Each retained channel belongs to a disjoint
summand; compatible Smith bases across cuts are not needed.

When `r=4`, each sector has dimension 24 and consists of all permutations
of its four colors. For `r<4`, exact color multiplicities can vary within
a sector; this is why the parity formulation applies beyond scalar clones.
Only the common entry and sector separation are used, not a common full
minor in all cuts.

## 2. Reuse of the existing exact polynomial certificate

The [existing certificate](cube-scalar-color-clone-floor.md#2-an-exact-all-size-polynomial-certificate)
proves, for every integer `d>=33`,

```math
Z_d(u)>(1+d u^{4/3})^3\qquad(0<u<1).                  \tag{3}
```

Its ten transformed coefficients have degree at most four in `d`.
Their Newton expansions at `d=33` have 50 nonnegative integer
coefficients, with strictly positive constants. The unchanged
[checker](../tools/check_cube_scalar_color_clone_floor.py) reconstructs
that table and additional exact controls; the
[receipt](../evidence/cube-scalar-color-clone-floor.json) retains its scope.
The polynomial identity itself is independent of a tensor construction.
Equation (2) supplies the new hypothesis needed to reuse it here.

Let the common physical type have entropy `H(p)>0` and mean occupancy
`E`. Define its pointwise profile objective

```math
\Phi_p(t)=
 \frac{\frac13\sum_a\log Z_a(t)-8E\log t}{4H(p)}.       \tag{4}
```

This is the common-scale objective minimized by the unweighted profile
method. The physical-letter generating function is `W(v)=1+d v`, and
`H(p)+E log(v)<=log W(v)` for any color distribution. Set
`u=t^2` and `v=u^(4/3)=t^(8/3)`. Equations (2)--(3) imply

```math
\frac13\sum_a\log Z_a(t)
 >3\log W(v)\ge3H(p)+8E\log t.
```

Thus `Phi_p(t)>3/4` for every `p,t` in the stated range; taking the
infimum preserves the floor `rho_prof>=3/4`. Finally `3*3^2>4^2`
proves `3/4>1/sqrt(3)`. No numerical optimization or large matrix is used.

## 3. Exact controls and evidence limits

```sh
make verify-color-distinct-sector
```

The target replays the existing all-size scalar polynomial certificate
and the [new structural controls](../tools/check_color_distinct_sector.py).
The latter verify 1728 Gaussian sector endpoint identities, 3816
zero-coupling identities and 108 complementary-minor entries symbolically.
For this Gaussian specialization the four-color minimum is exactly eight
in every axis and the complementary minor is `diag(1,x^2y^2z^2 I_d)`.
These are finite controls, rather than an exhaustive check of arbitrary
parity-preserving tensors; the general theorem is the ordinary proof above.

An asymmetric monomial tensor `T(a)=x^(3w(a))*delta_antipodal(a)` at
`q=5,t=1/2` checks why an individual-axis partition bound can fail even
though the three-axis geometric-mean bound holds. The
[structural receipt](../evidence/color-distinct-sector.json) records this.

The inherited source hashes match. Selected independent regressions
replayed the complete binary specialized profile, the ternary profile
through order six, and saved graded comparisons. Larger specialized
screens were not independently rerun. The historical `q4_full_profile`
key contains only 250 of 256 orders through 39. Specialization bounds
cumulative counts through each order, not individual layer multiplicities.
No specialized rank extrapolation is used by this theorem.

No new Lean or independent mathematical review is claimed. Full default
verification and the integration follow-up are recorded at their exact
sources in the review ledger.
