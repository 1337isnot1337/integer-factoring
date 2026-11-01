# A four-axis Gaussian profile floor from shared endpoint minors

## Result and scope

Every homogeneous polynomial Gaussian kernel on the sixteen vertices of a
four-dimensional cube, with the antipodal endpoint, has

```math
\rho_{\mathrm{profile}}(p,t)>\frac12
\qquad(0<p<1,\ 0<t<1)
```

for the common-scale, unweighted Smith-profile criterion defined below.
Its infimum is therefore at least `1/2`. This criterion cannot supply a
four-axis saving exponent above one half.

There is no Cayley, translation-invariance, monomial, coordinate-symmetry,
or pure-singleton assumption. The proof uses one shared four-by-four
endpoint minor for each pair of axes, followed by convexity and the binary
entropy identity. It is a floor for this rank-estimation method, not for
actual coefficient ranks or integer factoring. Weighted physical criteria,
different restrictions, nonhomogeneous kernels, and non-Gaussian tensors
remain outside its scope.

## 1. Kernel and numerical criterion

Work over a field `F` of characteristic two. Index a sixteen-by-sixteen
alternating polynomial matrix `K(x_0,x_1,x_2,x_3)` by the binary vertices
`v in {0,1}^4`. Every nonzero entry has the same positive homogeneous
degree `D`. Let `bar v=v xor (1,1,1,1)` and require the endpoint

```math
K_{u,v}(1,1,1,1)=[v=\bar u].                     \tag{1}
```

For an axis `a`, let `B_a` be the eight-by-eight crossblock with rows
`u_a=0` and columns `v_a=1`. Its determinant is nonzero by (1). Treat the
other three variables as elements of their rational function field, and
write its ordered `x_a`-adic Smith valuations as

```math
0\le s_{a,1}\le\cdots\le s_{a,8},\qquad
\Delta_4(a)=\sum_{j=1}^4s_{a,j},\qquad
Z_a(t)=\prod_{j=1}^8(1+t^{s_{a,j}}).              \tag{2}
```

These are generic valuations. Setting the other variables to one can
introduce cancellations and does not define these orders.

The physical Gaussian tensor has entries `Pf K[S,S]`. Its full cut
factorization has the subset sums of the one-particle orders, as in the
[Gaussian cube construction](cube-homogeneous-gaussian.md). Restricting
each of the sixteen physical words to a common occupied fraction `p`
makes every nonzero entry at word length `m` homogeneous of degree
`8Dpm`. The corresponding common-scale profile expression is

```math
\rho_{\mathrm{profile}}(p,t)=
\frac{\sum_{a=0}^3\log Z_a(t)-8Dp\log t}{32H(p)},
\qquad H(p)=-p\log p-(1-p)\log(1-p).             \tag{3}
```

Equation (3) is the numerical estimate being bounded. No effective
multivariate Smith factorization is inferred merely from the existence
of the valuations.

## 2. A common minor belongs to two cuts

Fix distinct axes `a,b`. Define

```math
V_{ab}=\{v\in\{0,1\}^4:v_a=v_b=0\}.
```

List its four vertices in increasing binary order, and use their antipodes
in the **same** order for the column list. Form

```math
M_{ab}=(K_{u,\bar v})_{u,v\in V_{ab}},\qquad
F_{ab}=\det M_{ab}.                              \tag{4}
```

This very same matrix is a minor of both `B_a` and `B_b`: every row has
both coordinates zero, and every column has both coordinates one. An
ambient row or column ordering may introduce a permutation, which does
not change a valuation. Its endpoint is the identity matrix by (1), so
`F_ab(1,1,1,1)=1`. In particular, the minor is not the zero polynomial.

Every term in its determinant has total degree `4D`; hence `F_ab` is a
nonzero homogeneous polynomial of that degree. If `v_a` denotes coordinate
valuation, then

```math
v_a(F_{ab})+v_b(F_{ab})\le4D.                    \tag{5}
```

Indeed, every monomial is divisible by both indicated coordinate powers,
whose combined degree cannot exceed its total degree. This remains valid
when determinant terms cancel.

The determinantal-divisor characterization of Smith orders says that
`Delta_4(a)` is the least `x_a` valuation among all four-by-four minors
of `B_a`. Applying it to the shared minor in (4) gives

```math
\Delta_4(a)+\Delta_4(b)
\le v_a(F_{ab})+v_b(F_{ab})\le4D.                \tag{6}
```

Sum (6) over all six unordered axis pairs. Each axis occurs three times,
so

```math
\sum_{a=0}^3\Delta_4(a)\le8D.                  \tag{7}
```

The same polynomial minor is essential. Independent specializations of
different minors would not justify the shared degree budget (5).

## 3. Convexity and entropy give the half-score obstruction

For fixed `0<t<1`, the function

```math
\phi(s)=\log(1+t^s),\qquad s\ge0,
```

is decreasing and convex. The sixteen orders with `j<=4`, across all four
cuts, have total at most `8D` by (7). Jensen's inequality and monotonicity
therefore give

```math
\sum_{a=0}^3\sum_{j=1}^4\phi(s_{a,j})
\ge16\phi\!\left(\frac{\sum_a\Delta_4(a)}{16}\right)
\ge16\log(1+t^{D/2}).                           \tag{8}
```

Every one of the other sixteen terms is strictly positive. Thus

```math
\sum_a\log Z_a(t)>16\log(1+t^{D/2}).             \tag{9}
```

Set `u=t^(D/2)`. The elementary entropy inequality

```math
\log(1+u)-p\log u\ge H(p)                       \tag{10}
```

follows, for example, by weighted AM–GM applied to `u/p` and `1/(1-p)`
with weights `p` and `1-p`. Substituting (9) into (3) now yields

```math
\rho_{\mathrm{profile}}(p,t)
>\frac{16(\log(1+u)-p\log u)}{32H(p)}
\ge\frac12.                                    \tag{11}
```

The pointwise strict inequality alone need not give a uniform gap above
one half. The conclusion needed here is that no choice of these parameters,
or of a homogeneous kernel in this class, can give `rho<1/2` by (3).

## Checks and attribution

Run:

```sh
python3 tools/check_four_axis_gaussian_profile_floor.py \
  --output evidence/four-axis-gaussian-profile-floor.json
```

The checker uses exact four-variable polynomials over `F2`. It checks the
shared row/antipodal-column indexing at every axis pair, computes each shared
determinant independently by its permutation formula, and verifies its
endpoint, homogeneous degree, and two-coordinate valuation budget on small
kernels, including kernels that are not Cayley. It also replays the
triangular augmentation minors that first suggested the argument, and
checks the cleared-denominator entropy inequality on rational parameters.

These are finite controls of the identities and indexing. The general
conclusion follows from the shared-minor, convexity, and entropy arguments
above. The non-Cayley shared-minor strengthening was identified during the
parallel research review. The Gaussian profile setup comes from the cited
cube note. This work used substantial OpenAI assistant support; no Lean
verification or improved achieved factoring exponent is claimed.
