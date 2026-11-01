# Symplectic blocks for symmetric cube cut operators

This note gives an exact reduction for multiplication by a symmetric
element of a characteristic-two exterior algebra. With `r` generators,
the full multiplication matrix has order `2^r`. One constant change of
basis splits it into at most `floor(r/2)+1` distinct blocks, each of
order at most `2^(floor(r/2)+1)`, with explicit multiplicities. The
change of basis works for arbitrary symmetric coefficients over every
commutative characteristic-two ring.

Applied to the natural higher-dimensional version of the
[homogeneous Cayley cube kernel](cube-homogeneous-gaussian.md), this
makes complete one-variable Smith profiles inexpensive enough to
compute through twelve axes. Exact all-scale certificates then prove
that this particular family's common-scale profile method misses
`1/sqrt(d)` in every dimension tested, including nine through twelve.
This is a scoped exclusion, not a bound on all cube tensors or on
factoring. The proved factoring score remains `1/2`.

The computation uses a specialization only for rejection. It does not
identify specialized orders with generic multivariable orders. The
direction of the resulting inequality is proved in Section 4.

## 1. A constant canonical form for all symmetric exterior elements

Let

\[
 A_r=R[\varepsilon_1,\ldots,\varepsilon_r]/
           (\varepsilon_1^2,\ldots,\varepsilon_r^2),
 \qquad \mathrm{char}(R)=2,
\]

and let `e_j` be the elementary symmetric polynomial of degree `j`
in the generators. Put `e_0=1`, `e_j=0` for `j>r`, and

\[
 u=e_1,\qquad \Omega=e_2=\sum_{i<j}\varepsilon_i\varepsilon_j.
\]

Every elementary form of even degree is a Pfaffian divided power of
`Omega`. Explicitly, its coefficient on any `2j` distinct generators
is the Pfaffian of the all-ones off-diagonal alternating matrix on
those generators. There are `(2j-1)!!` perfect matchings, an odd number,
so that coefficient is one. Thus the divided power is exactly `e_(2j)`.
These divided powers commute with a linear change of generators: expand
the Pfaffian by matchings and then expand every new generator linearly;
terms with a repeated generator vanish and the remaining minors give
the same coefficients in either order.

Also, directly counting the possible removed generators gives

\[
 u e_{2j}=(2j+1)e_{2j+1}=e_{2j+1}.                       \tag{1}
\]

Over `F_2`, the matrix of `Omega` is `I+J`, where `J` is all ones.
If `r=2s` is even it is nonsingular, with inverse itself. Choose a
symplectic basis `v_1,w_1,...,v_s,w_s` with `v_1=u`. If `r=2s+1` is
odd, its image is the even-coordinate hyperplane and its rank is `2s`;
the vector `u` is outside that hyperplane. Choose symplectic pairs
there and keep `u` as the extra generator. These choices are ordinary
binary Gaussian elimination and give a matrix in `GL_r(F_2)`, hence a
valid change of generators over every such `R`.

Writing `z_i=v_i w_i`, the divided-power identity gives

\[
 e_{2j}=e_j(z_1,\ldots,z_s),\qquad
 e_{2j+1}=u e_j(z_1,\ldots,z_s).                         \tag{2}
\]

Each `z_i` squares to zero. In the even case `u=v_1`, so `u z_1=0`.
The [implementation](../tools/cube_symmetric_exterior_blocks.py)
constructs the binary basis, rather than assuming a supplied change of
coordinates.

## 2. Exact blocks and their multiplicities

Consider arbitrary coefficients `a_0,...,a_r in R` and the operator of
multiplication by

\[
 f=\sum_{j=0}^r a_j e_j.                                 \tag{3}
\]

For one symplectic pair, multiplication by `z_i=v_i w_i` preserves
the two-dimensional active space `span{1,z_i}` and kills each of the
two spectator lines `span{v_i}`, `span{w_i}`. Splitting pairs according
to these invariant spaces gives the following complete direct sums.

### Odd number of generators

For `r=2s+1`, choose `k` active pairs and choose one of the two spectator
lines on every other pair. The resulting block is multiplication on
`R[u,z_1,...,z_k]/(u^2,z_1^2,...,z_k^2)` by

\[
 f_k=\sum_{j=0}^k a_{2j}e_j(z)
       +u\sum_{j=0}^k a_{2j+1}e_j(z).                    \tag{4}
\]

Its order and multiplicity are

\[
 2^{k+1},\qquad \binom sk 2^{s-k},\qquad 0\le k\le s.   \tag{5}
\]

### Even number of generators

For `r=2s`, retain the first pair in full because `u=v_1`. Split the
other `s-1` pairs. On a block with `k` active remaining pairs, write
`v=v_1,w=w_1`; its multiplier is

\[
 f_k=\sum_{j=0}^k a_{2j}e_j(z)
   +v\sum_{j=0}^k a_{2j+1}e_j(z)
   +vw\sum_{j=0}^k a_{2j+2}e_j(z).                       \tag{6}
\]

Its order and multiplicity are

\[
 2^{k+2},\qquad \binom{s-1}k2^{s-1-k},
 \qquad0\le k\le s-1.                                   \tag{7}
\]

The dimensions in (5) sum to `2^(2s+1)` and those in (7) to `2^(2s)`.
Every summand is retained; there is no rank truncation or discarded
state. The decomposition is conjugacy by a constant binary matrix.
It consequently preserves complete Smith factors over any scalar
discrete valuation ring, as well as ordinary ranks after specialization.

For coefficient polynomials modulo `x^p`, direct elimination on these
distinct blocks costs `2^(3r/2) poly(r,p)` bit operations under a naive
cubic matrix bound. The multiplicities are integers of `O(r)` bits;
the repeated blocks need not be recomputed. The unreduced cubic bound
would have the factor `2^(3r)`. This comparison concerns the explicit
matrix problem, not the global factoring runtime.

## 3. The natural homogeneous multi-axis Cayley family

For `d>=2`, put `D=lcm(1,...,d)` and work in the group algebra of
`(Z/2)^d`, with `epsilon_a=1+g_a`. Set

\[
 c_S=\prod_{a\in S}x_a^{D/|S|}\quad(S\ne\varnothing),
 \qquad c_\varnothing=\sum_{S\ne\varnothing}c_S,
 \qquad f=\sum_S c_S\varepsilon_S.                       \tag{8}
\]

Write `f=sum_h k_h g^h` and define the alternating Cayley kernel
`K_(v,w)=k_(v xor w)`. The chosen empty coefficient makes `k_0=0`.
At `(1,...,1)` every augmentation coefficient is one, so
`f=product_a(1+epsilon_a)=product_a g_a`. Thus the endpoint is precisely
the antipodal matching. Every entry is homogeneous of degree `D`.
The principal-Pfaffian tensor therefore has the required antipodal
equality endpoint.

At axis `a`, its one-particle crossblock is multiplication by
`sum_T c_(T union {a}) epsilon_T` on the remaining `r=d-1` generators.
Set every variable other than `x_a=x` to one. The resulting multiplier is

\[
 F_d(x)=\sum_{j=0}^{d-1}x^{D/(j+1)}e_j.                   \tag{9}
\]

This is exactly (3), so (4)--(7) compute the complete specialized
Smith profile. Since every nonconstant exterior monomial squares to
zero in characteristic two,

\[
 F_d(x)^2=x^{2D}I.                                       \tag{10}
\]

The multiplication matrix is triangular with diagonal `x^D`, so its
determinant is `x^(mD)`, where `m=2^(d-1)`. Its orders are nonnegative,
sum to `mD`, are at most `2D`, and pair to `2D`. The last statement
follows by comparing the orders of the inverse with
`F_d^(-1)=x^(-2D)F_d`.

Thus absolute series precision `2D+1` determines every order. The
checker performs elimination at that precision, retains every pivot,
and rejects a missing pivot rather than silently interpreting it as a
large order. At `d=9` the full matrix has order 256, while its largest
block has order 32; at `d=12` the corresponding orders are 2048 and 64.

## 4. Why specialization proves a rejection, not a favorable certificate

Let `sigma_(a,1)<=...<=sigma_(a,m)` be the generic orders at axis `a`
over the rational function field in the other variables, and let
`s_1<=...<=s_m` be the computed orders of (9). A specialization can
annihilate a minor's first coefficient, but cannot create a coefficient
at a smaller power of `x`. The determinantal-ideal characterization
therefore gives

\[
 \sum_{i=1}^k\sigma_{a,i}\le\sum_{i=1}^k s_i
 \quad(1\le k\le m).                                    \tag{11}
\]

Both full determinants equal `x^(mD)`, so equality holds at `k=m`.
For `0<t<1`, the function `s -> log(1+t^s)` is convex and decreasing.
The prefix inequalities imply

\[
 \sum_i\log(1+t^{\sigma_{a,i}})
 \ge\sum_i\log(1+t^{s_i}).                               \tag{12}
\]

One may prove (12) by successive transfers from larger coordinates to
smaller coordinates while preserving the total; these transfers reduce
a convex sum. Thus the specialized profile is an optimistic lower
bound for the generic counting expression. The same specialized list
occurs at every axis by the symmetry of (8).

The full Gaussian cut profile consists of all subset sums of its
one-particle orders. With common occupied fraction `p` on the binary
physical letters, homogeneity puts the total coefficient degree at
`m D p` times the word length. Set `x=-D log t` and
`H(p)=-p log p-(1-p)log(1-p)`. The declared common-scale profile pressure
is consequently bounded below by

\[
 R_d(p,x)=
 \frac{m^{-1}\sum_i\log(1+e^{-x s_i/D})+xp/d}{H(p)}.
                                                               \tag{13}
\]

Only `0<p<1` matters; the endpoint types have zero alphabet entropy.
Equation (12) is essential: using the specialized orders to claim an
upper bound on the generic pressure would reverse the valid inference.

## 5. Exact all-scale floors through twelve axes

The complete profiles from (4)--(9) give these rigorous floors for
(13), for every `x>=0` and every positive-entropy binary type:

| Axes `d` | Strict profile floor | Comparison with the half-score threshold |
| --- | --- | --- |
| 2 | `3/4` | `2*3^2 > 4^2` |
| 3 | `3/5` | `3*3^2 > 5^2` |
| 4 | `11/20` | `4*11^2 > 20^2` |
| 5 | `1/2` | `5 > 2^2` |
| 6 | `19/40` | `6*19^2 > 40^2` |
| 7 | `9/20` | `7*9^2 > 20^2` |
| 8 | `43/100` | `8*43^2 > 100^2` |
| 9 | `2/5` | `9*2^2 > 5^2` |
| 10 | `2/5` | `10*2^2 > 5^2` |
| 11 | `19/50` | `11*19^2 > 50^2` |
| 12 | `3/8` | `12*3^2 > 8^2` |

These are proofs for this particular family and profile estimate, not
numerical optima. In particular the natural nine-axis family cannot
realize the favorable abstract pressure below `1/3` allowed by the
shared-minor relaxation. The specialized nine-axis profile already
has six orders equal to 360 after its initial orders 280 and 350;
the relaxed target instead allocates its next six orders as two 405s
and four 480s. The full list, not this one prefix alone, is used by the
pressure certificate.

Here is the all-scale certification argument. For a proposed rational
floor `R`, binary entropy duality says

\[
 \max_{0\le p\le1}\{R H(p)-xp/d\}
   =R\log(1+e^{-x/(dR)}).
\]

It is therefore enough to certify

\[
 \sum_i\log(1+e^{-x s_i/D})
   >mR\log(1+e^{-x/(dR)})\quad(x\ge0).                  \tag{14}
\]

On `[0,64]`, the checker uses the existing
[integer exponential and logarithm enclosures](../tools/check_multiaxis_gaussian_profile_floor.py)
at 128-bit fixed-point scale. Both sides of (14) decrease. On each
rational interval `[a,b]`, a lower enclosure for the left side at `b`
is strictly larger than an upper enclosure for the right side at `a`.
Exact dyadic bisection supplies complete covers with respectively

\[
 69,\ 32,\ 39,\ 39,\ 47,\ 48,\ 55,\ 30,\ 80,\ 44,\ 77
\]

intervals for dimensions two through twelve. Every interval and its
positive integer gap is independently replayed; no floating optimizer
or logarithm is part of the certificate.

For `x>=64`, the smallest order is `D/d`. Since
`log(1+z)>=z/2` and `log(1+z)<=z` on `0<=z<=1`, it suffices that

\[
 \exp\!\left(64\frac{1-R}{dR}\right)>2mR.                \tag{15}
\]

The checker proves (15) by a 65-term positive rational Taylor partial
sum. The ratio only increases beyond 64, so this covers the entire
tail. Combining (11)--(15) proves the stated generic-profile floors.

## 6. Replay, scope, and the remaining construction problem

```sh
python3 tools/check_cube_symmetric_exterior_blocks.py \
  --output evidence/cube-symmetric-exterior-blocks.json
```

The [checker](../tools/check_cube_symmetric_exterior_blocks.py) and
[receipt](../evidence/cube-symmetric-exterior-blocks.json) retain:

- 96 complete elementary-form identities after constructed binary
  basis changes, on one through twelve generators.
- 28 complete comparisons between unreduced and block Smith profiles
  for independently chosen symmetric polynomial coefficients.
- 64 independent binary Toeplitz-rank checks of the series orders.
- Every specialized profile and every rational all-scale certificate
  in the table, including the full infinite tail argument.

The block theorem applies to arbitrary symmetric coefficients and is
not limited to the monomial choice (9). The exclusion, however, is
limited to that natural Cayley family, the common binary physical
type, and the common-scale unweighted Smith estimate through twelve
axes. It does not exclude non-Cayley kernels, coupled physical modes,
proper physical restrictions, or better estimates of actual coefficient
ranks. It makes no all-dimension exclusion claim.

The [fixed-dimension transfer](multiaxis-cube-transfer.md) now states
exactly how a successful family below `1/sqrt(d)` would improve the
factoring score. This computation eliminates one direct candidate for
realizing the nine-axis relaxed target and supplies smaller exact
operators for further coefficient searches. No new factoring exponent
has been established.

The new block reduction, implementation, and certificates were
developed with substantial OpenAI assistant support. They are ordinary
proofs and exact finite checks, without Lean formalization or external
mathematical review.
