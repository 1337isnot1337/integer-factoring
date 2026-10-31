# A uniform all-scale comparison for an eight-parameter reflection core

## Result and scope

The bounded-core reduction makes it possible to exhaust an explicit
eight-parameter family over `F4` on a long reflection-compressed interval.
At each of six directly checked interval sizes, all 65,536 coefficient
rows give only 234 distinct pairs of complete active Smith profiles.

Every such pair is dominated by the Gaussian in a stronger, exact sense:
its sorted pairwise sums of zero and infinity orders are coordinatewise
at most the Gaussian sums. Thus the two-cut partition product cannot
improve at **any real scale `s>1`**, at the natural gauges. A bivariate
polynomial degree argument extends the calculation at `M=128` to every
larger power-of-two interval. Together with the finite cases, this proves
the assertion for **every power of two `M>=4`**. No numerical pressure
comparison is needed for that conclusion.

This is a uniform result for the stated coefficient family. It does not
exclude other coefficient polynomials, physical state gauges, larger
cores, or other crossing families. The factoring score remains `alpha=1/2`.

## 1. The completely enumerated family

Use the [reflection hot-core construction](reflection-hot-core-defects.md)
with one selected logical mode on each leg, in physical variable order
`x,y,u,v`. Put `P=1`. Starting from the Gaussian, multiply by

```math
1+\sum_{H\in\{xy,xu,yu,xv,yv,uv\}}g_H f_H(w)X_H
 +\bigl(a(w^{-4}+1)+b(1+w^4)\bigr)xyuv,
```

where all eight coefficients belong independently to `F4`, and

```math
f_{xy}(w)=w^{-4}+w^{-2},\qquad
f_H(w)=1+w^2\quad(H\ne xy).
```

Every nonconstant exterior coefficient vanishes at `w=1`. Common
monomial clearing therefore gives a polynomial crossing with the
original exact endpoint. Its scalar shifts cancel between the two cuts.
The reduction retains three active modes, with orders
`2,2M+2,2M+6`; every remaining Gaussian mode is included as a spectator.

The complete recorded intervals are `M=4,8,16,32,64,128`. The raw
records are [M4](../evidence/reflection-hot-grid-M4.json),
[M8](../evidence/reflection-hot-grid-M8.json),
[M16](../evidence/reflection-hot-grid-M16.json),
[M32](../evidence/reflection-hot-grid-M32.json),
[M64](../evidence/reflection-hot-grid-M64.json), and
[M128](../evidence/reflection-hot-grid-M128.json).

The same conclusion holds after every common positive integer
substitution `w -> w^P`, with all natural heights multiplied by `P`.
Every Smith order and physical energy is multiplied by `P`, so this
only replaces `s` by `s^P` in the comparison.

## 2. A finite exact order certificate implies every scale

Let `a_1,...,a_8` and `b_1,...,b_8` be the active orders at zero and
infinity, before the cancelling scalar shifts. Let
`c_1<=...<=c_64` be the sorted multiset of their pairwise sums. Then

```math
Z_0(s)Z_\infty(s)=\sum_{k=1}^{64}s^{-c_k}.
```

Write `c_k^G` for the corresponding Gaussian list. The exact integer
comparison in every recorded class is

```math
c_k\le c_k^G\qquad(1\le k\le64).
```

For `s>1`, each summand on the candidate side is therefore at least the
Gaussian summand. The active partition product is at least the Gaussian
product. Multiplying by the identical spectator products preserves the
inequality. Both crossings use the same physical alphabets and height
gauges, so their physical type partitions are identical as well.

This proves the claimed no-improvement statement for all scales,
including scales outside the original search grid. Seven profile classes
can improve a zero-cut order, but their infinity costs prevent a gain in
the combined comparison.

## 3. Reproduction and independent checks

The [search driver](../tools/search_reflection_hot_core_filters.py)
enumerates `F4^8` literally. Each source component is a complete Laurent
polynomial matrix. The exact native Smith window expands as needed;
unresolved cases prevent the record from being marked complete.
The record stores all distinct profile pairs, their multiplicities,
representative coefficients, and source digests.

The [comparison checker](../tools/check_reflection_hot_core_grid.py)
rebuilds every representative with the separate literal polynomial
constructor, compares its two profiles, and verifies all 64 exact order
inequalities. Its optional exhaustive pass recomputes every coefficient
row and checks the entire profile-pair histogram. The
[verification record](../evidence/reflection-hot-grid-comparisons.json)
states precisely which replays were performed.

```sh
python3 tools/check_reflection_hot_core_grid.py --exhaustive \
  --output /tmp/reflection-hot-grid-comparisons.json
```

Floating optimized roots in the raw search records are exploration data.
The theorem in Section 2 uses only integer profiles and their exact
coordinate comparisons.

## 4. One finite interval determines every larger interval

At either local place put `c=1+z^2` and introduce an independent
variable `T`. It will later be specialized to `T=z^(2M)`. The three
active Gaussian one-particle entries are

```math
z^2,\quad z^2T,\quad z^6T.
```

The only interval-dependent geometric sum in the transported physical
operators is

```math
g_M=\sum_{j=0}^{M/2-2}z^{4j}
   =\frac{z^4+T}{z^4c^2}.
```

In the notation of the hot-core proof, the right-first operator is

```math
r_1=z^{-2}A_0+cT^{-1}z^{-2}A_1
       +\frac{z^4+T}{Tz^6c}A_2.
```

Every physical exterior monomial uses this operator at most once.
Multiply the entire active matrix by `c`. A component containing `r_1`
can use the denominator-free operator

```math
cr_1=cz^{-2}A_0+c^2T^{-1}z^{-2}A_1
        +(z^{-2}T^{-1}+z^{-6})A_2.
```

Every other component receives the scalar `c` directly. The other
three physical operators have only monomial coefficients. Thus each
entry is an explicit Laurent polynomial in `z,T`, independently of `M`.
After applying the specified physical heights and scalar coefficient
polynomials, direct symbolic expansion gives the uniform entry bounds

```math
0\le\deg_z\le20,\qquad 0\le\deg_T\le2.
```

These are bounds on the union of the supports of all nine source
matrices at both places, so they hold for every coefficient row. No
common monomial clearing is required. The factor `c` is a unit in each
local ring and changes no local Smith order. This is a local matrix
representation; it does not multiply the original global crossing by
an uncharged polynomial.

Every `k`-minor, for `1<=k<=8`, consequently lies in
`F4[z,T]` with `z`-degree at most `20k<=160`. If `T=z^L` with
`L>160`, different `T`-degree bands have disjoint `z`-degree ranges.
Specialization therefore preserves nonzero minors and orders their
leading terms lexicographically, first by `T` exponent and then by `z`
exponent. The same global range applies to **all** `k`-minors, which is
needed to compare their minimum orders with each other.

At `M=128`, one has `L=256`. If `d_k` is the sum of the first `k`
Smith orders, the leading bidegree of the minimum `k`-minor is exactly

```math
(t_k,z_k)=\bigl(\lfloor d_k/256\rfloor,\ d_k\bmod256\bigr).
```

For every larger interval, that determinantal minimum is therefore
`2M*t_k+z_k`. Subtracting consecutive minima gives each individual
Smith order as an explicit affine function of `2M`. Equal profiles at
`M=128` cannot split into different affine profiles later: the bounded
residue range makes every decoded bidegree unique.

The [bivariate checker](../tools/reflection_hot_core_uniform.py) constructs
these matrices independently of the geometric-sum implementation and
compares 4,096 complete component entries at four finite intervals.
It then decodes all 234 profile classes. In each class, sort the 64
pairwise affine orders lexicographically by slope and intercept, and
pair them with the similarly sorted Gaussian orders. Each difference
`Gaussian minus candidate` has nonnegative slope and nonnegative value
at `M=128`. The difference stays nonnegative for every larger `M`.
This supplies a fixed term-by-term comparison of the partition products
for every `s>1`, as in Section 2.

The complete affine profiles and exact margins are stored in
[the uniform certificate](../evidence/reflection-hot-core-uniform.json).
The finite comparisons at `M=4,8,16,32,64` cover the remaining
power-of-two intervals. This proves the all-interval statement in the
opening section.
