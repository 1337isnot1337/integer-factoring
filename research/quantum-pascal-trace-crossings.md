# Quantum-Pascal trace crossings and a uniform profile obstruction

## Construction and scope

A polynomial divided-power ladder gives explicit crossing tensors over
finite fields, outside the copied-coordinate register architecture.
The default vacuum uses Gaussian binomial coefficients at a primitive
root of unity. Its complete degree is `d(d-1)` for odd alphabet size `d`.

There is also a uniform obstruction. For every polynomial diagonal
vacuum in this trace family with first entry one and a valid cold
endpoint, the bare two-place Smith criterion has `rho>=3/4`. The
unweighted multipoint Hermite criterion also cannot give `rho<3/4`,
even after extending the coefficient field to include further places.
These statements concern the displayed sufficient profile criteria.
Physical word types, physical valuation gauges, and collective rank
analysis remain outside the obstruction. No factoring exponent above
the existing `1/2` is established.

The multiplicative ladder is related to the usual quantized Weyl
relation; [Levitt and Yakimov, arXiv:1606.02121v2](https://arxiv.org/abs/1606.02121v2)
provide primary algebraic context. The construction and proofs below
use explicit finite polynomial matrices.

## 1. Polynomial divided powers and the crossing endpoint

Work over a field `F`. Let `G_(k,y)(t)=[k choose y]_t` be the Gaussian
binomial polynomial. One definition uses

```math
G_{0,0}=1,\qquad
G_{k,y}=G_{k-1,y-1}+t^yG_{k-1,y},                          \tag{1}
```

with out-of-range entries zero. In particular `G_(y,y)=1` and
`deg G_(k,y)=y(k-y)`.

Choose polynomial vacuum entries `P_0(t)=1,P_1(t),...,P_(d-1)(t)`.
Let `X e_k=e_(k+1)` for `k<d-1`, and `X e_(d-1)=0`. Define the
divided lowering operator by

```math
D_t^{[y]} e_k=G_{k,y}(t)e_{k-y}\quad(k\ge y),
```

and zero otherwise. These are polynomial operators. With
`P(t)=diag(P_0,...,P_(d-1))`, put

```math
O_{xy}(t)=X^xP(t)D_t^{[y]},\qquad
K(t;x,y,u,v)=\operatorname{Tr}(O_{xy}(t)O_{uv}(t)).          \tag{2}
```

If `P_i(a)=0` for all `i>0`, then `O_xy(a)=e_x e_y^T`.
Consequently

```math
K(a;x,y,u,v)=[x=v][y=u].                                  \tag{3}
```

For `a!=0`, substitution `t=az` puts the endpoint at one without
changing the polynomial degree or local Smith orders at corresponding
places. Thus (2) gives the complete physical endpoint directly.

For the default family take `P_i(t)=G_(d,i)(t)`, `0<=i<d`, and let
`a` be a primitive `d`th root in a fixed finite extension, with
`char(F)` not dividing `d`. The polynomial identity

```math
G_{d,i}(t)[i]_t![d-i]_t!=[d]_t!
```

proves `P_i(a)=0` for `0<i<d`: the factors below `d` are nonzero
at `a`, while `[d]_a=0`. In characteristic two every fixed odd `d`
is available in a finite extension. Divided powers are used throughout;
normalization by a parameter-dependent quantum factorial is already
performed by the polynomial `G_(k,y)`, with no uncharged rational pole.

## 2. Actual degree and exact grade blocks

The ladder entry before the final shift is

```math
A_{k,y}=G_{k,y}P_{k-y}.
```

Nonzero tensor entries satisfy the integer charge relation
`x+u=y+v`. The constructor writes them as the finite polynomial sum

```math
K(x,y,u,v)=\sum_{\substack{v\le k<d\\y\le k-v+u<d}}
                  A_{k,v}A_{k-v+u,y}.                     \tag{4}
```

Both cuts split into exact square blocks of order at most `d`.
The first cut is indexed by `x-y=v-u`; the second by `x+u=y+v`.
These are literal support blocks, so the decomposition loses no rank.

For the default vacuum, writing `i=k-y` gives
`deg A_(k,y)=i(d+y-i)`. Since `i+y<=d-1`, this is at most
`floor(d(d-1)/2)` when `d` is odd. Hence the tensor degree is at most
`d(d-1)`. Equality holds at `x=y=u=v=(d-1)/2`, whose unique largest
path has `i=(d-1)/2` and monic leading coefficient. Thus

```math
\deg K=d(d-1)\qquad(d\text{ odd}).                         \tag{5}
```

Every seed has fixed degree, field, and alphabet. Tensor powers,
coefficient tables, and both exact factorizations therefore retain
the existing `2^O(lambda)` construction cost at binary label width
`lambda`. The [polynomial transfer](polynomial-crossing-interpolation.md)
and [multipoint transfer](multipoint-hermite-crossings.md) apply when
their strict numerical criteria hold.

## 3. A unit determinant in the first cut

Vectorize the `d^2` operators `O_xy` in the matrix-unit basis. Group
them by the difference `x-y`. In each group, the coefficient at input
column `k` is zero for `k<y`; its diagonal value at `k=y` is
`G_(y,y)P_0=1`. The change-of-basis matrix is therefore block triangular
with unit diagonal. Its determinant is one.

The bilinear form `(A,B) -> Tr(AB)` is a permutation pairing on
matrix units. The first cut of (2) is its Gram matrix in the above
polynomial basis. Hence

```math
\det K_{xy\mid uv}=\pm1.                                 \tag{6}
```

This holds for every field and all the stated polynomial vacua,
independently of their degrees. In particular the finite-place first-cut
Smith orders are all zero. At infinity, the sum of the raw orders is
zero, so its Laurent partition sum satisfies `U(s)>=d^2` for `s>=1`
by arithmetic--geometric mean.

## 4. Constant entries in the second cut

For every `0<=y<d`, the second-cut grade `d-1+y` has the literal entry

```math
K(d-1,y,y,d-1)=1.                                        \tag{7}
```

Only the top input path contributes: both divided-power diagonal
coefficients are one and both vacuum positions are zero. These `d`
entries lie in distinct grade blocks. At every finite place, therefore,
at least `d` second-cut Smith orders are zero.

Combining (6)--(7) proves the all-scale bare bound

```math
B(s)=U_{\rm row,\infty}(s)V_{\rm col,0}(s)\ge d^3,
\qquad \rho\ge\frac34.                                  \tag{8}
```

The opposite orientation has the finite first-cut sum `d^2`, while
the other raw infinity sum is at least `d^2` because its determinant
is a nonzero polynomial. It gives no smaller parameter.

## 5. The unweighted multipoint criterion has the same obstruction

At any finite interpolation point, (6)--(7) give local partition sums
at least `d^2` and `d`. A nonnegative local degree allocation only
increases the corresponding pressure. Thus no finite place can satisfy
the multipoint strict criterion below `rho=3/4`.

If infinity is the only selected place, it must receive the complete
degree budget `D`. Its normalized Smith orders have mean at most `D`
in both cuts: their sums are `d^2 D` minus the nonnegative degrees of
the respective determinants. Arithmetic--geometric mean then makes
the pressure at allocation `D` at least `4 log d`, giving `rho>=1`.
This proves the asserted `3/4` floor for the unweighted multipoint
criterion. Field extensions introduce more finite places, but (6)--(7)
hold at every one of them, so the conclusion is unchanged.

The argument does not apply the uniform determinant mean after a
physical word-type restriction. Such a restriction changes the
operator space, entropy, and available grade blocks. Its complete
cost must be analyzed separately using the
[type-and-gauge transfer](type-optimized-crossing-pressure.md).

## 6. Exact controls

```sh
make verify-quantum-pascal
```

The standard-library checker verifies complete polynomial endpoints at
3, 5, 7, 9, and 15 letters, quantum-factorial identities, the exact
degree, every displayed constant-one block entry, and independent
literal operator traces. At three and five letters it computes the
unit block determinants directly and recovers both bare spectra with
fraction-free polynomial minors. Rational scale comparisons replay (8).

The [constructor](../tools/quantum_pascal_crossing.py),
[checker](../tools/check_quantum_pascal_crossings.py), and
[record](../evidence/quantum-pascal-crossings.json) distinguish the
complete finite identities from the uniform proofs in Sections 3--5.
This construction and analysis used substantial OpenAI Codex assistance.
No new Lean work or stronger factoring score is claimed.
