# Local alphabet bounds for observable word-rank transfers

The [linear detector lemma](linear-rank-detectors-for-observable-tables.md)
gives two useful restrictions on the observable transition interface.
First, every scalar entry space has a local-alphabet lower bound.
Second, entries in the span of the three standard pairings cannot
give a score-improving spectral transfer in characteristic two,
for any external alphabet size or auxiliary dimension.

These are bounds on the individual-word rank spectral method.
They do not exclude additional cancellations or direct rank
calculations on grouped pieces.

## 1. A bound for every scalar entry space

Let the external alphabet have size `d`, and put `n=d^2`. For a
scalar four-leg tensor `t`, let its two cut ranks be `r,c`. When
`t` is nonzero, `1<=r,c<=n`, and

```math
(n+1)^2rc-n(r+c)^2=(nr-c)(nc-r)\ge0.
```

Consequently

```math
\sigma(t)=\sqrt{rc}\ge\frac{d}{d^2+1}(r+c).     \tag{1}
```

The same inequality holds at zero. The field-linear map

```math
L(t)=\operatorname{diag}(t_{\rm row},t_{\rm col})
```

has rank `r+c`; on the crossing it has rank `2d^2`. The detector
lemma therefore proves, over every coefficient field,

```math
\boxed{\operatorname{spr}(G_{\rm live})
 \ge\frac{2d^3}{d^2+1}.}                         \tag{2}
```

This applies to any finite observable matrix table satisfying
the all-word crossing endpoint, regardless of its scalar entry
space. For `d=2,3,4`, the bound is respectively
`16/5`, `27/5`, and `128/17`. Each exceeds `d^(sqrt(2))`.
Thus a score-improving observable word-rank transfer must use a
local alphabet of size at least five. A binary construction with
base below four could still be a useful diagnostic improvement,
but its base cannot fall below `16/5` by this method.

For an exact comparison, use `sqrt(2)<3/2` for `d=2,3` and
`sqrt(2)<10/7` for `d=4`. The latter comparison reduces to
`2^29>17^7`.

## 2. The coefficient of the crossing in the three-pairing space

Now work in characteristic two, with arbitrary `d>=2`, and define

```math
E=[x=v][y=u],\quad R=[x=y][u=v],\quad
C=[x=u][y=v].
```

These tensors are linearly independent. Let
`t=aE+bR+cC`, with `a` nonzero, and scale to `a=1`.
Put `q=d(d-1)/2`, `n=d^2`, and `delta=d mod 2`.

In the row cut, the matrix is

```math
K+brr^T+cI,
```

where `K` swaps the two `d`-dimensional tensor factors and `r`
is the vectorization of the identity matrix. Thus `Kr=r` and
`r^Tr=delta`.

If `c=1`, `I+K` has rank `q`: each unordered off-diagonal pair
contributes a rank-one block. The vector `r` lies in the diagonal
coordinates, outside that image. Hence the row rank is
`q+[b!=0]`.

If `c!=1`, `K+cI` is invertible. Its rank-one update loses one
rank precisely when `delta=1` and `b=1+c`; otherwise it stays
invertible. Indeed

```math
r^T(K+cI)^{-1}r=\frac{\delta}{1+c}.
```

The column formula exchanges `b,c`. These formulas leave the
following possibilities for the minimum geometric mean with
`a!=0`:

```math
\min\left(q+1,\sqrt{q(n-\delta)},n-\delta\right)=q+1.
                                                               \tag{3}
```

For `d=2`, this is immediate. For `d>=3`, even the stronger
comparison with `n-1` gives

```math
q(n-1)-(q+1)^2=q(n-q-3)-1>0,
```

and `n-1>=q+1`. Equality in (3) is attained by `E+R+C`, which
has both cut ranks `q+1`. Therefore

```math
\boxed{a\ne0\quad\Longrightarrow\quad
 \sigma(aE+bR+cC)\ge q+1.}                       \tag{4}
```

The scalar linear map `L(aE+bR+cC)=a` has crossing rank one and
satisfies the detector inequality with multiplier `q+1`.
Every observable table in this scalar space consequently has

```math
\operatorname{spr}(G_{\rm live})\ge q+1.         \tag{5}
```

## 3. All local alphabets are covered

For `d>=6`, `(d-1)^2>=4d`, so
`q>=d^(3/2)>d^(sqrt(2))`. For `d=5`, the bound is `q+1=11`;
the exact inequalities `sqrt(2)<10/7` and `11^7>5^10` prove that
it exceeds `5^(sqrt(2))`. Section 1 covers `d=2,3,4`.

**Corollary.** No finite observable transition table with scalar
entries in `span(E,R,C)` gives a score-improving word-rank
spectral transfer in characteristic two, for any local alphabet
`d>=2` and any auxiliary matrix size.

This includes noncommuting coefficient matrices and strict
observable quotients with no common crossing eigenvector.
Increasing only the auxiliary size cannot escape the corollary.
Constructive searches must enlarge the scalar entry space or
use an estimate stronger than the individual-word spectral bound.

## 4. Exact controls

`check_local_alphabet_spectral_barriers.py` verifies the rational
bound (1) on all rank pairs through `d=12`, the direct flattening
formulas over F2 and F4 for several even and odd alphabets, the
attaining tensor `E+R+C`, and the integer comparisons used above.
The universal claims are proved algebraically in this note.
