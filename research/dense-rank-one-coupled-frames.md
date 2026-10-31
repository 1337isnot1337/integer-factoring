# Dense rank-one coupled frames on complete physical alphabets

## Construction and scope

Every nonzero rank-one square-zero matrix over `F4` gives an explicit
coupled Laurent deformation of the four-letter reflection crossing.
The pole-one frame is polynomial for every such matrix.  Higher poles
are accepted only after the complete four-leg tensor has been added.

An exhaustive exact catalogue contains 5,355 generators and 7,176
distinct polynomial tensors.  All have the original degree 24, the
strict crossing endpoint, and no nonconstant scalar content.  This
enlarges the coordinate-transvection family in
[coupled-rational-physical-frames.md](coupled-rational-physical-frames.md).
It does not establish a pressure improvement or a new factoring
exponent.  The separate bounded search recomputes both complete
weighted profiles and uses one common pair of physical types.

## 1. The complete rank-one family

Let `G=u v^T` over a field of characteristic two.  For nonzero vectors
`u,v`, the matrix has rank one, and

```math
G^2=(v^Tu)G.
```

Thus `G^2=0` precisely when `v^Tu=0`.  Normalize the first nonzero
entry of `u` to one.  The image of `G` then determines this vector
uniquely, and `G` determines `v` uniquely.

Over `F4` in dimension four there are

```math
\frac{4^4-1}{4-1}=85
```

normalized choices of `u`.  Its nonzero annihilator vectors number
`4^3-1=63`.  Consequently this parameterization enumerates every
nonzero rank-one square-zero matrix exactly once, giving
`85*63=5355` matrices.

The generator is allowed to be dense.  It acts on the complete
physical alphabet; it need not be induced by a transformation on
individual exterior modes.

## 2. Restore the original coordinates before charging physical types

Write the tensor legs as `x,y,u,v`; vector names in this section are
written `a,b` to distinguish them from those leg labels.  The original
reflection tensor has constant coefficient

```math
T(0;x,y,u,v)=Q_{xy}(Q^{-1})_{uv}.
```

Applying `Q^-1` on `x` and `Q^T` on `v` changes this to
`[x=y][u=v]`.  In these Bell coordinates take `G=a b^T` with
`b^T a=0` and `A=I+fG`.  The four coupled maps are

```math
A_x,\qquad A_y^{-T},\qquad A_u,\qquad A_v^{-T}.
```

They preserve the Bell product and the crossing endpoint.  Because
`G^2=0` in characteristic two, `A^-1=A`.

After returning to the original physical coordinates, the four
rank-one generators, in leg order, are

```math
(Qa)(Q^{-T}b)^T,\quad ba^T,\quad ab^T,\quad
(Q^{-T}b)(Qa)^T.                               \tag{1}
```

The implementation applies (1) directly.  This retains the intended
physical letters when the subsequent pressure calculation assigns
heights and common types.  The constant coordinate changes are
included in the seed; their effect on a typed restriction is not
discarded.

For the four-letter seed the exact factors are

```math
Q=\begin{pmatrix}
0&0&0&1\\0&1&1&0\\0&1&0&0\\1&0&0&1
\end{pmatrix},\qquad
Q^{-1}=\begin{pmatrix}
1&0&0&1\\0&0&1&0\\0&1&1&0\\1&0&0&0
\end{pmatrix}.
```

## 3. An exact contraction algorithm and complete pole range

Expand the four leg actions in an auxiliary scalar `f`:

```math
T_f(w)=\sum_{r=0}^{4} f^r T_r(w),\qquad T_0=T.
                                                               \tag{2}
```

Applying `a b^T` to one leg consists of contracting against `b` and
expanding against `a`.  It uses `O(q^4)` field additions and scalar
multiplications on a tensor with `q` letters per leg.  Updating the
auxiliary powers in descending order takes ten such contractions
across all four legs.  If the starting degree is `D`, the complete
expansion therefore uses `O(q^4 D)` bit operations over fixed `F4`
with the direct coefficient representation.  A proposed pole `k`
can then be checked with `O(q^4(D+k))` bit operations.

Specifically, substitute `f=w^-k`, multiply by `w^(4k)`, add all
five tensors, and check divisibility by `w^(4k)`.  Only this final
sum is required to be polynomial.  The independent checker applies
the four Laurent matrices sequentially with exponent dictionaries;
it also permits poles in every intermediate leg action.

Every coefficient tensor `T_r` in (2) has degree at most `D`.  If
`k>D`, the intervals

```math
[-rk,D-rk],\qquad r=0,1,2,3,4
```

are disjoint.  All four nonzero-order intervals are negative.
Therefore an accepted pole greater than `D` is possible only when
`T_1=T_2=T_3=T_4=0`, meaning the frame acts trivially for every `f`.
The exhaustive catalogue checks that no listed generator has this
property.  Testing `1<=k<=D=24` is consequently a complete test of
all positive pole orders for this family.

The accepted-pole counts are:

| Pole order | Generators accepted |
|---:|---:|
| 1 | 5,355 |
| 2 | 1,755 |
| 3 | 21 |
| 4 | 21 |
| 5 | 6 |
| 6 | 6 |
| 7 | 6 |
| 8 | 6 |
| Every order at least 9 | 0 |

There are 3,600 generators admitting only order one; 1,734 admitting
exactly orders one and two; 15 admitting exactly orders one through
four; and six admitting exactly orders one through eight.

For a polynomial result, the coefficient at degree `D` comes only
from `T_0`: every introduced monomial has strictly negative degree.
The degree therefore stays exactly 24.  The local maps at `w=1` are
finite inverse-transpose pairs, so the endpoint stays the strict
crossing.  The unchanged invertible column-cut determinant at zero
precludes common positive `w`-content.  At any other finite place the
frames are local units, so the primitive source cannot acquire another
scalar factor.  No scalar content is removed in this construction.

## 4. Degree-preserving constant companions

For the same square-zero `G`,

```math
I+(f+c)G=(I+cG)(I+fG).
```

Hence every accepted `f=w^-k` also gives a polynomial tensor with
`f=w^-k+c`, for every fixed `c` in `F4`.  The added factor is a
constant invertible coupled physical map.  It preserves the exact
degree and both endpoint pairings.  When `c=1`, the four frames are
individually identity at `w=1`.

These companions can have different typed pressure profiles.  A
constant physical change of basis need not preserve a fixed diagonal
physical type subspace.  Their pressure requires a separate complete
calculation; polynomiality alone is not a rank improvement.

## 5. Exact checks and numerical boundary

The generator checker uses separate `F4` multiplication to verify all
5,355 nonzero matrices, square-zero identities, vanishing two-by-two
minors, and uniqueness.  Selected dense cases are then expanded by
an independent sequential Laurent calculation in Bell coordinates,
followed by an explicit restoration.  The controls compare 8,192
whole tensor entries and 5,376 endpoint entries, and independently
reject 11 pole choices.

The exhaustive replay reconstructs every accepted tensor and compares
its full coefficient digest, degree, and monomial content.  All 7,176
polynomial tensors are distinct.  Its input and source digests are
retained in the evidence.

```sh
python3 tools/dense_coupled_reflection_frames.py \
  --interval 4 --max-pole 24 \
  --output /tmp/dense-rank-one-catalogue.json
python3 tools/check_dense_coupled_reflection_frames.py \
  --catalogue evidence/dense-rank-one-coupled-frames.json \
  --output /tmp/dense-rank-one-check.json
```

The pressure search uses exact polynomial coefficients and complete
native local profiles, followed by floating Gibbs and gauge
optimization.  It retains every resolved best call, including when
an optimizer reaches a time limit.  A threshold candidate would still
require independent exact profile and pressure certification.

Source and evidence:

- [Dense constructor](../tools/dense_coupled_reflection_frames.py)
- [Independent checker](../tools/check_dense_coupled_reflection_frames.py)
- [Complete catalogue](../evidence/dense-rank-one-coupled-frames.json)
- [Exact replay](../evidence/dense-rank-one-coupled-frames-check.json)
- [Bounded pressure search](../tools/search_dense_coupled_reflection_frames.py)
