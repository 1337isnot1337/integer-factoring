# Bounded-weight source duals and primitive crossing directions

## Result and scope

The exact neutral compiler accepts every regular operator monomial of
weight at most the physical state weight. This adds lower-weight dual
parameters to the earlier [homogeneous completion](homogeneous-dual-completion.md).
At physical weight 16, the full parameter space has dimension 33 and
the physical alphabet has five states. The 140 formal dual parameters
produce an **exactly 50-dimensional affine tensor family**.

Every tensor in the family has the literal crossing endpoint over every
characteristic-two coefficient field. Its coefficients, full source
normalization, and relations are retained exactly. The family includes
both the homogeneous duals and the original integral inhomogeneous duals.

The resulting source seeds do not currently have competitive pressure.
Their useful additional output is a library of endpoint-zero directions.
Dividing each direction by a common polynomial factor and restoring one
factor `1+w` gives primitive polynomial directions that can be added to
any strict crossing on the same physical alphabet. This last operation
preserves the endpoint; it does not retain a source-operator interpretation.

## 1. The general source angle factor

Let the physical states `y,v,u` have weight `h`, and let an operator
monomial `a` have weight `w<=h`. Its source angle operation is

\[
 a^{\langle i\rangle}=\binom{h+w}{i}a.
\]

In the complete piece-generating series, put `r=s-h+i`. Then

\[
 \begin{aligned}
 \sum_s z^s Q_s(a,b)
 &=\sum_i(-1)^i\binom{h+w}{i}z^{h-i}
       \sum_r z^r a[r]b[-r]\\
 &=z^{-w}(z-1)^{h+w}\sum_r z^r a[r]b[-r].                \tag{1}
 \end{aligned}
\]

This is an identity of the source coefficient series; the source support
bound supplies the finite polynomial after cancellation. In characteristic
two, `z-1=1+z`.

Choose any power of two `Q` exceeding every part used. In the
[neutral Wick kernel](neutral-reflected-crossings.md), replace each part
`m` of `a` by an `A` node indexed `Q-m`. The product of the `A` gauges
is `z^-w(1+z)^w`; the `B` gauges contribute `(1+z)^h`. Their product is
exactly (1). The physical `Y,U` parts remain unchanged. The same Wick
proof therefore evaluates the complete source entry without assuming
`w=h` or equal particle counts for the parameter and states.

Odd total field count gives zero automatically. This is a real source
parity cancellation, not a deletion of a nonzero parameter.

## 2. All bounded right inverses

Let

\[
 V_{\le h}=\bigoplus_{w=0}^h W_{2,w},\qquad
 G(a)_j=[1]a[-h]b_j.                                    \tag{2}
\]

Suppose the homogeneous weight-`h` pairing is nonsingular and choose its
dual basis `a_0,...,a_(d-1)`. For every regular monomial `m` of weight
strictly less than `h`, put

\[
 k_m=m-\sum_jG(m)_j a_j.                                 \tag{3}
\]

These vectors form a basis of `ker G`: their lower-weight coordinates
are the distinct monomial coordinates, while the weight-`h` pairing is
invertible. Consequently every bounded dual is uniquely

\[
 \widetilde a_x=a_x+\sum_{|m|<h}c_{x,m}k_m.               \tag{4}
\]

Linearity of the source crossing in its first parameter makes (4) an
affine polynomial tensor family. Its endpoint is unchanged by (2)--(3).
The proof works over every characteristic-two extension field.

At `h=16`, `dim W_(2,16)=5` and `dim V_(<=16)=33`. There are 28 kernel
vectors per dual, hence `5*28=140` formal parameters. At `h=8`, the
corresponding dimensions are two and nine, giving fourteen formal
parameters. These counts are dimensions of parameter spaces, not cut ranks.

## 3. The complete tensor image

For a fixed first logical label, a kernel vector gives a slab indexed
by `(y,u,v)`. Pack every polynomial coefficient into a binary vector
and perform exact elimination. The complete slab ranks are:

| Weight | Physical alphabet | Formal kernel vectors per dual | Nonzero independent slabs | Complete tensor directions |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 2 | 7 | 2 | 4 |
| 16 | 5 | 28 | 10 | 50 |

For the retained weight-16 ordering, the ten nonzero slab indices are

\[
 0,3,5,7,8,11,15,16,21,23.                               \tag{5}
\]

They are independent; every other slab is identically zero. Inserting
each at any of the five first-label values gives disjoint coefficient
supports, so the total affine tensor dimension is fifty. Binary
independence remains independence over every characteristic-two field.

The complete weight-16 family has common lowest Laurent exponent `-32`
and verified exponent spacing two. In the compressed polynomial variable
its degree is at most 32. The larger spacing of the homogeneous reference
alone must not be applied to the new directions.

The original reduced integral duals expand in the same space. The checker
recovers their complete tensors from (4), using their actual lower-weight
coefficients, rather than merely checking their endpoints.

## 4. Primitive polynomial directions

Let `Delta(w)` be any nonzero polynomial tensor with `Delta(1)=0`, and
let `g(w)` be the monic greatest common divisor of all its nonzero scalar
entries. Then

\[
 \Delta_{\rm prim}(w)=(1+w)\,\Delta(w)/g(w)               \tag{6}
\]

is a polynomial tensor with zero endpoint. If `K(1)` is a strict crossing,
then `K+c Delta_prim` has the same endpoint for every field scalar `c`.
Multiplication of (6) by any additional polynomial also preserves this
endpoint. The degree is at most `deg Delta-deg g+1`.

This statement concerns actual coefficient divisibility. It does not
treat the common divisor as a free operation on cut profiles: all
profiles of the new tensor must be recomputed. It permits endpoint-zero
directions from this source family to be tested around a stronger,
unrelated reference crossing.

The scalar parameter `a=1` has a particularly transparent interpretation.
Its only mode is `a[0]=I`, and (1) gives

\[
 \Delta(z)=(1+z)^h B,
 \qquad B(y,u,v)=[u]b_v[0]y.                             \tag{7}
\]

After the weight-16 family's common translation and `w=z^2` compression,
this is `w^16(1+w^8)B`. The slab has twenty nonzero entries, all with
that same polynomial. Formula (6) is therefore exactly `(1+w)B`.
This degree-one direction follows directly from the identity mode.

The separate saturated-direction implementation can similarly reduce
six other independent slabs to degree 17. The literal bounded-source
data retain their original factors, so the source family and these
new endpoint-preserving directions can be audited separately.

## 5. Reproduction and finite scope

```sh
python3 tools/check_bounded_dual_completion.py
```

The [constructor](../tools/bounded_dual_completion.py) and
[checker](../tools/check_bounded_dual_completion.py) reconstruct the
complete coefficient data at weights eight and sixteen. The checker
compares 2,376 coefficients with the independent original two-operator
sum, verifies 154 kernel vacuum entries, checks 5,128 F4-specialized
endpoint entries, and proves the tensor-image ranks and all relations.
It also checks independence of the auxiliary Wick period.

The complete [weight-16 coefficient data](../evidence/bounded-dual-completion-h16.json)
and [exact verification record](../evidence/bounded-dual-completion-checks.json)
are retained. The original integral reference and the homogeneous
reference are both included. None of these checks asserts a favorable
Smith pressure or a stronger factoring exponent.

This construction and implementation used substantial OpenAI assistant support.
