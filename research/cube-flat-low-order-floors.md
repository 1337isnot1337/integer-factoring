# Two low-order obstructions for flat cube tensors

## Result

Consider an eight-leg binary tensor of common homogeneous polynomial
degree `D`, with the antipodal equality endpoint at `(1,1,1)`. Odd
physical coefficients are allowed, and the coefficient field is
arbitrary. Let `s_(a,1)<=...<=s_(a,16)` be the generic Smith orders at
coordinate cut `a`. Its common-degree pressure is

```math
R(t)=\frac{\frac13\sum_{a=1}^3\log\left(\sum_{i=1}^{16}t^{s_{a,i}}\right)
 -(D/3)\log t}{4\log2},\qquad 0<t<1.
\tag{1}
```

Remove the common coordinate monomial by setting

```math
m_a=s_{a,1},\qquad D_*=D-\sum_a m_a,\qquad
u_{a,i}=s_{a,i}-m_a.
\tag{2}
```

If `D_*=0`, the pressure is one. If `D_*>0`, either of the following
conditions implies an all-scale obstruction:

| Condition on the normalized generic orders | Consequence |
|---|---|
| `u_(a,2)=0` for at least one axis | `R(t)>29/50` for every `0<t<1` |
| `sum_a u_(a,2)<=D_*` | `R(t)>29/50` for every `0<t<1` |

Since `29/50>1/sqrt(3)`, any flat tensor that crosses the cube target
must have **rank one on every normalized coordinate face**, and must
also satisfy

```math
\frac13\sum_a u_{a,2}>D_*/3.
\tag{3}
```

These restrictions use only the determinant degree. They apply to
parity-breaking tensors as well as even tensors. The separate
[parity floor](cube-flat-parity-floor.md) excludes the whole even
family. Neither statement is a lower bound for integer factoring.

The [standard-library checker](../tools/check_cube_flat_low_order_floors.py)
verifies two small rational entropy certificates by exact integer
arithmetic. Its [evidence](../evidence/cube-flat-low-order-floors.json)
records the weights and comparisons.

## 1. Normalization and determinant degree

For each axis, the first Smith order is the minimum coordinate
valuation among all tensor coefficients. Thus every coefficient is
divisible by `product_a x_a^m_a`. Dividing this monomial preserves the
endpoint, replaces the degree by `D_*`, and replaces the orders by
`u_(a,i)`. Substitution in (1) shows that the pressure is unchanged.

All normalized orders are nonnegative and `u_(a,1)=0`.
If `D_*=0`, all normalized coefficients are constant, so the tensor is
its endpoint. Each cut is an invertible constant matrix, and (1) is one.

Suppose `D_*>0`. Each full sixteen-by-sixteen cut determinant is a
nonzero homogeneous polynomial of degree `16D_*`, because its endpoint
is an invertible permutation matrix. Therefore

```math
\sum_{i=1}^{16}u_{a,i}\le16D_*\qquad\text{for every axis }a.
\tag{4}
```

Only (4), the zero first orders, and the extra hypothesis of each case
are needed below. The same argument works for coefficients of degree
at most `D` when that degree cap is charged in (1).

## 2. Any coordinate face of rank at least two is too costly

Suppose `u_(B,2)=0` at one axis `B`, and call the other axes `A_1,A_2`.
Use these probability vectors on the sorted orders:

| Axis class | Weights |
|---|---|
| Either `A` axis | `(50,2,2,...,2)/80`, with fifteen copies of two |
| Axis `B` | `(33,33,1,...,1)/80`, with fourteen copies of one |

For an `A` axis the first order is zero. Equation (4) bounds its
weighted mean order by `(2/80)*16D_*`.
For `B`, the first two orders are zero, so the bound is
`(1/80)*16D_*`. Averaging the three bounds gives exactly

```math
\frac13\left(2\frac2{80}+\frac1{80}\right)16D_*=D_*/3.
\tag{5}
```

For any probability vector `w`, the log-sum inequality is

```math
\log\sum_i t^{u_i}\ge H(w)+\log t\sum_iw_i u_i.
\tag{6}
```

Because `log t<0`, (5) and (6) bound the pressure below by
`(2H_A+H_B)/(12 log 2)`. The exact comparison

```math
80^{1200}
 >2^{2784}(50^{50}2^{30})^{10}(33^{66})^5
\tag{7}
```

proves that this constant exceeds `29/50`. Its approximate value is
`0.5832782415`; the proof uses the integer comparison alone.

## 3. A small average second order is also too costly

Now suppose `sum_a u_(a,2)<=D_*`, without a face-rank assumption.
At every axis use the weights

```math
w=(33,12,1,\ldots,1)/59,
\tag{8}
```

where one appears fourteen times. The first order is zero. By (4),
the average weighted mean order is at most

```math
\frac{16D_*}{59}
 +\frac{11}{59}\left(\frac13\sum_a u_{a,2}\right)
 \le\frac{16D_*}{59}+\frac{11D_*}{177}=D_*/3.
\tag{9}
```

Equation (6) therefore gives `R(t)>=H(w)/(4 log 2)`.
The second exact comparison is

```math
59^{1475}>2^{3422}(33^{33}12^{12})^{25}.
\tag{10}
```

It proves `H(w)/(4 log 2)>29/50`; the value is approximately
`0.5830143388`. Finally `3*29^2=2523>2500=50^2`, which certifies the
comparison with `1/sqrt(3)`.

## 4. Use in constructive searches

The tests apply to the actual generic profiles. A minimum that appears
only after fixing the other two coordinates need not define a common
multivariate monomial, so specialized minima cannot be substituted into
the normalization step without an additional certificate.

For common degree six with no common coordinate factor, a candidate
needs a unique zero order at each cut and an average second order
strictly greater than two. The reusable
[flat seeds](../tools/cube_flat_seeds.py) provide an even reference and a
genuinely odd reference. Both have generic second order two at every
cut, so the second obstruction already excludes them as target-crossing
tensors. They remain useful for constructing exact local frames and
searching larger coefficient families.

Run the exact certificate with:

```sh
python3 tools/check_cube_flat_low_order_floors.py
```

The checker also exposes `low_order_obstructions(profiles, degree)` for
already certified generic profiles. It reports the applicable conditions;
it neither computes nor certifies a Smith profile.
