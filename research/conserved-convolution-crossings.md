# Conserved convolution crossings and a two-cut obstruction

## Result and scope

Let `G` be a finite abelian group of order `d>1`, over any coefficient
field. Consider tensors of the form

```math
B(x,y,u,v)=[x+y=u+v]\,f_{u+v}(x-v),\qquad x,y,u,v\in G. \tag{1}
```

Put `m=|2G|`, where `2G={2x:x in G}`. Every nonzero tensor in this
linear space has column-cut rank at least `m`. Consequently an exact
decomposition of the crossing into `k` tensors of form (1) satisfies

```math
\max_i r(B_i)c(B_i)\ge\frac{d^2m}{k}.                  \tag{2}
```

The class is closed under tensor products and arbitrary coefficient
recombination. Thus on the full label domain `G^t`, any
`exp(o(t))`-piece decomposition within this class has normalized
contraction at least

```math
\rho\ge\frac12+\frac{\log m}{4\log d}-o(1).             \tag{3}
```

For a cyclic group of even order, `m=d/2`. At `d>=64` the lower
bound is at least `17/24>1/sqrt(2)`. For the cyclic powers of two
`d=4,8,16,32`, a separate exact certificate below excludes a
half-rate improvement from the usual two-cut Smith envelope.

These statements concern this conserved-convolution class on its
**full group label domain**. They are not factoring lower bounds.
Arbitrary linear label restrictions and binary prefixes need not
preserve its blocks. For small fixed groups, (3) also leaves open a
sharper coefficient-rank analysis than the Smith envelope.

## 1. Exact crossings and their blocks

The crossing is
`E_G(x,y,u,v)=[x=v][y=u]`. It belongs to (1), with
`f_w(delta)=[delta=0]` for every `w`.

More generally, choose group-algebra polynomials `F_w(z)` with
`F_w(1)=1`, where `1` is the group-algebra identity. Then

```math
P(x,y,u,v;z)
 =[x+y=u+v]\,[X^{x-v}]F_{u+v}(z)                      \tag{4}
```

satisfies `P(1)=E_G`. Its coefficient tensors lie in (1).
Thus a polynomial curve in this class is a direct exact crossing
constructor, independently of the source vertex operators.

In the row cut `(x,y)|(v,u)`, fix `w=x+y=u+v`. The block, with
row `x` and column `v`, has entry `f_w(x-v)`. It is a convolution
matrix. In the column cut `(x,u)|(v,y)`, put

```math
j=x-u=v-y.
```

The conservation identity guarantees that different `j` values are
separate matrix blocks. The block's entry at row `x` and column `v` is

```math
                         f_{x+v-j}(x-v).              \tag{5}
```

Unlike the elementary-abelian case, (5) need not itself be a convolution
matrix.

## 2. Every nonzero kernel occupies at least `|2G|` column blocks

Suppose `f_w(delta)` is nonzero. For every `x in G` set

```math
v=x-\delta,\qquad y=w-x,\qquad u=w-x+\delta.
```

This is a nonzero tensor entry. Its column-block index is

```math
                         j=2x-w-\delta.               \tag{6}
```

As `x` varies, these indices form one coset of `2G`, of size `m`.
Thus at least `m` separate column blocks are nonzero. Each has rank
at least one, proving `c(B)>=m` over every field.

In a decomposition `sum_i B_i=E_G`, row-rank subadditivity gives
`sum_i r(B_i)>=d^2`. Some nonzero summand therefore has
`r(B_i)>=d^2/k`, and its column rank is at least `m`. This proves (2).

Tensoring (1) on `G` and `H` gives the same formula on `G×H`,
with kernel equal to the product of the two kernels. Linear
recombinations also retain (1). On `G^t` the label dimension is
`d^t` and the doubling-image size is `m^t`. Apply (2), take logarithms,
and divide by `4t log d` to obtain (3).

For `G=C_d` with even `d`, (3) becomes

```math
\rho\ge\frac34-\frac{\log2}{4\log d}-o(1).
```

At `d=64` its constant is `17/24`. The exact comparison
`2·17^2>24^2` proves `17/24>1/sqrt(2)`. The constant increases
with `d`. For odd cyclic groups `m=d` and the lower bound is `3/4`.

## 3. An exact Smith-envelope certificate for small cyclic seeds

Let `P(z)` be a nonzero degree-`D` polynomial curve from (4) with
`P(1)=E_G`. Remove a common power of `z` if necessary, and suppose
`D>=1`. Set `N=d^2`. The row Smith valuations have at least one zero;
those of the reversed column polynomial `z^D P_col(1/z)` have at
least `m` zeros by Section 2. The sum of either valuation list is at
most `ND` by its determinant degree.

For a normalized tilt `0<x<=1`, convexity therefore gives

```math
F_a(x)\ge1+(N-1)x^{N/(N-1)},\qquad
F_b(x)\ge m+(N-m)x^{N/(N-m)}.                           \tag{7}
```

Here the powers are the valuations divided by `D`. The usual Smith
bound is the infimum of `F_a(x)F_b(x)/x`, as in
[finite-seed-amplification.md](finite-seed-amplification.md),
Section 4. If there are more than the stated number of zero valuations,
Jensen applied after discarding just the stated zeros still proves (7).

Since `m>=1`, put `e=N/(N-m)` and `y=x^e`. Replacing the smaller
exponent in the first expression of (7) by `e` weakens its lower
bound. With `t=(N-m)/N`, this yields

```math
\frac{F_a(x)F_b(x)}x
\ge \frac{m+B_0y+Cy^2}{y^t},
\quad B_0=N(m+1)-2m,\quad C=(N-1)(N-m).               \tag{8}
```

For even cyclic `d`, `m=d/2`. Apply weighted AM--GM to the three
terms of (8), with weights

```math
w_0=\frac{d+2}{4d},\qquad
w_1=\frac{2d-2}{4d},\qquad
w_2=\frac14.
```

They sum to one and have weighted `y` exponent zero. Hence (8) is
at least a constant `L_d` satisfying the exact rational-power identity

```math
L_d^{4d}
=\left(\frac{2d^2}{d+2}\right)^{d+2}
 \left(\frac{d^2(d^2+2d-2)}{d-1}\right)^{2d-2}
 \left(2d(2d-1)(d^2-1)\right)^d.                       \tag{9}
```

For each `d=4,8,16,32` the checker proves, using integers only,

```math
                    (L_d^{4d})^7>d^{80d}.             \tag{10}
```

Thus `L_d>d^(20/7)`, and this Smith envelope has
`rho>5/7>1/sqrt(2)`. The last comparison is `2·5^2>7^2`.
The binary case is already excluded by the general two-label Smith
bound. Together with Section 2 for `d>=64`, this closes a half-rate
search using this Smith interface on cyclic groups of order a power
of two. Equation (10) does not assert a lower bound on the actual
coefficient ranks for the four small domains.

## 4. A concrete cyclic curve

For `d=2^a`, work in
`F_2[C_d]=F_2[X]/(X^d-1)` and put

```math
F_w(z)=z^{d-1-w}\bigl(1+(1+z)X\bigr)^w,
\qquad 0\le w<d.                                     \tag{11}
```

This is polynomial of degree `d-1` and has `F_w(1)=1`.
It therefore gives the exact crossing (4). At `d=8` its row
valuation multiplicities are

```text
1,2,3,4,5,6,7,8,7,6,5,4,3,2,1
```

and its reversed-column multiplicities are four in every degree
zero through fourteen, with eight in degree seven. The finite
floating-point optimum of that profile is about `0.885706`, well
above the target; it is recorded only as a diagnostic. The uniform
obstructions are (2)--(3) and the integer certificate (10).

## Verification

The checker constructs the full matrices from group kernels, checks
the column blocks independently, verifies exact crossing sums and
tensor-product closure, and tests (10) by integer arithmetic.
It also computes the Smith multiplicities in Section 4 using exact
truncated polynomial-module ranks.

```sh
python3 tools/check_conserved_convolution.py \
  --output evidence/conserved-convolution-crossings.json
```
