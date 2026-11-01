# Physical types in same-place Hermite reconstruction

## Result and scope

Physical letter weights extend the [same-place Hermite criterion](multipoint-hermite-crossings.md),
but they do not remove its obstruction at `rho = 3/4`. This note gives the
weighted criterion, including independent paired types and independent cut
scales, and proves that every strict certificate from that criterion has
`rho > 3/4`. The theorem permits arbitrary fixed alphabets and arbitrary
integer leg heights at every selected place.

This is a restriction on a sufficient rank estimate. It is not a lower
bound for actual coefficient ranks, arbitrary interpolation functionals,
or factoring algorithms. In particular, a piece with its two rank
estimates taken from different sets of places is outside the theorem.
The established factoring exponent remains `1/2`.

## 1. The weighted criterion

Let `K(z;x,y,u,v)` be a polynomial tensor of degree at most `D` over a
fixed finite field. The alphabets on `x,v` are `X`, and those on `y,u`
are `Y`. Assume the literal crossing endpoint

```math
K(1;x,y,u,v)=[x=v][y=u].                                      \tag{1}
```

Choose distinct places `a_j` on the projective line over the field,
all different from one. At a finite place put `K_j(t)=K(a_j+t)`.
At infinity put `K_j(t)=t^D K(t^{-1})`. These are the same homogeneous
sections used in the unweighted Hermite construction.

Fix rational types `pi_X,pi_Y`, with entropy sum

```math
H=H(pi_X)+H(pi_Y)>0.                                         \tag{2}
```

At each place and each cost cut independently, choose integer heights
`h_(j,c,l)(a)` for every physical letter. For `c=0,1`, multiply the
local tensor entry by the corresponding power

```math
t^{h_{j,c,x}(x)+h_{j,c,y}(y)+h_{j,c,u}(u)+h_{j,c,v}(v)}.
                                                                  \tag{3}
```

Let `nu_(j,c,i)` be the restored Smith orders of this Laurent matrix
in cut `c`; negative orders are allowed. Define

```math
\begin{aligned}
U_j(s)&=\sum_i s^{-\nu_{j,0,i}},&
V_j(r)&=\sum_i r^{-\nu_{j,1,i}},\\
A_{j,c}&=\langle pi_X,h_{j,c,x}+h_{j,c,v}\rangle
       +\langle pi_Y,h_{j,c,y}+h_{j,c,u}\rangle .
\end{aligned}                                                     \tag{4}
```

Choose nonnegative rational allowances `d_j` with sum `D`, and real
scales `s_j,r_j >= 1`. Put

```math
F_j=\log U_j(s_j)+\log V_j(r_j)
       +(d_j+A_{j,0})\log s_j+(d_j+A_{j,1})\log r_j.          \tag{5}
```

If a rational `1/2 <= rho < 1` satisfies

```math
\max_j F_j<2rho H,                                           \tag{6}
```

the construction supplies effective square crossing pieces with that
width parameter. The pair-swap construction and bit-cost transfer are
those of [independent paired types](paired-type-crossing-pressure.md).

For a proof, take tensor powers of admissible length `m` and restrict
the paired legs to the fixed types. A height in (3) then has the same
sum `m A_(j,c)` on every retained entry. Thus an original Hasse jet
of index `k` is exactly a coefficient of index `k+m A_(j,c)` in the
gauged tensor. It does not become a symbol-dependent linear combination
of jets. This scalar shift is valid at a finite place as well as at
infinity.

Use Hermite multiplicities `m d_j+O(1)` totaling `Dm+1`. The local
Smith coefficient estimate gives

```math
r(A_{j,k})c(A_{j,k})\le poly(m)\exp(mF_j).                  \tag{7}
```

The polynomial prefactor also covers negative local orders: clear a
fixed scalar power of `t` in the seed, apply the ordinary estimate,
and restore the corresponding coefficient shift. The fixed clearing
power contributes only `O(m)` to the coefficient range.

Hermite reconstruction expresses the crossing as the sum of these
`O(m)` jets. Tensor with the pair-swapped copy, keeping its jet index
separate. The resulting equal alphabet has size
`exp(mH+O(log(m+1)))`; both cut ranks multiply. Equation (6), with
its strict margin, gives the required width bound. The complete tables,
Hermite coefficients, and rank factorizations can be constructed in
`exp(O(m))=2^O(lambda)` bit operations. All seed data are fixed.

## 2. A profile inequality for a constant restriction

Let `M(t)` be a polynomial matrix over a field. Let `B(t)` be a
submatrix selected by constant row and column restrictions. More
generally, constant left and right linear maps are allowed. Suppose
the resulting matrix is square and generically invertible.
For every `s>=1`, its local Smith profile satisfies

```math
\sum_i s^{-\nu_i(B)}\le\sum_i s^{-\nu_i(M)}.                \tag{8}
```

Zero generic channels of `M`, if any, contribute zero. Here is a
direct proof that does not assume an individual invariant-factor
interlacing formula. At precision `L`, multiplication by `M` on
truncated power series is represented by its lower block-Toeplitz
matrix. Applying the constant restrictions on both sides cannot
increase its rank. For nonnegative Smith orders,

```math
R_M(L)=\sum_i(L-\nu_i(M))_+ .                              \tag{9}
```

For `x=1/s` with `0<x<1`, summing the rank inequalities and using

```math
(1-x)^2\sum_{L>=1}(L-a)_+x^{L-1}=x^a                     \tag{10}
```

proves (8). The case `s=1` is the generic-rank inequality. A common
scalar clearing power proves the same result for Laurent matrices.

In particular, let `P_m` be the complete type-restricted tensor power
of `K`. Its local unweighted profiles in the two cuts have partition
sums `U_(j,m)^type(s)` and `V_(j,m)^type(r)`. Applying (8) to the
gauged tensor power and restoring its scalar type shift gives

```math
\begin{aligned}
U_{j,m}^{type}(s)&\le s^{mA_{j,0}}U_j(s)^m,\\
V_{j,m}^{type}(r)&\le r^{mA_{j,1}}V_j(r)^m.
\end{aligned}                                                     \tag{11}
```

These inequalities include all physical type costs. Omitting either
scalar shift can reverse the claimed bound.

## 3. The `3/4` floor survives arbitrary physical types

The two type-class sizes are `D_X(m),D_Y(m)`. Put

```math
n_m=D_X(m)D_Y(m),\qquad \log n_m=mH+O(\log(m+1)).            \tag{12}
```

The complete tensor `P_m` is a polynomial crossing of degree at most
`Dm`, with both cut matrices of size `n_m`. Apply the unweighted
same-place floor to this tensor itself, with allowances `m d_j`.
That proof only needs two square cut matrices of the same size and
the common entry-degree budget; equal physical alphabet sizes are
unnecessary. Homogeneous scalar content and an oversized degree bound
are already included in its proof.

For

```math
f_n=\log(2n)-\frac{n+1}{2n}\log(n+1),
```

the unweighted floor and (11) imply

```math
m\max_jF_j\ge\log n_m+f_{n_m}.                            \tag{13}
```

As `m` tends to infinity through admissible lengths,
`f_(n_m)/log n_m` tends to `1/2`. Divide (13) by `m` to obtain

```math
\boxed{\max_jF_j\ge\frac32 H.}                            \tag{14}
```

Consequently the strict inequality (6) requires `rho>3/4`.
Independent type choices, independent physical heights, independent
cut scales, or additional same-place jets cannot make this criterion
reach `1/sqrt(2)`.

This does not extend the floor to arbitrary upper bounds on actual
jet ranks. The argument compares Smith-profile estimates with the
complete type-restricted Smith profiles. A stronger estimate of the
coefficient ranks could lie below both.

## 4. Exact controls

Run

```sh
python3 tools/check_weighted_multipoint_hermite.py \
  --output evidence/weighted-multipoint-hermite.json
```

The checker builds complete finite type restrictions, computes their
local profiles by independent block-Toeplitz ranks, and compares (11)
at rational scales. It includes negative heights, unequal paired type
classes, finite places, infinity, and a negative control that omits
the scalar type shift. These checks support the transcription of the
uniform argument; they are not a stronger factoring bound.

The construction builds on the repository's Hermite and paired-type
transfers. The present proof and checker used OpenAI assistant support.
