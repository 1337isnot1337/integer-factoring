# A Smith-certificate floor for sum-preserving pair conjugations

## Result and scope

Conjugating the binary elementary pencil by an invertible matrix on
ordered pairs can change its column ranks while keeping every row rank.
The basic Smith certificate nevertheless cannot improve the half-rate
factoring bound if the matrix preserves each pair-sum block and commutes
with exchange of the two labels.

For this whole family, at every positive binary width, the certificate has

```math
\rho\ge \frac14\log_2(64/9)
       =\frac32-\frac12\log_2 3
       =0.7075187496\ldots>\frac1{\sqrt2}.
```

The statement limits the sufficient two-spectrum envelope in
[the Smith transfer](smith-profile-crossing-transfer.md). It is not a
lower bound on actual coefficient rank sums. The allowed block matrices
are arbitrary invertible matrices, including pair permutations. Changing
the polynomial parameter separately in different blocks is outside the
statement.

## 1. The conjugated family

Work over any field of characteristic two. Let `b>=1`, `q=2^b`,
`n=q^2`, and label the four physical legs by `G=F_2^b`. Define

```math
p(z)=I_4+z(E_2+I_4),\qquad P_b(z)=p(z)^{\otimes b}.
```

For each `g in G`, the row block with `x+y=g` has free coordinate `x`.
Let `T_g` be the permutation `x -> x+g`. Choose any constant invertible
`q by q` matrix `S_g` that commutes with `T_g`, and put

```math
S=\bigoplus_{g\in G}S_g,\qquad S_gT_g=T_gS_g.      \tag{1}
```

A special case chooses a permutation `phi_g` with
`phi_g(x+g)=phi_g(x)+g`, giving
`S(x,x+g)=(phi_g(x),phi_g(x)+g)`. Put

```math
Q_b(z)=S P_b(z)S^{-1}.
```

This is a polynomial of degree at most `b`, with `Q_b(0)=I` and
`Q_b(1)=E_q`. Both of its tensor cuts have full generic rank.

The block structure in (1) preserves conservation support
`x+y=u+v`. It does not require the resulting entries to depend only on
`x+y` and `x+u`. Thus this family can leave the smaller convolution-kernel
space even though the following obstruction still applies to its Smith
certificate.

## 2. The unchanged row spectrum

The one-bit row degree spectrum at infinity is `(1,0,0,-1)`.
Tensor products add Smith orders: the tensor products of the invertible
outside factors remain invertible over the valuation ring. Conjugation
by any constant invertible matrix also preserves all Smith orders. Thus
the row Laurent sum is exactly

```math
U(s)=(s+2+s^{-1})^b,\qquad s\ge1.                   \tag{2}
```

## 3. One forced low-order channel in each column block

Arrange the column cut in blocks indexed by `w=x+u=y+v`, using `x,y`
as the free coordinates inside a block. Write this `q by q` block as

```math
C_w(z)_{x,y}=Q_b(z;x,y,x+w,y+w).                    \tag{3}
```

Each block is full rank at `z=1`, since at that endpoint it is the
permutation `y=x+w`.

The entire original row block indexed by `g` has degree at most
`|g|`, where `|g|` is Hamming weight. At an equal bit its two-state
block is the identity; only unequal bits contribute a degree-one factor.
Constant block conjugation preserves this degree bound. At `z=1`,
commutation with `T_g` preserves the exchange endpoint. Thus

```math
f_{g,x}(z):=Q_b(z;x,x+g,x+g,x),\qquad
\deg f_{g,x}\le |g|,\qquad f_{g,x}(1)=1.            \tag{4}
```

Set `g=w` in (4). In every column block,

```math
C_w(z)_{x,x+w}=f_{w,x}(z)\ne0,\qquad
\operatorname{ord}_0 f_{w,x}\le |w|.                \tag{5}
```

The smallest Smith order of a polynomial matrix is the smallest
valuation of any entry. If `nu_{w,1},...,nu_{w,q}` are the orders of
`C_w`, (5) gives `min_j nu_{w,j}<=|w|`. Thus the full column sum obeys

```math
V(s)=\sum_{w,j}s^{-\nu_{w,j}}
 \ge\sum_{w\in G}s^{-|w|}
 =(1+s^{-1})^b.                                   \tag{6}
```

Independently, the determinant of the complete `n by n` column matrix
has degree at most `nb`. Its nonzero determinant therefore gives
`sum_{w,j}nu_{w,j}<=nb`. Arithmetic--geometric mean yields

```math
V(s)\ge n\,s^{-b}=(4/s)^b.                         \tag{7}
```

For pair permutations the stronger equality `f_{w,x}(z)=z^{|w|}`
also holds. It is unnecessary: the degree and nonzero endpoint in (4)
already prove the valuation bound for arbitrary block matrices.

These two lower bounds need not be attained by the same spectra. Both
hold for every allowed conjugation, which is enough for their maximum.

## 4. Exact minimization and the strict threshold

Equations (2), (6), and (7) show

```math
B(s)^{1/b}\ge
(s+2+s^{-1})\max(1+s^{-1},4/s).                    \tag{8}
```

For `1<=s<=3`, the second term in the maximum is larger, and the
right side is `4(1+1/s)^2`, a decreasing function. For `s>=3`, the
first is larger, and the right side is `(s+1)^3/s^2`. Its logarithmic
derivative is `(s-2)/(s(s+1))`, so it increases on this interval.
Both values at `s=3` equal `64/9`. Therefore

```math
B(s)\ge(64/9)^b,\qquad
\rho=\frac{\log B(s)}{4\log q}
\ge\frac14\log_2(64/9).                            \tag{9}
```

There is a short exact integer certificate for the strict comparison.
The inequalities

```math
2^{317}>3^{200},\qquad 2\cdot283^2>400^2
```

give, respectively,

```math
\frac14\log_2(64/9)>\frac{283}{400}
>\frac1{\sqrt2}.
```

Exchanging the two Smith cuts cannot help. The row matrix at zero is
the identity, so its zero orders all vanish. Its corresponding Laurent
sum is `n`. The sum of the column degree spectrum is the degree of its
nonzero determinant, which is nonnegative. Arithmetic--geometric mean
then makes the other Laurent sum at least `n`, giving `B>=n^2` and
`rho>=1` for this orientation.

## Verification and attribution

Run `python3 tools/check_sum_preserving_smith_barrier.py`. The checker
generates independent sum-preserving, exchange-commuting permutations
and invertible block matrices;
constructs their exact polynomial tensors; checks both endpoints and all
forced entries in (5); and computes both Smith spectra over `F_2`.
It also checks selected Smith orders independently with block-Toeplitz
ranks, exact rational instances of (6)--(8), and the two integer
threshold comparisons. The proof establishes the all-width statement;
these finite checks are controls.

The family and obstruction were developed with substantial OpenAI
assistant support during the search for stronger crossing constructions.
No Lean formalization or improved factoring exponent is claimed.
