# The ideal Hankel row profile cannot have a rank-one opposite endpoint

## Result and conventions

Work over any field of characteristic two. Let `d>=3`, `N=d^2`, and
index `N` by `0,...,N-1`. Define the polynomial matrix

```math
H(z)_{ij}=z^{i+j}[i+j<N],\qquad D=N-1.                    \tag{1}
```

For an arbitrary constant invertible `N` by `N` matrix `U`, consider

```math
F(z)=U H(z)H(1)^{-1}U^{-1}.                              \tag{2}
```

Interpret its entries as a four-label tensor, using the row cut
`(x,y)|(v,u)`, where every label ranges over `d` values. Its other cut
is `(x,u)|(v,y)`. Since `F(1)=I_N`, the coefficient tensors sum to
the crossing `[x=v][y=u]`.

The row Smith valuations of (1), and hence (2), are exactly

```math
0,\ \underbrace{N,\ldots,N}_{N-1\text{ times}}.          \tag{3}
```

However, the leading coefficient `[z^D]F(z)` has **opposite-cut rank
at least two**. In particular the reversed opposite-cut polynomial
cannot have only one zero Smith valuation. Thus this ideal Hankel
realization of the row determinant/Jensen profile cannot also realize
a rank-one opposite endpoint.

The statement covers every **constant** left/right normalization of
`H` whose value at one is the crossing: if `A H(1) B=I`, then
`B=H(1)^(-1)A^(-1)`, so that normalization is (2). It does not cover
parameter-dependent gauges, arbitrary matrices with the same row
valuation list, or general crossing constructions. It is not by itself
a lower bound on the best attainable Smith envelope.

The restriction `d>=3` is necessary for this endpoint argument.
Section 4 gives a rank-one opposite endpoint at `d=2`.

## 1. The exact row valuations

The upper-left entry of `H(z)` is one. Eliminate the other entries in
its first row and column by polynomial elementary operations. On the
remaining indices `1<=i,j<N`, the Schur complement is

```math
z^{i+j}[i+j<N]-z^{i+j}
 =-z^{i+j}[i+j\ge N].                                  \tag{4}
```

It is `z^N C(z)`, where the constant matrix `C(0)` is minus the
reversal matrix of size `N-1`. Hence `C(z)` is invertible over the
formal power series ring. These unit operations give (3).

In particular the sum of the valuations is `N(N-1)=ND`, the largest
value permitted by the determinant degree. All the nonzero valuations
are equal. This is the precise sense in which (1) realizes the ideal
row determinant/Jensen profile. Constant invertible multiplication
does not change those valuations.

## 2. Normalization makes the leading coefficient regular unipotent

Let `R` be the reversal matrix of size `N`, and let `S` be the lower
shift with ones immediately below the diagonal. Then

```math
H(1)=R(I-S)^{-1},\qquad [z^D]H(z)=R.                    \tag{5}
```

The inverse in (5) is the finite sum `I+S+...+S^(N-1)`.
Consequently

```math
[z^D]\bigl(H(z)H(1)^{-1}\bigr)
 =R(I-S)R=I-S^{\mathsf T}.                              \tag{6}
```

This is one regular unipotent Jordan block of size `N`.
Conjugating by `U` preserves its Jordan type and its multiplicative
order. In characteristic two that order is exactly

```math
2^{\lceil\log_2 N\rceil}.                              \tag{7}
```

Indeed, `(I-S^T)^(2^a)=I+(S^T)^(2^a)`, and the shift has nilpotency
index exactly `N`.

## 3. An invertible rank-one opposite cut has smaller unipotent order

Let `L` be any invertible row matrix whose opposite-cut rank is one.
There are two `d` by `d` matrices `P,Q` with

```math
L_{(x,y),(v,u)}=P_{xu}Q_{yv},\qquad
L=(P\otimes Q)\,\mathsf{Swap},                          \tag{8}
```

where `Swap(e_v tensor e_u)=e_u tensor e_v`.
Conversely, (8) has opposite-cut rank one when `P,Q` are nonzero.
Since `L` is invertible, both `P` and `Q` are invertible. Squaring gives

```math
L^2=PQ\otimes QP.                                      \tag{9}
```

Suppose also that `L` is unipotent. Choose a power of two `2^a`
with `L^(2^(a+1))=I`. Equation (9) gives

```math
(PQ)^{2^a}\otimes(QP)^{2^a}=I.                          \tag{10}
```

If two square matrices have tensor product equal to the identity,
they are reciprocal nonzero scalar matrices. This follows directly
from their matrix blocks: off-diagonal blocks force the first matrix
to be diagonal, and its diagonal blocks force the second to be scalar.
The matrices in (10) are also similar, since `PQ` and `QP` are similar.
Their two scalars must therefore be equal. Their square is one, and
in characteristic two the only such scalar is one.

It follows that `(PQ)^(2^a)=I`, so `PQ` is unipotent. Its nilpotent
part has index at most `d`; therefore

```math
(PQ)^{2^{\lceil\log_2 d\rceil}}=I.
```

The same holds for `QP`. Using (9) once more shows

```math
\operatorname{ord}(L)
\le 2^{1+\lceil\log_2 d\rceil}.                         \tag{11}
```

For every integer `d>=3`,

```math
d^2>2^{1+\lceil\log_2 d\rceil}.
```

At `d=3`, this says `9>8`. For `d>=4`, use
`2^ceil(log2 d)<2d` and `d^2>=4d`. Hence (7), with `N=d^2`,
is strictly larger than (11). A regular unipotent leading coefficient
from (6) cannot have rank one in the opposite cut. This proves the
result.

## 4. The two-label endpoint really can occur

For `d=2`, let

```math
P=I_2,\qquad Q=\begin{pmatrix}1&1\\0&1\end{pmatrix},
\qquad L=(P\otimes Q)\mathsf{Swap}.
```

Its opposite-cut rank is one by (8). Put `T=L-I_4`.
In the ordered basis `e_00,e_01,e_10,e_11`,

```math
T^3e_{11}=e_{00}\ne0,\qquad T^4=0.
```

Thus `L` is a regular unipotent block of size four. More explicitly,
the columns

```math
e_{00},\quad e_{00}+e_{01}+e_{10},\quad e_{10},\quad e_{11}
```

form an invertible matrix `U` with
`L=U(I_4-S^T)U^(-1)`. Applying this `U` in (2) gives an actual
normalized Hankel curve whose leading opposite cut has rank one.
This example verifies the boundary of the order argument; it does
not claim a favorable complete Smith envelope.

## Verification

The checker verifies (3), (5)--(6), the normalization conventions,
the two-label endpoint example, random twisted-swap unipotent order
bounds, and the exact order separation for many integer dimensions.
All matrix operations use exact binary elimination.

```sh
python3 tools/check_ideal_hankel_endpoint.py \
  --output evidence/ideal-hankel-endpoint.json
```

These finite checks test the algebraic transcription. The proof above
establishes the claim for all `d>=3` and all fields of characteristic two.
