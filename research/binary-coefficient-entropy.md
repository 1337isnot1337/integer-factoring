# Exact ranks of the elementary binary coefficient family

## Result

The full binary tensor-power family has a precise limiting rank rate.
It is larger than the general support lower bound proved in
[coordinate-support-rank.md](coordinate-support-rank.md).

Over `F_2`, let

```math
E(x,y,u,v)=[x=v][y=u],\qquad
R(x,y,u,v)=[x=y][v=u],\qquad C=E+R,
```

where each index is one bit. For `t`-bit labels, define

```math
A_i=[z^i](R+zC)^{\otimes t},\qquad 0\le i\le t.
\tag{1}
```

The pieces sum to the crossing. Let `r_i,c_i` be their ranks in the
row `(x,y)|(v,u)` and column `(x,u)|(v,y)` flattenings. Define

```math
H_2(a)=-a\log_2a-(1-a)\log_2(1-a),
```

with the endpoint values defined by continuity. Then

```math
\begin{aligned}
\log_2 r_i
 &=2tH_2\!\left(\frac{i}{2t}\right)+O(\log(t+2)),\\
c_i&=r_{t-i},\\
\max_i\log_2(r_ic_i)
 &=4tH_2(1/4)+O(\log(t+2)).
\end{aligned}
\tag{2}
```

The error is uniform in `0<=i<=t`. Thus the best leading coefficient
for these literal pieces is

```math
\rho=H_2(1/4)=0.8112781244\ldots.
\tag{3}
```

This rules out extrapolating their attractive very small examples to an
improvement on `1/sqrt(2)`. It does not give the same `0.811...` lower
bound after arbitrary coefficient recombinations or restrictions of the
label set. The separate support theorem gives `rho>=3/4` for those
broader operations under its stated piece-count hypothesis.

The only imported rank result below is Wilson's inclusion-matrix formula.
A primary source with a proof is Peter Frankl,
[*Intersection Theorems and Mod p Rank of Inclusion Matrices*,
Theorem 4.1](https://www.renyi.hu/~pfrankl/j62.pdf), JCTA 54 (1990).
The application to the crossing coefficients and the entropy estimate
are proved here.

## 1. A nilpotent block calculation

In the row flattening, `E` is the identity of order four. The matrix
`R` has rank one and satisfies `R^2=0`, because its only nonzero
two-by-two block is the all-ones matrix over `F_2`.
It is therefore similar to one nilpotent Jordan block of order two
and two zero blocks of order one. Consequently

```math
R+zC=zI+(1+z)R
\sim \bigl(zI_2+(1+z)N\bigr)\oplus z\oplus z,
\qquad N^2=0,\quad\operatorname{rank}N=1.
\tag{4}
```

Tensoring the fixed changes of basis preserves every coefficient rank.
Choose the two-dimensional block in `k` of the `t` positions. There
are `binom(t,k)2^(t-k)` copies of the resulting block. Identify its
underlying vector space with the regular module of

```math
B_k=\mathbf F_2[x_1,\ldots,x_k]/(x_1^2,\ldots,x_k^2).
```

Each `N` is multiplication by one variable. The coefficient of order
`i` in that block is multiplication by

```math
\sum_{j\ge t-i}\binom{j}{t-i}e_j(x_1,\ldots,x_k)
=\left(\prod_{a=1}^k(1+x_a)\right)e_{t-i}(x_1,\ldots,x_k).
\tag{5}
```

Here `e_j` is an elementary symmetric polynomial and is zero for
`j>k`. The equality follows by extracting a coefficient from
`prod_a(1+(1+w)x_a)` and using `x_a/(1+x_a)=x_a` in `B_k`.
The prefactor is a unit. If `u(k,r)` denotes the rank of multiplication
by `e_r` on `B_k`, equation (5) gives the exact formula

```math
r_i=\sum_{k=t-i}^{t}\binom tk2^{t-k}u(k,t-i).
\tag{6}
```

The column flattening of `R` is a permutation matrix `P` consisting
of one transposition and two fixed points. Put `T=I+P`; it again
has rank one and square zero. In that flattening,

```math
R+zC=I+(1+z)T=z\bigl(z^{-1}I+(1+z^{-1})T\bigr).
```

The row and column expressions have similar nilpotent parts. Reversing
the polynomial after taking the tensor power proves `c_i=r_(t-i)`.

## 2. Exact inclusion ranks

Multiplication by `e_r` sends degree `j` to degree `j+r` in `B_k`.
Its matrix is the zero-one inclusion matrix between `j`-subsets and
`(j+r)`-subsets of a `k`-set. Write

```math
b=\min(j,k-r-j).
```

Complementation and transposition reduce the matrix to the usual range
`b<=k-(b+r)`. Wilson's formula gives

```math
\operatorname{rank}_{\mathbf F_2}W_{j,j+r}(k)
=\sum_{a=0}^{b}
 \left[\binom{b+r-a}{r}\text{ is odd}\right]
 \left(\binom ka-\binom k{a-1}\right).
\tag{7}
```

Set `binom(k,-1)=0`. The spaces of different input degree have disjoint
output degree, so

```math
u(k,r)=\sum_{j=0}^{k-r}
\operatorname{rank}_{\mathbf F_2}W_{j,j+r}(k).
\tag{8}
```

Equations (6)--(8) are an exact finite formula, including the small
characteristic cancellations. They require only integer binomial
coefficients and parity tests, with no large matrix construction.

For a simpler estimate, the last term `a=b` in (7) always survives.
The matrix dimension supplies the upper bound. Thus

```math
\frac1{k+1}\binom kb
\le \binom kb-\binom k{b-1}
\le\operatorname{rank}W_{j,j+r}(k)
\le\binom kb.
\tag{9}
```

The first inequality uses
`binom(k,b)-binom(k,b-1)=binom(k,b)(k-2b+1)/(k-b+1)`
and `b<=k/2`.

## 3. The binomial tail and its entropy

Define

```math
B(t,i)=\sum_{a=0}^{i}\binom{2t}{a}.
```

Combining (6), (8), and (9) gives the useful integer bounds

```math
\boxed{\frac{B(t,i)}{t+1}\le r_i\le2B(t,i).}
\tag{10}
```

To see the counting step explicitly, each value of
`min(j,k-r-j)` occurs once or twice as `j` ranges. Put `r=t-i`
and `a=t-k`. The one-copy sum in (6) is

```math
\sum_{\substack{a,b\ge0\\a+2b\le i}}
 \binom ta2^a\binom{t-a}{b}
=\sum_{v=0}^{i}[z^v](1+2z+z^2)^t
=B(t,i).
\tag{11}
```

Replacing `1/(k+1)` by `1/(t+1)` proves the lower side of (10).
The possible second copy proves its upper side.

For `0<=i<=t`, the elementary binomial entropy bounds imply

```math
\log_2B(t,i)=2tH_2(i/(2t))+O(\log(t+2)).
```

For example, the largest summand lies between
`2^(2tH_2(i/(2t)))/(2t+1)` and `2^(2tH_2(i/(2t)))`,
and there are at most `2t+1` summands. This proves the first line of (2).
The function
`H_2(a)+H_2(1/2-a)` is concave and symmetric about `a=1/4`.
Its maximum there, with the nearest integer coefficient for odd `t`,
proves the last line of (2).

## 4. Checks and scope

Run `python3 tools/check_binary_coefficient_entropy.py`.
The checker compares (7) with direct finite-field elimination, compares
(6) with complete small crossing tensors, and checks (10) using exact
integers over a wider parameter range. It also records large coefficient
profiles without constructing their exponentially large matrices.

The [symmetric three-matching theorem](symmetric-three-matching-barrier.md)
extends the same leading lower bound to arbitrary recombinations in
the full permutation-symmetric span of all three local pairings.
That extension needs two further rank filtrations and a walk-count
argument; it does not follow from (2) alone. The proofs used substantial
OpenAI assistant support and no new Lean formalization.
