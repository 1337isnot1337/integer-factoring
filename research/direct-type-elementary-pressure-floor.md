# Every direct binary type of the elementary crossing stays above the half-rate target

## Result and scope

Take a tensor power of the elementary binary crossing, then restrict all
four physical legs to the full Hamming layer of one fixed weight.
For every block length and every nontrivial weight, its two-cut Smith
pressure satisfies

$$\boxed{\rho\ge\frac34>\frac1{\sqrt2}.}\tag{1}$$

This is a bound after the physical type restriction. It does not replace
the restricted tensor by a bound for the original full word alphabet.
It rules out this entire direct type family as a way to pass the current
half-rate threshold through the stated Smith evaluator.

The theorem concerns the elementary binary seed and full common Hamming
layers. It does not cover other seeds, further codes inside a layer,
different physical restrictions on different legs, or added auxiliary
constructions. It is not a lower bound on integer factoring.

The proof needs only one unit Smith factor from each explicit Johnson
block. A complete formula for the other Smith factors is unnecessary.
An additional exact characteristic-two identity in Section 5 explains
the symmetry observed in the balanced blocks.

## 1. Seed, types, and pressure

All polynomial matrices are over `F_2[z]`. For binary labels define

$$C_{xyuv}=[x=u][y=v],\qquad E_{xyuv}=[x=v][y=u],$$

and use

$$T(z)=C+z(E+C).\tag{2}$$

At `z=1` this is the exact crossing `E`. Let `T_(t,k)` be its
`t`-fold tensor power restricted on each of its four legs to words
of weight `k`. Its physical alphabet has size

$$q=\binom tk.$$

Bit complementation exchanges `k` and `t-k`, so assume
`1<=k<=t/2`. The polynomial degree is `d=2k`: only positions at
which the first two words differ contribute a positive degree, and
two weight-`k` words can differ in at most `2k` positions.

Use the row cut `(x,y)|(u,v)` and column cut `(x,u)|(y,v)`.
For a matrix of degree at most `d`, its infinity profile means the
Smith orders at zero of its reversal `z^d M(1/z)`.
If `e_i` is the row infinity profile and `f_i` the column zero
profile, the [Smith pressure](block-smith-crossings.md) at `s>=1` is

$$P(s)=\left(\sum_i s^{d-e_i}\right)
       \left(\sum_i s^{-f_i}\right),\qquad
\rho(s)=\frac{\log P(s)}{4\log q}.\tag{3}$$

All these profiles have `q^2` finite entries. Indeed, the endpoint
at `z=1` is a permutation matrix in either cut.

## 2. Explicit row and column blocks

For a pair of weight-`k` words `x,y`, suppose their common ones
occupy `k-j` positions and their symmetric difference has size `2j`.
The positions and values where they agree are preserved by every
nonzero entry of the row matrix. Thus the row matrix splits into
blocks indexed by those fixed positions. The number of blocks is

$$R_j=\frac{t!}{(k-j)!(t-k-j)!(2j)!},\qquad 0\le j\le k.\tag{4}$$

Inside one block, `x` and `u` specify `j`-subsets `S,U` of the
`2j` differing positions. If `a=|S\setminus U|`, the entry is

$$z^{2a}(1+z)^{2j-2a}.$$

After reversal at the global degree `2k`, this becomes

$$z^{2(k-j)}(1+z^2)^{|S\cap U|}.\tag{5}$$

At zero, the matrix remaining after the displayed monomial is the
all-ones matrix. It has exactly one unit Smith factor. Therefore
each of the `R_j` blocks contributes an infinity order `2(k-j)`.
In the first sum of (3), these factors alone give

$$\sum_i s^{2k-e_i}\ge\sum_{j=0}^k R_j s^{2j}.\tag{6}$$

For the column cut, first fix `D=x\mathbin\triangle u`, of size
`2j`, and the orientation of the `j` ones of `x` inside `D`.
At these positions, the nonzero entry forces `y=u`, `v=x`, and
contributes `z^(2j)`. Outside `D`, both `x=u` and `y=v` have
`k-j` ones. The remaining matrix has entries

$$z^{2j}(1+z^2)^{|S\setminus U|},\tag{7}$$

indexed by `(k-j)`-subsets of a `(t-2j)`-element set.
The number of blocks is

$$B_j=\binom t{2j}\binom{2j}j
    =\frac{t!}{j!^2(t-2j)!}.\tag{8}$$

After removing `z^(2j)`, each block again has an all-ones constant
matrix and exactly one unit factor. Hence

$$\sum_i s^{-f_i}\ge\sum_{j=0}^k B_j s^{-2j}.\tag{9}$$

These are exact block decompositions. They retain the full coefficient
matrix on the chosen physical alphabet, including all cancellations
inside each block.

## 3. A binomial inequality gives the uniform floor

Put

$$p_j=\binom kj\binom{t-k}j.$$

The block counts satisfy

$$R_j B_j=q\,p_j\binom t{2j}.\tag{10}$$

Vandermonde's identity gives both

$$\binom t{2j}\ge\binom kj\binom{t-k}j=p_j,
\qquad \sum_{j=0}^k p_j=\binom tk=q.\tag{11}$$

Apply Cauchy--Schwarz to the two positive sums in (6) and (9):

$$\begin{aligned}
P(s)
&\ge \left(\sum_j R_j s^{2j}\right)
      \left(\sum_j B_j s^{-2j}\right)\\
&\ge \left(\sum_j\sqrt{R_jB_j}\right)^2\\
&\ge \left(\sqrt q\sum_j p_j\right)^2
 =q^3.\tag{12}
\end{aligned}$$

This proves `rho(s)>=3/4` at every allowed evaluation point, without
optimizing a numerical parameter. Since `(3/4)^2=9/16>1/2`, it is
strictly above the half-rate target.

The opposite assignment of the two cuts cannot evade the conclusion.
The row matrix at zero is the identity, so its zero profile consists
of `q^2` zeros. For any nonsingular polynomial matrix `M` of degree
at most `d`, the reversed Smith orders `a_i` satisfy

$$\sum_i a_i=dq^2-\deg\det M\le dq^2.$$

Thus the arithmetic--geometric mean inequality gives
`sum_i s^(d-a_i)>=q^2` for `s>=1`. The opposite pressure is at
least `q^4`, and its normalized value is at least one.

## 4. Balanced types have identical good profiles

When `t=2k`, define the middle-slice matrix

$$H_j(w)_{S,U}=(1+w)^{|S\cap U|},\qquad
 |S|=|U|=j,\quad S,U\subseteq[2j],\tag{13}$$

with `H_0=[1]`. The reversed row blocks in (5) are exactly

$$z^{2(k-j)}H_j(z^2).$$

In the column block (7), the outside ground set has size `2(k-j)`.
Complementing its column label changes the distance exponent into
the intersection exponent. Thus its profile is that of
`z^(2j) H_(k-j)(z^2)`. Moreover, `B_j=R_(k-j)`.

The row infinity and column zero profiles are therefore identical
for every balanced type. This identity is stronger than an equality
of determinants or average orders.

## 5. Exact middle-slice duality in characteristic two

For every `j>=0`, the matrix in (13) satisfies

$$\boxed{H_j(w)^2=w^{2j}I.}\tag{14}$$

Here is a direct proof that uses no rational eigenbasis. Fix two
`j`-subsets `S,U` and write `a=|S\cap U|`. Summing over the middle
subset gives the coefficient

$$[x^j]\,(1+x(1+w)^2)^a
 (1+x(1+w))^{2(j-a)}(1+x)^a.$$

In characteristic two this is

$$[x^j]\,(1+xw^2+x^2(1+w)^2)^a
 (1+x^2(1+w)^2)^{j-a}.\tag{15}$$

Choose `h` linear terms `xw^2` from the first `a` factors. The
other `j-h` factors contribute through `1+x^2(1+w)^2`.
When `j-h` is positive and even, the required central coefficient
is `binom(j-h,(j-h)/2)`, which is even. When it is odd, the
coefficient is zero. Only `h=j` survives; this is possible precisely
when `a=j`, or `S=U`. Its value is `w^(2j)`, proving (14).

Let the Smith orders of `H_j` over `F_2[[w]]` be
`h_1<=...<=h_m`, where `m=binom(2j,j)`. Equation (14) implies

$$h_i+h_{m+1-i}=2j.\tag{16}$$

Indeed, `H_j^(-1)=w^(-2j)H_j`; compare the orders of the inverse
with those of the scalar multiple. In particular every order lies
between zero and `2j`, and their sum is `jm`.

The finite profiles through `j=6` have an additional simple pattern:
an order `i<j` occurs `binom(2j-1,i)-binom(2j-1,i-1)` times, paired
with order `2j-i`, and order `j` is absent. That complete multiplicity
formula is not proved here and is not used in (1).

Characteristic zero diagonalization is insufficient for this question.
Already `H_1=[[1+w,1],[1,1+w]]` has orders `(0,2)` over `F_2[[w]]`.
The usual sum-and-difference eigenbasis has determinant two and is
not an invertible basis in this ring.

## 6. Exact replay and a source boundary

Run the [checker](../tools/check_direct_type_elementary_floor.py):

```sh
python3 tools/check_direct_type_elementary_floor.py
```

It verifies the restricted tensors and both block decompositions entry
by entry on small types, computes their exact Smith profiles, checks
the rational pressure lower bounds, and checks the integer identities
in (10)--(11) across a larger grid. It separately replays (14) and
the finite middle-slice profiles. The complete counts are saved in
[`direct-type-elementary-floor.json`](../evidence/direct-type-elementary-floor.json).

For continued work on the remaining multiplicity formula, the integral
Bier basis is a relevant primary tool. Ducey et al.,
[*Integer diagonal forms for subset intersection relations*,
Theorem 1.7](https://arxiv.org/html/2310.09227v2),
give its triangular action on intersection matrices. Their later
small-block diagonalization requires `n>=3k-1` and cannot be imported
directly for the middle layer `n=2k`. The elementary proofs above
do not rely on that extension.
