# Rank-sensitive Fourier contraction with digit restrictions

## Result and limits

This note extends the [arithmetic Fourier blocks](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/logarithmic-wire-network.md)
on PR #23, pinned at commit `147652574e7ce1ff68b2c378712846cbf438da34`.
One complete input-digit family can be included if the opposite input
family is left unrestricted at the block boundary. The resulting query
is a weighted digit-box sum on an explicitly computable affine subspace.
Both its rank and its nullity give algorithms.

The rank computations can be performed without storing a dense matrix.
For an `m` by `m` matrix over a prime field of characteristic greater
than `m`, a scalar-entry oracle suffices to compute determinants, rank,
and individual solution coordinates in `2^O(m)` oracle calls and
`O(m + log(m+2) log(mq+2))` additional working bits. This is useful when
`m = O(n/log n)` and the arithmetic prime `q` is polynomial in `n`.

There is no new unconditional factoring exponent here. For the actual
convolution matrix, almost all ambient output frequencies have full
Hankel rank. Even after pruning frequencies with zero carry weight,
`exp(Omega(n))` full-rank frequencies remain. Explicitly visiting those
frequencies therefore costs `exp(Omega(n))`, even if every full-rank
block query is inexpensive. The last section proves this restriction
without assuming an exact Hankel-rank counting formula.

All statements below are ordinary proofs, not Lean formalizations.

**Sharper local implementation.** The later
[scalar Berkowitz and pivot-mask construction](small-space-berkowitz-rank.md)
improves the rank/determinant query time here to
`m^O(log m) poly(log q)`, with the same linear auxiliary-space scale,
and removes the restriction on the characteristic. The trace formulas
in Section 3 remain valid as a separate elementary implementation.
The outer frequency obstruction is unchanged.

## 1. Absorb one input box into a bilinear block

Let `q` be an arithmetic prime. Let `A` be a fixed finite commutative
ring of characteristic `p`, where `p != q`, and use

```math
R=A[Z]/(1+Z+\cdots+Z^{q-1}),\qquad \zeta=Z.
```

The scalar `q` is a unit. For each `c in F_q`,

```math
\sum_{y\in\mathbb F_q}\zeta^{cy}
=q[c=0].                                                     \tag{1}
```

Indeed `zeta^q=1`; for nonzero `c`, multiplication by `c` permutes
the exponents, and their sum is zero in `R`. This proof permits zero
divisors in `A` and in `R`.

Let `x` have `s` digit coordinates and `y` have `t` unrestricted field
coordinates. Fix a matrix `J in F_q^(s by t)`, vectors `alpha,beta`,
digit sets `D_i subset F_q`, and weights `w_i:D_i -> A`. Define

```math
Z=\sum_{x\in\prod_iD_i}\sum_{y\in\mathbb F_q^t}
 \left(\prod_iw_i(x_i)\right)
 \zeta^{\alpha^Tx+\beta^Ty+x^TJy}.
```

Summing `y` first with (1) gives the exact identity

```math
Z=q^t\sum_{\substack{x\in\prod_iD_i\\J^Tx=-\beta}}
       \left(\prod_iw_i(x_i)\right)\zeta^{\alpha^Tx}.    \tag{2}
```

This is a contraction of a mixed input/arithmetic block. It does not
require the digit sets to be subgroups. All boundary-transform
normalization factors are retained separately under the same convention
as the source Fourier network.

A binary prefix restricts the base-`B`, `B=2^d`, digits to fixed digits,
at most one partial interval, and unrestricted digits in `[0,B)`.
Substitute fixed digits into (2), including their weights and constant
phase, before taking rank. Thus `s` can count only the remaining active
digits. Bit-product weights and one coordinate tag remain products of
independent digit weights, as required by (2).

The same argument applies to the bivariate network. The absorbed input
digits are first mapped by its explicit linear evaluation map `E`.
After summing the free opposite-side inputs, the affine system is
`J^T E x = -beta`. Its entries are computable from the boundary
frequencies. One must take the rank of this composed matrix, rather
than the rank of each middle gate separately.

## 2. Complementary rank and nullity algorithms

Write the affine system in (2) as `C x=b`, with rank `r`. First test
consistency. If it fails, the answer is zero. If it holds, retain
independent rows and columns `I,J` with invertible `r by r` minor.

### Nullity enumeration

Enumerate the digit values of the `s-r` nonpivot coordinates. For each
choice, solve the pivot equations by Cramer's rule, check every recovered
pivot digit against its set, and multiply its digit weight. The chosen
rows suffice because consistency was tested. Alternatively, checking
every equation again gives a direct negative control.

If all active digit sets have at most `B` elements, this uses at most
`B^(s-r)` candidates. For unequal set sizes, the exact number is the
product of the sizes of the nonpivot sets. A pivot-column basis maximizing
the product of its digit-set sizes minimizes this enumeration count;
the uniform upper bound does not depend on optimizing that basis.

Full column rank gives at most one candidate. Rank zero gives a product
of independent digit sums, which should be evaluated directly rather
than by enumerating the whole box.

### Rank enumeration

Use only the `r` independent equations, now written `C_0x=b_0`.
Orthogonality gives

```math
\sum_{C_0x=b_0}\!w(x)\zeta^{\alpha^Tx}
=q^{-r}\sum_{v\in\mathbb F_q^r}\zeta^{-v^Tb_0}
 \prod_{i=1}^s\left(\sum_{a\in D_i}
 w_i(a)\zeta^{a(\alpha_i+(C_0^Tv)_i)}\right).           \tag{3}
```

For prefix digit sets and independent bit weights, each inner sum has
a product formula. If the high bits are fixed and low bits indexed by
`L_i` are free, it is a fixed monomial and fixed weight times

```math
\prod_{j\in L_i}(1+u_{ij}\zeta^{2^jc}),
\qquad c=\alpha_i+(C_0^Tv)_i.                         \tag{4}
```

No division by `1-zeta^c` is used, so the formula remains valid when
that element is zero or a zero divisor. Equation (3) takes `q^r`
outer candidates and polynomially many elementary operations per
candidate if ring elements are materialized. Together, the two routes
give a candidate count

```math
\min\{B^{s-r},q^r\}.                                  \tag{5}
```

There is also an ordinary meet-in-the-middle option: split the digit
coordinates into two halves; enumerate each half's vector `C_0x`, phase,
and weight; aggregate equal vectors; then pair complementary vectors.
Sequential sorting avoids assuming constant-time dictionary access.
This takes `B^ceil(s/2) poly(s,q)` time and potentially exponentially
many stored signatures. It is a different memory tradeoff, not a
sublinear-workspace consequence of (5).

## 3. Determinants and rank without dense matrix storage

Here the matrix entries are supplied by a read-only oracle. Its time
and workspace are charged separately. Suppose `q>m` and `U` is an
`m by m` matrix over `F_q`. Put `p_k=tr(U^k)`. The finite-degree
Newton identity is

```math
\det U=(-1)^m
\sum_{\ell=1}^m\frac{(-1)^\ell}{\ell!}
\sum_{\substack{k_1+\cdots+k_\ell=m\\ k_i\ge1}}
 \prod_{i=1}^{\ell}\frac{p_{k_i}}{k_i}.                 \tag{6}
```

For `m=0`, return one. All displayed denominators are units because
`q>m`. To prove (6), use the formal identity through degree `m`

```math
\det(I-zU)=\exp\!\left(-\sum_{k=1}^m p_kz^k/k\right)
       \pmod{z^{m+1}}.
```

It follows, for example, by differentiating both sides as formal power
series, using the adjugate identity for the derivative of a determinant,
and matching constant terms. The coefficients `1,...,m` are invertible,
so this differential equation determines the coefficients uniquely.
Expanding the truncated exponential gives (6).

Enumerate the `2^(m-1)` ordered compositions using an `m-1` bit separator
mask. Only the current mask, field accumulators, and a few counters are
stored. For one requested `p_k`, compute matrix-power entries recursively:

```math
(U^k)_{ab}=\sum_c
 (U^{\lfloor k/2\rfloor})_{ac}
 (U^{\lceil k/2\rceil})_{cb}.
```

Use the entry oracle when `k=1`, and the identity matrix when `k=0`.
The recursion depth is `O(log(m+2))`. A frame needs
`O(log(m+2)+log q)` bits. The number of entry queries is
`m^{O(log(m+2))}` per trace, which is `2^O(m)` after increasing
an absolute constant. Thus (6) gives

```math
T_{\det}=2^{O(m)}\operatorname{poly}(\log q)\cdot(1+T_{entry}),
\qquad
S_{\det}=S_{entry}+O(m+\log(m+2)\log(mq+2)).           \tag{7}
```

Field inverses may be computed by repeated squaring to power `q-2`;
they require no factorization or large lookup table.

For a rectangular matrix with at most `m` rows and columns, enumerate
row and column subsets in decreasing common cardinality and test their
determinants. There are at most `4^m` pairs of subsets. The first
nonzero minor supplies the rank and a pivot minor. The two subset masks
have `O(m)` bits, so (7) remains valid. Consistency is checked by
comparing the ranks of `C` and `[C|b]`, taking the harmless dimension
increase into account. Cramer coordinates are ratios of two determinants
of the pivot minor, again with the same bounds.

These are deliberately exponential-in-dimension algorithms. They are
useful here because the dimension is a number of base-`B` digits, while
each digit has `Theta(log n)` bits. If `m=O(n/log(n+2))` and
`q=n^O(1)`, (7) is

```math
T_{\det}=\exp(O(n/\log(n+2)))(1+T_{entry}),
\qquad S_{\det}=S_{entry}+O(n/\log(n+2)).              \tag{8}
```

The formulas also clarify the extra storage of the mixed block:
nullity enumeration uses `(s-r) log B` digit-counter bits, while
rank enumeration uses `r log q` frequency-counter bits. The input
matrix and frequency vector remain read-only. Neither formula claims
that an arbitrary `O(n)`-bit boundary can be stored for free when the
query is later inlined into a larger evaluator.

## 4. Stream coefficients of the root algebra

For (2), a fixed digit assignment contributes a base-ring scalar times
one power `zeta^a`. In the degree-less-than-`q-1` representation, this
power has one nonzero coefficient if `a mod q < q-1`; if its residue
is `q-1`, every coefficient is `-1`. Thus a requested output
coefficient can be accumulated directly, rerunning the candidate loop.
All `q-1` coefficients may be emitted sequentially without retaining
the polynomial. This adds only a factor `q` to time.

For (3)--(4), a scalar coefficient of a balanced product can be obtained
by recursive convolution. At a product node, enumerate the two input
coefficient indices; reduce first modulo `Z^q-1`, then subtract the
coefficient at `q-1` from each lower coefficient. This is exactly
reduction modulo `1+...+Z^(q-1)`. A product with `u` factors uses
`q^{O(log(u+2))}` time and
`O(log(u+2) log(q+2))` additional bits over the fixed base ring.
In the intended polynomial parameter range this is quasipolynomial
time and polylogarithmic space, below (8).

This coefficient-by-coefficient option avoids a hidden array of `q`
root-algebra coordinates in the local-space claim. It does not make
arbitrary later arithmetic on stored root-algebra values free.

## 5. Actual convolution frequencies and their full-rank obstruction

For `m` base-`B` digits, ordinary convolution has coefficients
`c_k=sum_(i+j=k) a_i b_j`. Giving these coefficients frequencies
`gamma_0,...,gamma_(2m-2)` produces the bilinear phase

```math
\sum_k\gamma_kc_k=a^TH_\gamma b,
\qquad (H_\gamma)_{ij}=\gamma_{i+j}.                    \tag{9}
```

Thus the matrix is Hankel, with the actual arithmetic prime from the
network. The polynomial `det H_gamma` has total degree `m` and is
nonzero over every field: its monomial `gamma_(m-1)^m` has coefficient
`(-1)^(m(m-1)/2)`, coming from the unique reverse permutation.

The elementary polynomial root bound implies

```math
\#\{\gamma:\operatorname{rank}H_\gamma<m\}
 \le m q^{2m-2},\qquad
\#\{\gamma:\operatorname{rank}H_\gamma=m\}
 \ge(1-m/q)q^{2m-1}.                                 \tag{10}
```

For completeness, prove the root bound by induction on the number of
variables. Regard a nonzero degree-`D` polynomial as a univariate in
the last variable with highest nonzero coefficient of degree at most
`D-j` in the earlier variables. Outside the zeros of that coefficient
there are at most `j` roots; inside it allow all `q` values. Induction
gives at most `Dq^(k-1)` zeros in `k` variables.

The network chooses `q>m(B-1)^2`, so `q>m` and, in the logarithmic-digit
regime, the full-rank fraction in (10) tends to one. Since
`m=Theta(n/log n)` and `log q=Theta(log n)`, even just these full-rank
frequencies number `exp(Omega(n))`. A routine that explicitly visits
each ambient frequency therefore cannot prove the half-rate bound,
let alone improve it, solely by accelerating individual rank strata.

### Nonzero carry weights still leave exponentially many full-rank terms

Assume `1<=N<B^m`, as in the factoring instance. Let `F_N(c)` be the
indicator of the full legal carry chain on the `2m-1` coefficient
variables, each represented in `[0,q)`, with zero initial and final
carry. Its equalities telescope to

```math
\sum_{i=0}^{2m-2}c_i B^i=N.                            \tag{11}
```

All `c_i` are nonnegative. Consequently a legal setting has `c_i=0`
for every `i>=m`. Conversely the ordinary base-`B` digits of `N`
give a legal setting. The zero coefficient vector is illegal because
`N>=1`. Thus the induced function of the first `m` coefficient
coordinates is nonconstant.

Its Fourier transform is independent of
`gamma_m,...,gamma_(2m-2)`. By Fourier inversion using (1), a
nonconstant function cannot have its transform supported only at the
zero frequency. Therefore there is a nonzero lower-frequency vector

```math
a=(\gamma_0,\ldots,\gamma_{m-1})\ne0
```

for which the carry Fourier value is nonzero. This statement remains
valid in the root algebra with zero divisors: the Fourier matrix is
invertible because `q` is a unit. Every extension of this fixed `a`
to the remaining frequencies has the same nonzero carry value.

For any such nonzero `a`, the Hankel determinant is a nonzero
polynomial in its remaining `m-1` coordinates. Here is an explicit
certificate. Let `j` be the last index with `a_j!=0`. If `j=m-1`,
set every remaining frequency to zero. The matrix is anti-triangular
with nonzero anti-diagonal `a_j`. If `j<m-1`, set
`gamma_(m+j)=1` and all other remaining frequencies to zero.
Split the row and column indices into `[0,j]` and `[j+1,m-1]`.
The off-diagonal blocks are zero, the first diagonal block is
anti-triangular with anti-diagonal `a_j`, and the second is an
anti-identity matrix. The determinant is nonzero.

Apply the same degree-`m` root bound to these `m-1` free variables.
For `m>=2`, at least

```math
(1-m/q)q^{m-1}                                      \tag{12}
```

extensions have full Hankel rank and nonzero carry Fourier weight.
In the logarithmic-digit parameters this is `exp(Omega(n))`.
Hence a method that visits each nonzero carry-frequency term
individually still has exponential work. Faster rank-sensitive local
queries and zero-weight pruning alone do not remove that cost.

This restriction concerns explicit term-by-term Fourier enumeration.
It does not constrain the existing crossing evaluator, a collective
sum of high-rank frequencies, or further cancellation exploited by an
algorithm. It is not a lower bound for factoring.

## Verification

`check_hankel_digit_box.py` independently compares the composition
determinant with ordinary elimination, both affine-box algorithms
with direct weighted sums, and small Hankel full-rank counts with
(10). It also tests the completion certificate for every nonzero
initial row in its small cases, and checks actual nonzero carry
support against (12). These tests verify finite transcriptions.
Equations (1)--(12) provide the all-size arguments.
The Python checker uses ordinary lists and is not a measurement or
implementation of the packed bit-space machine from Section 3.
