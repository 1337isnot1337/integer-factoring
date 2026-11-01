# Base-prime digits, exact carries, and the cost of lifting Fourier queries

## Scope and result

This note tests a different arithmetic starting point for the half-rate
factoring evaluator. Base-`p` digits range over the whole field `F_p`.
That makes a quadratic Fourier sum easy to eliminate, but an integer
carry has information that its residue modulo `p` does not retain.

The exact conclusions are:

* Reducing the carry equations modulo `p` retains only the lowest-digit
  product condition. For an input not divisible by `p`, the relaxed
  system admits exactly `(p-1)p^(2m-2)` pairs of `m`-digit words.
* A polynomial-size prime-power lift makes all carries exact. A matrix
  of rank `r` modulo `p` then leaves at most `p^(m-r)` canonical digit
  candidates in a mixed affine block, regardless of the lift precision.
* The digit indicator at that precision has `p^ell-p+1` nonzero Fourier
  coefficients. More strongly, the actual carry tensor has at least
  `p^((ell-1)(m-1))` nonzero frequencies whose convolution Hankel matrix
  is invertible modulo `p`.
* The field-characteristic restriction in the determinant lemma can be
  removed using an integer trace identity and a short `p`-adic modulus.

The last two statements separate inexpensive individual block queries
from the cost of summing their boundary frequencies. Explicitly visiting
all surviving lift frequencies is still exponential in the input bit
length. There is no improved unconditional factoring exponent here.
No Lean formalization is added.

For the best local rank-query bound, use the
[division-free scalar Berkowitz construction](small-space-berkowitz-rank.md):
it gives quasipolynomial time in the matrix order and linear auxiliary
space in every characteristic. Section 5 below retains an independent
integer-trace implementation and its explicit precision accounting.

The arithmetic and carry conventions come from the
[Fourier network on PR #23](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/logarithmic-wire-network.md),
pinned at commit `147652574e7ce1ff68b2c378712846cbf438da34`.
The related mixed-block identity and prime-field resource accounting
are in [the digit-box note](hankel-digit-box-contraction.md).

## 1. What reduction modulo the digit base discards

Let `p` be prime and write two words with `m` base-`p` digits as

```math
a=\sum_{i=0}^{m-1}a_i p^i,\qquad
b=\sum_{i=0}^{m-1}b_i p^i,\qquad 0\le a_i,b_i<p.
```

Their convolution coefficients satisfy

```math
c_j=\sum_{i+k=j}a_i b_k,\qquad
0\le c_j\le C:=m(p-1)^2.
```

For an input `0<N<p^m`, extend its digits by zeros to convolution
length `L=2m-1`. The exact carry equations are

```math
c_j+t_j=N_j+p t_{j+1},\qquad t_0=t_L=0.                \tag{1}
```

Carries in a valid multiplication have

```math
0\le t_j\le T:=m(p-1),                                \tag{2}
```

because `(C+T)/p=T`. The exact equations telescope to `ab=N`.

Modulo `p`, every outgoing-carry term in (1) vanishes. The first
equation is `a_0b_0=N_0`. For each `j>=1`, independently choose

```math
t_j=(N_j-c_j)\bmod p\quad\hbox{in }[0,p).
```

These are within the permitted range (2). All remaining equations hold,
and the final carry does not occur with a nonzero coefficient. Therefore
the relaxed system has exactly the digit pairs whose lowest digits
satisfy the first equation.

If `p` does not divide `N`, there are `p-1` choices for the pair of
lowest digits, while the other `2m-2` digits are unrestricted:

```math
\#\text{relaxed pairs}=(p-1)p^{2m-2}.                 \tag{3}
```

For example, at `p=2,m=2,N=3`, the false pair `a=b=3` has
convolution `(1,2,1)` and relaxed carries `(0,1,1,0)`. Every equation
holds modulo two, while `ab=9`. The exact first carry should have been
zero. Counting the relaxed witnesses cannot distinguish factorization
from these spurious products.

Passing to an extension field `F_(p^d)` does not restore the missing
information: its characteristic is still `p`, so the same terms vanish.

### An arbitrary completion cannot recover information already identified

This loss is not specific to the relaxed carry tensor. No function of
the convolution vector reduced modulo `p` alone can equal the exact
integer-product predicate on all realizable digit inputs.

For `p>=3`, take `N=(p-1)^2`. The constant-digit pairs
`(p-1,p-1)` and `(1,1)` have the same reduced convolution, the constant
one. Their integer products differ, and both pairs are within `[1,N]`.
For `p=2`, the pairs `(1,5)` and `(3,3)` both have reduced polynomial
product `1+X^2`, but their integer products are five and nine.

The phenomenon is not confined to perfect powers. At `p=3`, the true
pair `(7,19)` for `N=133` and the false pair `(5,11)` have the same
reduced convolution. At `p=5`, the true pair `(11,13)` for `N=143`
and the false pair `(19,17)` have the same reduced convolution. All
four factors in either example are positive, odd, and at most `N`.
Thus a prospective completion must retain additional lifted data;
altering only its values away from the realizable reduced convolutions
cannot repair these collisions.

## 2. An exact prime-power lift and a useful rank bound

Choose the least `ell` with

```math
P=p^\ell>(mp+1)(p-1).                                \tag{4}
```

For digits and carries in the stated ranges, each integer discrepancy
in (1) has absolute value below `P`, because

```math
-(pT+p-1)\le c_j+t_j-N_j-p t_{j+1}\le C+T=pT.
```

Thus its congruence modulo `P` is equivalent to the exact equality.
Convolution itself has no wrap because `P>C`. The lift is effective
by successive integer powers. For polynomial-size `p` and `m`, `P`
is polynomial-size as a number, and its bit length is `O(log(mp))`.

The full-field digit domain is now gone: each digit is the subset
`{0,...,p-1}` of `Z/PZ`, not that whole ring. Nevertheless an exact
affine block can exploit the full base-field domain.

**Digit-lift lemma.** Given `D x=b (mod P)` with `s` canonical digits
`0<=x_i<p`, let `r` be the rank of `D mod p`. Choose a nonsingular
rank-`r` minor modulo `p`. Enumerate the `s-r` nonpivot digits,
solve for the pivot residues modulo `p`, and take their unique
representatives in `[0,p)`. Check the original equations modulo `P`.
This enumerates each exact solution once, using at most

```math
p^{s-r}                                               \tag{5}
```

candidates. The proof is immediate from the invertible pivot system:
every solution has one enumerated free-digit assignment, and its pivot
digits are already determined by their residues. In particular full
column rank modulo `p` leaves at most one candidate at any precision.
Higher-precision consistency is checked, not assumed from mod-`p`
consistency.

This gives a real mixed-block query method. It does not say that both
digit families can be summed by unrestricted-ring Gauss elimination.

## 3. Exact Fourier support of one lifted digit

Use a fixed coefficient field of characteristic different from `p`.
Let `R` be its quotient by

```math
\Phi_P(Z)=1+Z^{P/p}+\cdots+Z^{(p-1)P/p},\qquad\zeta=Z.
```

It can be a product of fields. The polynomial `Z^P-1` is separable,
and the roots of `Phi_P` all have exact order `P`. Consequently
`1-zeta^c` is a unit whenever `P` does not divide `c`. The normalized
Fourier matrix on `Z/PZ` is invertible, since its orthogonality scalar
`P` is also a unit.

The Fourier transform of the canonical digit indicator is

```math
U(c)=\sum_{d=0}^{p-1}\zeta^{cd}.                       \tag{6}
```

For `c=0`, it is the nonzero unit `p`. For `c!=0 (mod P)`,

```math
U(c)=\frac{1-\zeta^{pc}}{1-\zeta^c}.                 \tag{7}
```

The denominator is a unit. The numerator is zero exactly when
`P | pc`; otherwise it is a unit as well. Thus (6) has precisely
the `p-1` zeros at the nonzero multiples of `P/p`, and

```math
|\operatorname{supp}U|=P-p+1.                         \tag{8}
```

For `m` independent digit coordinates the support size is
`(P-p+1)^m`. When `ell=1`, (8) is one, recovering ordinary full-field
orthogonality. Every exact lift in (4) has `ell>=2`; the support is then
dense. For example, at `p=2,P=4`, the digit sum at frequency one is
`1+zeta`, not zero, whereas the sum over the whole ring is zero.

Formula (8) concerns an explicit expansion into all nonzero frequencies.
The digit transform still has a short geometric description. A different
algorithm could potentially sum many such descriptions together.

## 4. Carry-support pruning still leaves many invertible Hankel terms

Let `F_N(c)` be the full legal-carry indicator on `L=2m-1` coefficient
coordinates represented in `[0,P)`, with zero initial and final carry.
One may either enforce the natural carry range (2), or the larger range
`[0,P)`; the argument works with either convention. The exact equalities
telescope to

```math
\sum_{i=0}^{2m-2}c_i p^i=N<p^m.                       \tag{9}
```

Nonnegative coefficients imply `c_i=0` for every `i>=m`.
Therefore the carry Fourier value is independent of the upper
`m-1` frequencies.

Here the carry gate tests the literal integer equation, including the
chosen carry range. Replacing that test by a congruence modulo `P`
while allowing all carries in `[0,P)` is a different, unsound shortcut:
the no-wrap proof in Section 2 used the natural bound `T`. Its range
cannot be dropped when invoking that proof.

The induced function of the first `m` coordinates is not invariant
under adding `h=P/p` to coordinate `m-1` modulo `P`. Indeed its
ordinary base-`p` digit vector is legal; after this change it is
illegal because `h>=p` makes the left side of (9) greater than `N`.
There is no coordinate wrap in this example: the original highest
digit is at most `p-1`, and `h+p-1<P` for `ell>=2`.

Fourier inversion now gives a nonzero lower-frequency value with

```math
\gamma_{m-1}\not\equiv0\pmod p.                      \tag{10}
```

If every nonzero Fourier value had that coordinate divisible by `p`,
translation by `h` would leave the function unchanged, contrary to
the preceding paragraph. This implication uses invertibility of the
Fourier matrix and remains valid in the quotient root algebra.

Fix the first `m` coordinates of a frequency satisfying (10).
Allow each upper coordinate to be an arbitrary multiple of `p`
modulo `P`. There are

```math
(P/p)^{m-1}                                          \tag{11}
```

such extensions, all with the same nonzero carry weight. Their Hankel
matrix `H_ij=gamma_(i+j)` is anti-triangular modulo `p`, with the
unit `gamma_(m-1)` on its anti-diagonal. It is therefore invertible
modulo `p`, and hence over `Z/PZ`.

In logarithmic-digit parameters `m log p=Theta(n)` with `m>=2`,
(4) implies that (11) is `exp(Omega(n))`. For fixed `p`,
`m=Theta(n)` and the count can be even larger. Thus individually
processing only the nonzero carry frequencies, even after discarding
all noninvertible Hankel cases, still has exponential work.

This is a limitation of that explicit traversal. It is not a lower
bound for collective Fourier summation or for factoring algorithms.

### Factor-length splits and an algebraic support bound

Splitting the factor lengths in advance can remove the unused upper
coefficient positions. Thus (11) is a statement about the two full
`m`-digit words, not every possible Fourier layout. A separate bound
applies to the natural bounded-carry tensor after such a split.

Consider `L` coefficient coordinates in `Z/PZ`, a fixed incoming and
outgoing carry, and the exact carry equations with every intermediate
carry in `[0,T]`. Suppose the carry tensor `f` is nonzero. Its legal
coefficient vectors are determined by its `L-1` internal carries, since
`c_j=N_j+p t_(j+1)-t_j`. Consequently

```math
|\operatorname{supp}f|\le(T+1)^{L-1}.
```

For every nonzero function on a finite abelian group `G`, over a
splitting field whose characteristic does not divide `|G|`, the support
product inequality is

```math
|\operatorname{supp}f|\,|\operatorname{supp}\widehat f|
\ge |G|.
```

Here is an algebraic proof, so no complex absolute-value argument is
assumed. The translation matrix `M_(x,y)=f(y-x)` is diagonalized by
the invertible Fourier matrix, hence its rank is
`|supp fhat|`. Choose a row basis. The union of its supports covers
every column: otherwise that column is zero on all rows, contrary to
translation of one nonzero value of `f`. Each basis row has `|supp f|`
nonzero entries. Their union therefore has size at most
`rank(M)|supp f|`, proving the inequality.

Specialize the cyclotomic quotient to a splitting field to apply the
inequality; specialization cannot create additional nonzero Fourier
coefficients. For the bounded-carry function this gives

```math
|\operatorname{supp}\widehat f|
\ge\frac{P^L}{(T+1)^{L-1}}
\ge P p^{L-1},
```

where the last inequality uses the integer consequence
`P>=p(T+1)` of (4). This is `exp(Omega(n))` when
`L=Theta(n/log p)`, including the usual factor-length splits. A fixed
nonzero outgoing carry permits a product with one additional base-`p`
digit and does not change the count of internal carry choices.

This support argument concerns the natural bounded-carry tensor. It
does not cover every possible extension of the integer-product predicate
outside the range of actual convolution coefficients, and it does not
forbid a collective sum over its support.

## 5. Small characteristic is not a determinant-storage obstruction

The earlier composition/trace determinant formula divided by
`1,...,m`, so it needed characteristic greater than the matrix order.
That requirement can be removed for rank computations modulo any prime
`p` without retaining a dense matrix.

Lift the entries of an `m by m` matrix to nonnegative integers below
`p`, and write `t_k=tr(U^k)` over the integers. For a partition
`lambda=(1^(a_1)2^(a_2)...)` of `m`, put

```math
z_\lambda=\prod_k k^{a_k}a_k!,\qquad
\ell(\lambda)=\sum_k a_k.
```

The integer trace identity is

```math
m!\det U=
\sum_{\lambda\vdash m}(-1)^{m-\ell(\lambda)}
\frac{m!}{z_\lambda}\prod_k t_k^{a_k}.                 \tag{12}
```

It follows by grouping equal ordered-composition terms in the Newton
identity. Every coefficient `m!/z_lambda` is an integer: it counts
permutations with cycle lengths prescribed by `lambda`.

Let `v=v_p(m!)` and work modulo `Q=p^(v+1)`. Its bit length is

```math
(v+1)\log_2p\le m\frac{\log_2p}{p-1}+\log_2p
\le m+\log_2p.                                      \tag{13}
```

Strip factors of `p` from each factorial and each `k`. If `w=v_p(z_lambda)`,
then `w<=v`, and the coefficient in (12) can be formed modulo `Q` as

```math
p^{v-w}\,(m!)_{unit}\,(z_\lambda)_{unit}^{-1}\pmod Q.
```

Only units are inverted. If `R` is the accumulated residue on the
right of (12), then it is divisible by `p^v` as a canonical integer
in `[0,Q)`. The desired determinant modulo `p` is

```math
(R/p^v)\,((m!)_{unit}\bmod p)^{-1}\pmod p.             \tag{14}
```

No invalid division inside `F_p` has occurred.

### Avoid a logarithmic stack of long trace residues

It would be incorrect to retain an `O(m)`-bit modular accumulator at
every matrix-power recursion level and still claim linear workspace.
Instead compute bits of the nonnegative integer traces before their
reduction modulo `Q`.

Every `U^k` entry for `k<=m` has at most
`O(m log(mp+2))` bits. To query one bit of a product of two virtual
integers, run ordinary binary multiplication column by column up to
that position. The carry is bounded by the operand bit length, so its
storage is logarithmic in that length. A bit of a sum of `m` virtual
integers is obtained in the same way, with a carry bounded by `m`.

Apply these bit-query routines to balanced matrix powering. The power
parameter halves on each recursive level. Only indices, bit positions,
and small column carries are retained, using
`O(log(m+2) log(m log(p+2)+2))` bits. The bit-query time is
`(m+log(p+2))^{O(log(m+2))}`. Reading the trace bits from high to low,
one outer `O(m+log p)`-bit accumulator reduces them modulo `Q`.

Enumerate partitions using a separator mask among all `2^(m-1)`
compositions, retaining only nondecreasing part lists and rescanning
the mask when multiplicities are needed. The mask has `O(m)` bits;
no `m`-entry array of `O(log m)`-bit parts is necessary.
Equations (12)--(14) therefore give

```math
T=2^{O(m)}(m+\log(p+2))^{O(\log(m+2))}(1+T_{entry}),
S=S_{entry}+O(m+\log p+\log(m+2)\log(m\log(p+2)+2)).   \tag{15}
```

For polynomial-size `p` in the intended parameter ranges, this is
`2^O(m)` time and `O(m)` additional working bits. Enumerating row and
column masks and applying (15) gives rank and pivot minors in the same
scale. This is a resource-accounting lemma for standard exact linear
algebra, not a claim that determinant computation itself is new.

## 6. Why separate CRT counts do not enforce a common product

CRT reconstructs a residue of one fixed integer from its residues. It
does not reconstruct the size of an intersection of different witness
sets from their individual cardinalities.

For arithmetic moduli `M_i`, the statements `ab=N (mod M_i)` must hold
for the same pair `(a,b)`. Evaluating a separate count for each modulus
and applying CRT to those counts computes neither the intersection nor
the integer factor count. A correct conjunction retains the joint
residue conditions, equivalently the product modulus when the moduli
are coprime.

If factor values are represented by unrestricted CRT residue vectors,
their congruence product factors by coordinate, but their common integer
range restriction couples those coordinates. Without a range or a
no-wrap argument, solutions include `ab=N+k product(M_i)`.
For example `3*4=5 (mod 7)` is a false integer-factor pair even though
both candidates are at most five. This is different from the existing
CRT theorem for reconstructing an already exact divisor count modulo
small coefficient moduli; that theorem preserves the same integer
witness predicate in every query.

## Verification and unresolved step

The exact checker tests the relaxed carry count, false-witness examples,
the no-wrap bound, the digit-lift lemma, the dense Fourier support and
unit property, nonzero carry frequencies with invertible Hankel lifts,
and the small-characteristic integer trace formula. Its Python lists
are finite reference representations, not a measured implementation of
the packed workspace bound.

The unresolved constructive step is a collective sum of these lifted
digit and carry tensors that avoids visiting their exponentially many
nonzero frequencies. Small characteristic and CRT alone do not supply
that sum.
