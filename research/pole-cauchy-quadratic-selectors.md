# Quadratic poles give inverse residue-block sums

## Constructive result

For a binary-encoded integer `h>=0` and residue `r`, the rational weight

\[
 \boxed{\displaystyle
 \chi_{r,h}(a)=
   \left(\frac{2^h}{a-r+2^hw}\right)^3,
 \qquad w^2+w+1=0,}
 \tag{1}
\]

is integral at every integer `a` and reduces modulo two to the indicator
of `a=r mod 2^h`. Its representation and evaluation precision have cost
polynomial in the **numerical value of `h`**, which is the bit length of
the displayed coefficient `2^h`. This does not assign constant cost to
a binary-encoded exponential shift.

The exact Cauchy inverse contraction against a row weight and a column
weight of this form is computable in polynomial bit time on explicit
unions of binary-encoded odd progressions. Consequently, if `L` is the
effective inverse-integrality certificate, the algorithm computes

\[
 \sum_{\substack{a\in A,\ a=r\ (2^h)\\
                  b\in B,\ b=s\ (2^g)}}
     (2^L M^{-1})_{b,a}\pmod2
 \tag{2}
\]

without enumerating matrix nodes. Arbitrary numerical matrix
cardinality is allowed.

This is a new family of nonregular rational weights. Its poles lie in
the quadratic extension rather than at integer nodes. The
[regular selector precision lower bound](regular-analytic-selector-precision.md)
retains its unit-inverse gate restriction and is not contradicted.

The output in (2) concerns the **normalized inverse**. It does not count
the true factor entries of the original matrix. Section 7 identifies
which Smith directions its parity can retain. The established factoring
score is unchanged.

## 1. Exact residue selectors

Write `a-r=2^j u` with `u` odd. If `j<h`, the denominator in (1) has
valuation `j`, so the fraction has positive valuation. Its cube reduces
to zero.

If `a-r` is divisible by `2^h`, divide the denominator by `2^h`. Its
residue in the field with four elements is either `w` or `1+w`.
Both are nonzero, and every nonzero element of that field has cube one.
Thus (1) reduces to one exactly on the specified residue class. This
also covers `a=r`, and it covers every integer `a`, including negative
ones. No domain-specific maximum denominator is needed for the selector
identity itself:

\[
 v_2(a-r+2^hw)=\min(v_2(a-r),h)\le h.
 \tag{3}
\]

In particular, `2^h/(a-r+2^hw)` is integral everywhere.

For the computation below, use the canonical representative
`0<=r<2^h`. The associated pole is `c=r-2^h w`, which is outside every
integer node set. This choice also makes all cross-product factors and
the resolvent denominator nonzero; Section 5 checks that explicitly.

## 2. The rational resolvent extends to quadratic poles

Retain the matrix contract from
[external integer poles](pole-cauchy-rational-poles.md):

\[
 t=2^kw,\quad z=N+t,\quad M_{a,b}=\frac{t}{z-ab},
 \quad m=|A|=|B|,
 \tag{4}
\]

where `N` is positive and odd, the nodes are distinct positive odd
integers, and the displacement separates the nonzero integer
differences `N-ab`. Each domain is an explicit union of disjoint
step-two progressions.

The rational identity

\[
 F(c,d)=\sum_{a,b}\frac{(M^{-1})_{b,a}}{(a-c)(b-d)}
       =\frac z{t(cd-z)}(1-Q(c,d)),
 \tag{5}
\]

with

\[
 Q(c,d)=
  \frac{\prod_a(z-ad)\prod_b(z-bc)}
       {z^m\prod_a(a-c)\prod_b(b-d)},
 \tag{6}
\]

is algebraic over the quadratic field. It therefore applies to poles
`c,d` with integer coordinates in the basis `1,w`, provided the displayed
denominators are nonzero. The algorithm in this note additionally
requires every numerator factor at the base poles to be nonzero. That
extra condition permits unit normalization of each product jet. The
selector poles in Section 1 always satisfy it.

Equation (6) separates into a function of `c` times a function of `d`:

\[
 Q(c,d)=U(c)V(d),\qquad
 U(c)=\frac{\prod_b(z-bc)}{z^m\prod_a(a-c)},\quad
 V(d)=\frac{\prod_a(z-ad)}{\prod_b(b-d)}.
 \tag{7}
\]

This separation permits two univariate product jets, followed by a small
two-variable recurrence. There is no need to expand a full matrix or
construct a dense two-variable product of degree `m`.

## 3. Affine quadratic products have short valuation partitions

Consider a sequence

\[
 q_i=(A+Ci)+(B+Di)w,\qquad 0\le i<\ell.
 \tag{8}
\]

In the unramified quadratic extension,

\[
 v_2(q_i)=\min(v_2(A+Ci),v_2(B+Di)).
\]

Assume none of the factors is zero. Remove the common power of two from
the four affine coefficients. If both step coordinates are even, one
constant coordinate is odd, and the whole sequence is a unit progression
with an even step.

Otherwise split the index by parity. The step is nonzero modulo two in
the quadratic residue field. At most one parity makes both coordinates
of the new constant even. Every other parity branch is a unit progression
with an even step. Divide the possible even branch by two and continue
only that branch.

At each continuing split the number of indices at least halves. A
singleton is handled directly, including an isolated high valuation.
Thus only `O(log(1+ell))` strata are required. Each is an index
progression together with a common factor valuation. Summing
`count * valuation` gives the exact product valuation as a binary integer.
No factor list is generated.

The normalized unit products are evaluated with the existing
[modular interval product algorithm](modular-interval-jet-arithmetic.md).
The same partition works for short jets whose nonconstant coefficients
are divisible by a sufficiently large power of two. If `S` exceeds
every base-factor valuation, perturbing the affine sequence by
`2^S X` times another integral affine sequence leaves every normalized
jet integral. Its constant coefficient is a unit, and its full step jet
is even. The interval product algorithm therefore applies unchanged.

The [implementation](../tools/pole_cauchy_quadratic_selectors.py) records
the index progression of each stratum. It normalizes both the base
factor and the perturbation before entering the residue ring. Negative
coordinates and constant affine sequences are included. A zero factor
is rejected explicitly.

## 4. Exact jets evaluate higher quadratic poles

For positive numerical orders `r,s`, define

\[
 C_{r,s}(c,d)=\sum_{a,b}
       \frac{(M^{-1})_{b,a}}{(a-c)^r(b-d)^s}.
 \tag{9}
\]

Choose `S` larger than the valuations of all factors in the four products
of (6) and larger than `v2(cd-z)`. Set `delta=2^S`. The exact formal
identity

\[
 [X^{r-1}Y^{s-1}]F(c+\delta X,d+\delta Y)
       =\delta^{r+s-2}C_{r,s}(c,d)
 \tag{10}
\]

follows by expanding each separated denominator. This is coefficient
extraction at an explicitly normalized unit point. There is no
finite-difference approximation or factorial division in this version.

The affine strata compute univariate jets of `U(c+delta X)` through
degree `r-1` and `V(d+delta Y)` through degree `s-1`. Keep their possibly
large product valuations as separate binary integers. Their sum is

\[
 q=v_2(Q(c,d)).
\]

Only the normalized unit jets are multiplied. Large opposite valuations
are cancelled exactly before a residue precision is chosen.

Write `v0=v2(cd-z)`. After division by `2^v0`, the remaining denominator
has the shape

\[
 E(X,Y)=\alpha+\beta X+\gamma Y+\eta XY,
 \tag{11}
\]

where

\[
 \alpha=(cd-z)/2^{v_0},\quad
 \beta=\delta d/2^{v_0},\quad
 \gamma=\delta c/2^{v_0},\quad
 \eta=\delta^2/2^{v_0}.
\]

The constant `alpha` is a unit and all coefficients are integral. For
the coefficients `e_(i,j)` of `1/E`, start with `e_(0,0)=alpha^(-1)`
and use

\[
 e_{i,j}=-\alpha^{-1}
       (\beta e_{i-1,j}+\gamma e_{i,j-1}+\eta e_{i-1,j-1}),
 \tag{12}
\]

with missing indices interpreted as zero. This requires `O(rs)` ring
operations. Convolving this table with the product of the two
univariate unit jets at the single desired coefficient also takes
`O(rs)` operations.

### Precision accounting

Let `L` certify integrality of `2^L M^-1`, and let

\[
 u=\max_a v_2(a-c),\qquad v=\max_b v_2(b-d).
\]

These maxima are part of the affine-stratum computation. Applying the
entrywise inverse bound to the simple-pole contraction (5) gives

\[
 q\ge\min(0,k+v_0-L-u-v).
 \tag{13}
\]

The argument is the same as for integer poles: each summand in `F` has
valuation at least `-L-u-v`, and the prefactor in (5) has valuation
`-k-v0`. The quadratic pole coordinates do not change this estimate.

For requested absolute precision `p`, put

\[
 J=(r+s-2)S,\qquad T=p+J+k+v_0,\qquad e=\min(0,q).
 \tag{14}
\]

Compute the normalized expression modulo `2^(T-e)` and return it with
scale exponent

\[
 e-k-v_0-J.
 \tag{15}
\]

The sum of the residue precision and this scale is exactly `p`. The
`J` term explicitly pays for coefficient division in (10).

If `q<T`, each normalized product jet is needed only modulo `2^(T-q)`.
Multiplication by `2^(q-e)` restores the required precision `T-e`,
coefficient by coefficient. If `q>=T`, the entire normalized product
term vanishes to the required precision and its unit jets are omitted.
All these jets are integral by the choice of `S`, so the omission also
applies to their nonconstant coefficients.

Equation (13), the polynomial inverse-loss bound for explicit progression
unions, and the encoded bound on individual affine coordinates ensure
that every numerical precision is polynomial in the input lengths,
`k`, `p`, and the numerical orders. The number of affine strata is
polynomial in the binary interval lengths. The Newton divisions inside
the interval-product routine retain their previously proved precision
allowance. This proves polynomial bit cost for (9), independently of
the numerical dimension `m`.

## 5. Applying the cubic selectors

Use

\[
 c=r-2^h w,\qquad d=s-2^g w,
 \qquad 0\le r<2^h,\quad0\le s<2^g.
\]

Both poles have nonzero imaginary coordinate, so neither can equal an
integer matrix node. For positive nodes, the imaginary coordinates of
`z-ad` and `z-bc` are `2^k+a*2^g` and `2^k+b*2^h`, respectively, and
are positive. Those numerator factors cannot vanish. The imaginary
coordinate of `cd-z` is

\[
 -r2^g-s2^h-2^{h+g}-2^k<0,
\]

so this denominator cannot vanish either. Thus every selector request
satisfies the nonvanishing contract.

Evaluate `C_(3,3)(c,d)` using Section 4, and multiply the result by

\[
 2^{L+3h+3g}.
\]

The resulting sum is integral by (1) and the definition of `L`.
Reducing modulo two gives exactly (2). Its integral output is extracted
from the returned scale and residue; the implementation checks the
required divisibility when the intermediate scale is negative.

The whole-domain denominator maxima satisfy `u<=h` and `v<=g` by
(3). Numerical work therefore grows polynomially with the selected
residue precisions. There is no factor proportional to `2^h` or `2^g`
in the number of coefficient operations.

## 6. Why the regular-selector lower bound does not apply

At an odd residue center, the inverse gate in (1) has an even
denominator on part of the odd input domain. It is not a regular
unit-inverse gate of the earlier compiler. Its powers of two are
handled through exact affine valuations and the rational identity,
then restored by an explicitly charged scale.

The regular representation required an exponentially large external
dyadic denominator for a deep residue selector. Here a cubic rational
denominator in the quadratic field suffices, and the product identity
evaluates its contraction by a different algorithm. The distinction is
a precise extension of the accepted representation, not a claim that
the earlier lower bound was false.

## 7. What inverse block parity reveals

Let the integral matrix `M` have Smith orders
`sigma_1<=...<=sigma_m`. If `2^L M^-1` is integral, then
`L>=sigma_m`. Modulo two, the rank of the normalized inverse is

\[
 \#\{i:\sigma_i=L\}.
 \tag{16}
\]

This follows by inverting a Smith decomposition and multiplying by
`2^L`. If `L>sigma_m`, the normalized inverse is entirely zero modulo
two. If `L>0`, the original unit Smith directions `sigma_i=0` vanish
in that reduction. Thus block parity cannot be identified with the
original factor-pair incidence matrix.

A concrete control uses `N=15`, `A=B={1,3,5}`, and `k=5`. Its exact
Smith orders are `(0,0,4)`, and the inverse certificate gives `L=4`.
The singleton inverse block at `(a,b)=(1,1)`, which is not a factor
pair, has parity `1+w`. The singleton blocks at the true factor pairs
`(3,5)` and `(5,3)` have parity zero after this normalization.

The algorithm provides access to structured inverse block sums.
Recovering the original unit directions or aggregating the product
congruence graph `ab=N mod 2^k` still requires an additional argument.
No such argument or improved factoring exponent is asserted here.

## 8. Exact controls and an implicit example

The [checker](../tools/check_pole_cauchy_quadratic_selectors.py) verifies
affine valuation partitions and normalized product jets against literal
products. It checks selector parity at positive and negative integers,
compares repeated quadratic-pole contractions with dense exact
inversion, and compares residue-block outputs with literal sums of
normalized inverse entries. Integer simple poles are cross-checked
against the independent evaluator. Invalid contracts are rejected.

The large control uses a 402-bit input, dimension `2^200+13`, and
residue classes `5 mod 2^7` and `11 mod 2^9`. It computes the selected
block parity without listing the nodes, and checks the full underlying
cubic-pole contraction against a higher-precision replay. The
[evidence](../evidence/pole-cauchy-quadratic-selectors.json) records its
precision certificate and source hashes.

```sh
python3 tools/check_pole_cauchy_quadratic_selectors.py \
  --output evidence/pole-cauchy-quadratic-selectors.json
```

Use `--small-only` for the independent finite controls alone.
