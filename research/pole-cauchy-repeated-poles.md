# Higher external poles with explicit precision costs

## Constructive result

The [external-pole evaluator](pole-cauchy-rational-poles.md) extends to

\[
 C_{r,s}(c,d)=\sum_{a\in A,b\in B}
       \frac{(M^{-1})_{b,a}}{(a-c)^r(b-d)^s}
 \tag{1}
\]

for positive integer orders `r,s` and integer poles outside their
respective domains. The algorithm uses exactly `r*s` simple-pole
contractions. Its bit cost is polynomial in the **numerical orders**,
the encoded domains and poles, the displacement, and the requested
absolute precision. A binary-encoded enormous order is not treated
as a constant-cost input.

All divisions by powers of two are included in the precision budget.
The algorithm does not enumerate matrix nodes. It evaluates a larger
family of separated rational weights; the unresolved joint
factor-count pole is unchanged. This does not improve the established
factoring exponent.

## 1. An exact finite-difference identity

Write `Delta_h` for the forward difference in the pole variable:
`Delta_h f(c)=f(c+h)-f(c)`. For a positive integer `r`,

\[
 \frac{\Delta_h^{r-1}(a-c)^{-1}}
      {h^{r-1}(r-1)!}
   =\frac1{\prod_{j=0}^{r-1}(a-c-jh)}.
 \tag{2}
\]

This is an exact rational identity. For example, one forward difference
gives `h/((a-c)(a-c-h))`. Inductively, the two adjacent products share
all interior factors, and their subtraction supplies the next
integer multiple of `h`. Equivalently, partial fractions of the
right-hand side give the binomial coefficients in the left-hand side.

The right-hand side approaches `(a-c)^(-r)` when `h` tends to zero
2-adically. A finite choice of `h` suffices, with a uniform error bound
over the whole implicit domain.

## 2. Choosing the finite step

Use the matrix and domain contract of the simple-pole evaluator. In
particular, the positive odd node sets `A,B` have equal cardinality,
each is an explicit union of disjoint step-two progressions, and

\[
 M_{a,b}=\frac{2^kw}{N+2^kw-ab},\qquad w^2+w+1=0.
\]

Let `L` be the effective certificate for integrality of `2^L M^-1`.
The domain-precision result bounds `L` polynomially in the explicit
piece counts and `k`. Compute

\[
 u=\max_{a\in A}v_2(a-c),\qquad
 v=\max_{b\in B}v_2(b-d).
 \tag{3}
\]

The simple-pole algorithm computes these maxima without listing the
nodes. Since the poles are external integers, both maxima are finite
nonnegative integers bounded by the input lengths.

For requested absolute precision `p>=1`, choose

\[
 S=p+L+ru+sv+\max(u,v),\qquad h=2^S.
 \tag{4}
\]

Then `S>u,v`. Every shifted pole `c+jh` stays outside `A` and every
shifted pole `d+ell*h` stays outside `B`. Their denominator valuations
remain exactly those of the unshifted poles. Indeed, subtracting a
multiple of `2^S` cannot change a finite valuation strictly below `S`.

For `x=a-c`, factor the right-hand side of (2) as

\[
 x^{-r}\prod_{j=0}^{r-1}(1-jh/x)^{-1}.
\]

Each factor in the product is congruent to one modulo
`2^(S-v2(x))`. Therefore

\[
 v_2\left(
  \frac1{\prod_{j=0}^{r-1}(a-c-jh)}-\frac1{(a-c)^r}
 \right)\ge S-(r+1)u.
 \tag{5}
\]

The analogous column error is at least `S-(s+1)v`. The row and column
product weights themselves have valuations at least `-ru` and `-sv`.
Subtracting their two-variable product from the desired product thus
has valuation at least

\[
 S-ru-sv-\max(u,v).
\]

Multiplication by an inverse entry loses at most `L` bits. Summing
over any number of entries does not reduce the minimum 2-adic
valuation. Equation (4) consequently makes the full contraction
error at least `p`.

## 3. Exact division and its precision budget

Let `F(c,d)=C_(1,1)(c,d)` be the simple-pole contraction. Equations
(2)--(5) give the computable approximation

\[
 C_{r,s}(c,d)=
 \frac{\displaystyle
   \sum_{j=0}^{r-1}\sum_{\ell=0}^{s-1}
       (-1)^{r+s-2-j-\ell}
       {r-1\choose j}{s-1\choose\ell}
       F(c+jh,d+\ell h)}
      {h^{r+s-2}(r-1)!(s-1)!}
   +O(2^p).
 \tag{6}
\]

Put

\[
 V=(r+s-2)S+v_2((r-1)!)+v_2((s-1)!).
 \tag{7}
\]

Each call to the simple-pole evaluator requests **absolute precision
`p+V`**. The numerator in (6) is an integer linear combination, so
its approximation error remains divisible by `2^(p+V)`. Dividing
by the denominator then leaves absolute precision `p`.

The implementation keeps the denominator's power of two as a scale
exponent. Its remaining odd factor is inverted in the residue ring.
For each simple-pole return `F_j=2^(e_j)R_j+O(2^(p+V))`, choose
`E=min_j e_j`. Shift each residue by `e_j-E`, add with its signed
binomial coefficient, and invert the odd denominator modulo
`2^(p+V-E)`. The final scale is `E-V`. Thus

\[
 (p+V-E)+(E-V)=p,
\]

which explicitly accounts for all precision lost in (6). No inverse
of two is treated as a unit operation.

## 4. Bit cost and scope

The numerical value of `S` is polynomial in the numerical orders,
encoded input lengths, and `p`. The bit lengths of the shifted poles
are at most their original lengths plus `O(S+log(r+s))`. The integer
`V` is also polynomial in those parameters. Computing factorials and
binomial coefficients has polynomial cost in the numerical orders;
their binary lengths are `O((r+s)log(r+s))`.

Each of the `r*s` simple-pole calls has polynomial bit cost at this
enlarged precision and pole length. Their exact affine valuation
counts remain independent of the numerical matrix cardinality.
Combining the returned residues likewise has polynomial cost. This
proves the stated algorithmic bound.

Finite linear combinations of these contractions can be evaluated by
linearity, with the displayed number of terms charged and any
rational coefficient denominator included in the absolute precision
budget. This statement does not grant free evaluation of a short
circuit that expands into exponentially many distinct pole terms.

## 5. Independent exact checks

The [checker](../tools/check_pole_cauchy_repeated_poles.py) verifies
280 exact univariate identities through order seven. It also compares
36 complete contractions against dense rational inversion, using
orders one through three, zero and negative poles, progression
unions, and external poles whose denominator valuations reach 20.

For the 32 controls with a nontrivial finite difference, it separately
forms the exact product approximation in (2) and verifies its
valuation error. Four invalid order or precision contracts are
rejected. The [evidence](../evidence/pole-cauchy-repeated-poles.json)
records the precision parameters and source hashes.

```sh
python3 tools/check_pole_cauchy_repeated_poles.py \
  --output evidence/pole-cauchy-repeated-poles.json
```

The [implementation](../tools/pole_cauchy_repeated_poles.py) is separate
from the simple-pole evaluator and preserves that evaluator's input
contract and precision certificates.
