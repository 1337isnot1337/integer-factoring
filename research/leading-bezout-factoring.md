# One global leading Bezout coefficient suffices for factoring

## Reduction and scope

There is a deterministic polynomial-query reduction from smallest-prime
factoring to **one leading coefficient** of a globally regularized
Bezout inverse. The requested coefficient precision and displacement
exponent are linear in the original input bit length. The interval
length and clearing exponent are large integers with short binary
descriptions. No factor-count bound is required.

The new ingredient is extra interval padding. On every actual factor
row, the normalized Vandermonde derivatives become the same unit
modulo the requested precision. The leading coefficient then becomes
a common unit times a positive sum of proper divisors. Its residue
cannot cancel to zero unless the selected factor support is empty.

This is a more specific extraction interface, with an implemented
conditional wrapper. The global coefficient evaluator remains open.
Consequently the complete proved factoring score is still `1/2`.

## 1. The one-coefficient contract

Let `N>=3` be odd, with bit length `n`, and let `a_0` be an odd threshold
with `3<=a_0<=N`. Set

\[
 p=2n+1,\qquad u=n+p=3n+1,\qquad m=2^u,
 \tag{1}
\]

and use the ordinary progressions

\[
 A=\{a_0+2i:0\le i<m\},\qquad
 B=\{3+2j:0\le j<m\}.
 \tag{2}
\]

Both intervals extend beyond `N`. Their true factor rows are exactly

\[
 F_{a_0}=\{a:a\mid N,\ a\ge a_0,\ N/a\ge3\}.
 \tag{3}
\]

Choose `H` from the two product corners as in the
[fixed-scale interface](regularized-bezout-factor-moments.md), and put

\[
 K=H+p,\qquad E=2m-2-u,\qquad
 P=\prod_{a\in A}(X-a),\qquad
 Q_t=\prod_{b\in B}(N+2^K w-bX).
 \tag{4}
\]

Let `V_t` be the unique polynomial of degree less than `m` with
`V_t Q_t=1 mod P`. The requested scalar is

\[
 \boxed{c(N,a_0)=
 [X^{m-1}]\,w\,2^{K+2E}V_t(X)\pmod {2^p}.}
 \tag{5}
\]

The leading coefficient is unchanged by `X=1+Y`, so either coefficient
convention gives the same query. The theorem below makes (5) real
modulo `2^p`. No other coefficient of `V_t` is requested.

The encoded domain sizes have `O(n)` bits, `H=O(n)`, and `K=O(n)`.
The numerical clearing exponent `K+2E` is exponential in `n`, but its
binary representation has `O(n)` bits. A proposed provider must honor
this normalized-coefficient contract without expanding that power of two.

## 2. Derivative units are constant on every factor row

For the row at index `i`, write

\[
 \alpha_i=2^{-E}P'(a_0+2i),\qquad
 \alpha_0=(-1)^{m-1}\operatorname{odd}((m-1)!),
 \tag{6}
\]

where `odd(z)=z/2^(v_2(z))`. The factorial derivative formula gives

\[
 \frac{\alpha_i}{\alpha_0}
 =\frac{(-1)^i}{\binom{m-1}{i}},\qquad
 (-1)^i\binom{m-1}{i}
 =\prod_{\ell=1}^{i}\left(1-\frac m\ell\right).
 \tag{7}
\]

If `a` is a factor row, then `a<N` and its index satisfies `i<2^n`.
Every nonempty factor in (7) differs from one by a quantity of
valuation at least `u-n>=p`. Thus

\[
 \alpha_i\equiv\alpha_0\pmod {2^p}
 \quad\text{on every factor row.}
 \tag{8}
\]

The same proof applies to the matching column `b=N/a`, whose index
`j=(b-3)/2` is also less than `2^n`. Its derivative unit is congruent
to the same `alpha_0`, because the derivative formula depends on the
progression length and index, not on its first node.

Finally, every odd integer satisfies
`a^(2^u)=1 mod 2^(u+2)` for `u>=1`, by repeated squaring starting
with `a^2=1 mod8`. Hence

\[
 a^{-(m-1)}\equiv a\pmod {2^p}.
 \tag{9}
\]

## 3. A positive divisor sum prevents cancellation

The fixed regularization theorem identifies the scalar in (5) as

\[
 c(N,a_0)=\sum_{a\in F_{a_0}}
       \frac1{a^{m-1}\alpha_{i(a)}\alpha_{j(N/a)}}
       \pmod {2^p}.
 \tag{10}
\]

Apply (8)--(9). The complete coefficient satisfies

\[
 \boxed{\alpha_0^2 c(N,a_0)
       \equiv\sum_{a\in F_{a_0}}a\pmod {2^p}.}
 \tag{11}
\]

The sum on the right is a nonnegative integer less than `N^2`, hence
less than `2^(2n)<2^p`. Since `alpha_0` is odd,

\[
 \boxed{c(N,a_0)\ne0\pmod {2^p}
       \quad\Longleftrightarrow\quad F_{a_0}\ne\varnothing.}
 \tag{12}
\]

The decoder needs only the zero test. It does not have to compute
`alpha_0` or the divisor sum. If the exact sum is desired, its multiplier
is computable by the existing normalized factorial-unit algorithms.

The same argument sharpens the primitive-exponent result. If `F_(a_0)`
is nonempty, the residue polynomial `W` from the
[displaced extraction theorem](displaced-bezout-extraction.md) has

\[
 \mu=v_c(W)\le v_2([X^{m-1}]W)
 =v_2\left(\sum_{a\in F_{a_0}}a\right)\le2n-1.
 \tag{13}
\]

Here equality of the leading-coefficient valuations follows from (11)
and the strict bound below `2^p`. Therefore the standard choice `K=H+p`
already gives the exact primitive formula

\[
 r_K=K+2E-\mu>2E+H-1
 \quad\text{when }F_{a_0}\ne\varnothing,
 \tag{14}
\]

whereas empty support gives `r_K<=2E+H-1`. The primitive scale test now
uses a numerical displacement linear in `n`, independently of how many
proper divisors `N` has.

## 4. Recovering the smallest prime factor

At `a_0=3`, the support is nonempty exactly when the odd integer `N`
is composite. If it is empty, the smallest prime factor is `N` itself.

Otherwise, (12) is a monotone suffix query on the odd thresholds:
it is true exactly through the largest proper divisor `d_max`.
That divisor is at most `N/3`. Binary search therefore recovers
`d_max` with `O(n)` coefficient queries. Then

\[
 \operatorname{spf}(N)=N/d_{\max}.
 \tag{15}
\]

This follows because proper divisors correspond to complementary
divisors at least three, and maximizing one minimizes the other.
Even inputs have smallest prime factor two and need no query.

The [conditional wrapper](../tools/leading_bezout_factoring.py) implements
the contract construction, zero test, and threshold search. It receives
the global coefficient oracle explicitly as an argument. Its overhead
is polynomial in `n`, and every oracle input has `O(n)` encoded size
and requests one scalar with `O(n)` bits of precision.

Consequently a polynomial-bit-time evaluator for precisely the scalar
contract (5) would imply polynomial-time deterministic factoring.
More generally, a uniform evaluator at a qualifying improved
subexponential cost would transfer that bound through only polynomially
many queries. The wrapper does not provide such an evaluator.

## 5. The remaining computation

The single coefficient has the exact Lagrange expression

\[
 c(N,a_0)=w\,2^{K+2E}
       \sum_{a\in A}\frac1{P'(a)Q_t(a)}\pmod {2^p}.
 \tag{16}
\]

Each individual row can be normalized and evaluated by the existing
interval product methods. Equation (16) still has `m` rows. The
coefficient contract does not permit that row enumeration at polynomial
cost. Likewise, its short output does not itself provide a short
arithmetic circuit for the normalized sum.

The required aggregation retains the arithmetic correlation between
`a` and its reciprocal representative in `B`. The extra padding makes
the final factor weights simple, but it does not compute that matching.
The [displaced extraction note](displaced-bezout-extraction.md) states
the related full-jet and primitive-scale obligations explicitly.

Thus the remaining engine can be posed narrowly: evaluate (5), with the
specified normalization, on the ordinary progressions (2), at the
linear precision in (1), without expanding the interval. Either this
coefficient engine or an efficient evaluator for the primitive exponent
in (14) would complete a factoring reduction. Neither is supplied here.

## 6. Exact evidence

The [checker](../tools/check_leading_bezout_factoring.py) separates the
conditional wrapper from its test oracle. That oracle is explicitly
given the finite proper-divisor list of each small input. It forms the
actual derivative-unit weights from binomial coefficients and the
normalized factorial unit; it does not replace them by the desired
answer before checking (11).

The [receipt](../evidence/leading-bezout-factoring.json) records 128 odd
inputs, 460 finite known-factor oracle queries, 471 exact relative
factor-weight checks, and 60 independent literal factorial-unit controls.
The smallest-prime outputs agree with direct small-input division.
Even-input controls exercise the zero-query branch.

A separate 402-bit example constructs only the encoded contract: its
precision is 805 bits, its interval count has binary exponent 1207,
and its displacement exponent is 3222. No coefficient or factor is
computed for that example.

```sh
python3 tools/check_leading_bezout_factoring.py --output /tmp/leading-bezout.json
cmp /tmp/leading-bezout.json evidence/leading-bezout-factoring.json
```

The finite oracle tests check the reduction and its normalizations.
The mathematical argument establishes its unrestricted input range.
Neither substitutes for the global coefficient evaluator in section 5.
