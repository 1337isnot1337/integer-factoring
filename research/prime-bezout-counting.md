# Factoring from leading coefficients at prime-adic precision one

## Result and scope

The global inverse-coefficient interface can use **precision one** over
polynomially bounded odd primes. On each matched residue rectangle,
the normalized leading coefficient is exactly the number of integer
factor pairs modulo that prime. It has `O(log n)` output bits, where
`n` is the input bit length. Polynomially many such scalar queries
recover the smallest prime factor.

This sharpens the [one-coefficient reduction](leading-bezout-factoring.md).
The progressions need only extend past `N`; their derivative units are
already constant modulo the chosen prime on every row. No extra
interval padding or bound on the number of divisors is needed.

The [implemented wrapper](../tools/prime_bezout_counting.py) receives
the scalar provider explicitly. The global normalized coefficient
evaluator remains open, and the complete factoring score remains `1/2`.
In particular, low output precision does not supply the required
normalization or compute the correlated reciprocal matching.

## 1. A matched prime-adic rectangle

Let `ell` be an odd prime not dividing `N`. Fix positive integers
`a_0,b_0` with `a_0 b_0=N mod ell`, and let

\[
 m=\ell^u,\quad u\ge0,\qquad
 A=\{a_0+\ell i:0\le i<m\},\quad
 B=\{b_0+\ell j:0\le j<m\}.
 \tag{1}
\]

Every node is an `ell`-adic unit. The common derivative order is

\[
 E=m-1+v_\ell((m-1)!)
   =\frac{\ell(m-1)}{\ell-1}-u.
 \tag{2}
\]

Choose `H>=u+1` so that every nonzero difference `N-ab`, for
`a in A,b in B`, has absolute value less than `ell^H`. The two
extreme products determine such an `H` without listing either interval.
Set `K=H+1`.

Choose a quadratic nonresidue `d mod ell`. Work in the unramified
quadratic extension with valuation ring
`O=Z_ell[omega]`, where `omega^2=d`. Such a `d` can be found by
testing `2,...,ell-1`. This field choice is part of the query: a fixed
quadratic that ramifies for a selected prime would be invalid.

Define

\[
 z=N+\ell^K\omega,\qquad
 P(X)=\prod_{a\in A}(X-a),\qquad
 Q_z(X)=\prod_{b\in B}(z-bX).
 \tag{3}
\]

Let `V_z` be the unique inverse of `Q_z` modulo `P`, with
`deg V_z<m`. It exists because `z` is irrational and the nodes are
rational. The scalar query is

\[
 c_\ell(N;A,B)
 =[X^{m-1}]\,\omega\ell^{K+2E}V_z(X)\pmod\ell.
 \tag{4}
\]

The exponent `K+2E` is represented in binary. A provider must compute
the normalized coefficient before reduction; this definition does not
request an expanded power `ell^(K+2E)`.

## 2. Derivative units are constant on the full grid

For `a_i=a_0+ell i`, put `alpha_i=ell^(-E)P'(a_i)`. Directly,

\[
 P'(a_i)=\ell^{m-1}(-1)^{m-1-i}i!(m-1-i)!,\qquad
 \frac{\alpha_i}{\alpha_0}
       =\frac{(-1)^i}{\binom{m-1}{i}}.
 \tag{5}
\]

The base-`ell` digits of `m-1` are all `ell-1`. Lucas's identity,
or coefficient comparison in
`(1+T)^(ell^u-1)=(1+T^(ell^u))/(1+T)` over `F_ell`, gives

\[
 \binom{m-1}{i}\equiv(-1)^i\pmod\ell.
 \tag{6}
\]

Thus every binomial in (5) is a unit and `alpha_i=alpha_0 mod ell`.
This also proves the constant derivative order in (2).

Since `m` is odd, `alpha_0` is the `ell`-free part of `(m-1)!`.
Remove one factor of `ell` at a time. Its residue is the product,
for `h=1,...,u`, of all units below `ell^h`. Every nonzero residue
appears `ell^(h-1)` times in that level. Pairing field units with their
inverses shows that their product is `-1`; the repetition count is odd.
Consequently

\[
 \boxed{\alpha_i\equiv(-1)^u\pmod\ell
        \quad\text{for every }0\le i<m.}
 \tag{7}
\]

The column derivative units `beta_j` satisfy the same identity,
independently of the first column. In particular, `alpha_i beta_j=1`
modulo `ell`. Also `ell-1` divides `m-1`, so every unit `a` satisfies

\[
 a^{m-1}\equiv1\pmod\ell.
 \tag{8}
\]

## 3. The leading coefficient is the factor count

For each row `a`, exactly one column `b_*` satisfies
`ab_*=N mod ell^(u+1)`. The column interval is a full set of
representatives of its fixed nonzero residue at this modulus.
For every other column,

\[
 v_\ell(N-ab)=v_\ell(b_*-b)<u+1.
 \tag{9}
\]

The product of these off-match differences has valuation `E`.
If `ab_*=N`, the matched factor of `Q_z(a)` is `ell^K omega`.
If it is not zero, its valuation is below `H` and is unchanged by
the displacement in (3).

Lagrange interpolation gives the exact expression

\[
 \omega\ell^{K+2E}V_z(X)
 =\sum_{a\in A}
  \frac{\omega\ell^{K+2E}}{P'(a)Q_z(a)}
                         \frac{P(X)}{X-a}.
 \tag{10}
\]

Every scalar weight in (10) is integral. For a nonfactor row its
valuation is at least `K-H+1>=1`, so it vanishes modulo `ell`.
For an actual pair `ab=N`, the displacement changes every off-match
unit by zero modulo `ell`, because `K>=u+2`. The row weight becomes

\[
 \frac1{a^{m-1}\alpha_i\beta_j}\equiv1\pmod\ell
 \tag{11}
\]

by (7)--(8). Each polynomial `P(X)/(X-a)` is monic of degree `m-1`.
Taking the leading coefficient proves

\[
 \boxed{c_\ell(N;A,B)
       =\#\{(a,b)\in A\times B:ab=N\}\pmod\ell.}
 \tag{12}
\]

In particular, the requested scalar is real modulo `ell`. This
statement holds for the actual specialization `N`; it does not rely
on a generic determinant value or an assumed absence of factor pairs.

### The complete pair and the maximal primitive scale

Write `D` for the factor count in (12), and let the unique reduced
pair satisfy `U_z P+V_z Q_z=1`, with both degrees less than `m`.
Put `L=K+2E` and clear both polynomials by `omega ell^L`.
Every row is congruent to `a_0 mod ell`, so (10) gives the complete
polynomial reduction

\[
 \widehat V\equiv D(X-a_0)^{m-1}\pmod\ell.
 \tag{12a}
\]

The cleared `U` is integral: divide `omega ell^L-widehat V Q_z`
by the monic integral polynomial `P`. Furthermore, Frobenius and
`m=ell^u` give

\[
 P\equiv(X-a_0)^m,\qquad
 Q_z\equiv-b_0(X-a_0)^m\pmod\ell.
\]

Reducing the Bezout identity in the polynomial ring over the residue
field and cancelling its nonzero monic common power yields

\[
 \boxed{\widehat U\equiv b_0D(X-a_0)^{m-1},\qquad
        \widehat V\equiv D(X-a_0)^{m-1}\pmod\ell.}
 \tag{12b}
\]

Let `r` be the smallest nonnegative integer for which there are
integral polynomials `U,V` with `UP+VQ_z=ell^r`. Monic reduction
puts any such pair in the same degree ranges, so uniqueness gives
`r` by the coefficient content of the reduced pair. Since `omega`
is a unit, the known scale satisfies `r<=L`, and (12b) proves

\[
 \boxed{r=L\quad\Longleftrightarrow\quad D\ne0\pmod\ell.}
 \tag{12c}
\]

Thus a provider for the single equality `r=L` can replace each
coefficient query in the reduction below. This conclusion needs no
unknown-content stabilization threshold or factor-count bound. It
does not compute the primitive scale itself. For the suffix decision,
take the logical OR over all selected residue rectangles and primes.
If every equality test is false, each selected prime divides every
rectangle count and hence their total. Their product exceeds the
possible nonzero total, so the support is empty. Equivalently, the
wrapper can receive zero-or-one equality answers: at most `ell-1`
answers are summed for a prime, so that sum is nonzero modulo `ell`
exactly when at least one answer is one.

## 4. Disjoint residue rectangles give suffix counts

Fix an integer threshold `2<=q<=N`, and choose the least power
`m=ell^u>=N`. For each `r=1,...,ell-1`, start the row interval at
the least integer at least `q` congruent to `r mod ell`. Start its
column interval at the least integer at least two congruent to
`N/r mod ell`.

Every proper factor pair with `a>=q` lies in exactly one of these
rectangles. Both coordinates are below `N`, and the chosen interval
lengths contain all such coordinates. Trivial pairs involving one
are excluded by the starting points. Therefore

\[
 D_q:=\#\{a:a\mid N,\ q\le a,\ N/a\ge2\}
   \equiv\sum_{r=1}^{\ell-1}c_\ell(N;A_r,B_r)\pmod\ell.
 \tag{13}
\]

Each query in this sum requests one coefficient at prime-adic
precision one. No factor weights or factorial units need to be
computed by the suffix decoder.

## 5. An elementary polynomial bound for the prime list

Let `n` be the bit length of `N>=3`, so `n>=2`, and put
`B=16n^2`, `M=B/2`. Write `L=lcm(1,...,2M)` and let `R` be
the product of all primes at most `2M`.

For each prime, Legendre's formula for `binom(2M,M)` is a sum of
zeros and ones, with at most `floor(log_p(2M))` summands. Hence
this binomial divides `L`. Also

\[
 \binom{2M}{M}
 =\prod_{i=1}^M\frac{M+i}{i}\ge2^M,
 \qquad L\le R^{\log_2(2M)}.
 \tag{14}
\]

It follows that the product `R_odd=R/2` of odd primes satisfies

\[
 \log_2 R_{\rm odd}
 \ge\frac{M}{\log_2(2M)}-1
 =\frac{8n^2}{\log_2(16n^2)}-1
 \ge2n-1>n.
 \tag{15}
\]

The last inequality uses `log_2(16n^2)<=4n` for `n>=2`.
Thus the odd primes at most `16n^2` have product greater than `N`.
No prime distribution estimate is required.

Generate these primes in increasing order, stopping once their product
exceeds `N`. Trial-divide `N` by each selected prime first. A divisor
found during this increasing search is its smallest prime factor.
If none is found, every selected prime is a valid field for (12).

Now `0<=D_q<N`. If every residue in (13) vanishes, the prime product
divides `D_q` and exceeds it, so `D_q=0`. Conversely, zero support
gives zero in every field. This proves an exact nonempty-suffix test.
Binary search finds the largest proper divisor; dividing `N` by it
returns the smallest prime factor. Prime and even inputs are handled
by the same empty-support and direct-division branches as before.

## 6. Uniform cost contract

The number of selected primes, their magnitudes, and the number of
residue rectangles are polynomial in `n`. There are `O(n)` threshold
steps. Trial primality tests, finding a quadratic nonresidue, binary
contract construction, and modular summation all have polynomial bit
cost under these bounds.

For every query, `m<ell N`, so its degree and clearing exponent have
`O(n+log ell)` binary digits. The numerical height and displacement
exponents are `O(1+log_ell N)`. The output is one element of `F_ell`,
requiring only `O(log n)` bits.

An evaluator must be uniform in `ell`, the interval starts, `u`, and
`N`, and must account for the normalized coefficient in (4). A theorem
for each fixed prime with uncontrolled prime-dependent constants
would not give a uniform factoring bound. A provider whose total cost
is polynomial in `n` and `ell` would give deterministic polynomial-time
factoring through this reduction. No such provider is supplied.

The exact scalar still equals a global Lagrange sum over `m` rows.
Its low residue precision removes a precision obstacle; it does not
bound the width of a reciprocal-matching tree or evaluate that sum.

## 7. Exact controls

The [checker](../tools/check_prime_bezout_counting.py) evaluates complete
small rational Lagrange sums in `Q(sqrt(d))`, with the same selected
nonresidue as the query contract. These controls test normalization,
real output, derivative units on every row, and the factor-count value.
Additional displacement comparisons check that the retained coefficient
does not depend on extra padding.

The receipt records 75 complete rational leading coefficients, 966
literal derivative-unit checks, eight independent padding comparisons,
and 2,136 elementary primorial, factorial-unit, and Lucas controls.
Eighteen complete rational Bezout pairs also verify (12b), their
exact primitive scales, and the maximum-scale test (12c).

Separate wrapper controls use an explicitly supplied finite divisor
list as their oracle. They test residue partitioning, prime selection,
zero detection, and threshold recovery; they are not a large-domain
coefficient implementation. A 402-bit contract control only checks
encoded parameter sizes. The complete receipt is
[prime-bezout-counting.json](../evidence/prime-bezout-counting.json).
The wrapper runs cover 261 inputs, of which 60 use its supplied oracle,
with 2,096 explicit known-factor coefficient queries in total.

```sh
python3 tools/check_prime_bezout_counting.py --output /tmp/prime-bezout-counting.json
cmp /tmp/prime-bezout-counting.json evidence/prime-bezout-counting.json
```
