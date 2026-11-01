# Verified divisor recovery in a varying characteristic

## Result and scope

A weighted divisor evaluator need not use the same characteristic at every
input length. The existing polynomial-hash and additive-character recovery
method extends to this setting with the tight workspace bound below.

This note gives a uniform deterministic recovery algorithm. It works over
an explicitly represented finite field `K_n`, selected deterministically
from the input length, provided the weighted divisor sums below are
available. The field and all arithmetic costs are part of the algorithm.
It neither assumes a nonzeroness oracle nor treats a field element as a
unit-cost integer of unbounded length.

Let `A(n)` bound one weighted sum, coefficient-algebra setup, and one
scalar operation in bit operations. Let `S(n)` bound one weighted sum's
workspace, that setup workspace, and a constant number of scalar registers.
These costs include the supplied coefficient representation. The least
prime factor, complete factorization, sorted divisors, and exact
divisor-prefix counts can be obtained in

```math
T(n)=A(n)\,2^{O(n/\log(n+2))}\operatorname{poly}(n),\qquad
S_{\rm recovery}(n)=O\bigl(S(n)+n/\log(n+2)\bigr).       \tag{1}
```

The input and any prefix or interval endpoints are read-only; outputs are
streamed. The same conclusion applies to a finite commutative coefficient
algebra when it comes with the character constants and faithful constant
embedding used below. For finite fields those properties are constructed
explicitly, in every characteristic.

Consequently an effective crossing construction over finite fields whose
characteristic and degree vary with the original block length can retain
any supplied logarithmic saving exponent `0<alpha<=1`, once the streamed
evaluator and a sufficiently cheap common field for each input are built.
The global-block construction at the end explains that last requirement.
No crossing family with `alpha>1/2` is supplied here.

The prior weighted multiplication network and virtual-factor interface are
from [the fixed-characteristic construction](https://github.com/1337isnot1337/integer-factoring/blob/d94e4a99fec453768d69ae083038ba65cf31ae9e/research/weighted-prefix-constructor.md)
and [paired-descriptor recovery](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/sublinear-space-half-rate.md).
Polynomial evaluation hashes, binary and ternary coordinate characters,
and the irreducible-polynomial search are already proved in
[bit-tagged witness recovery, Sections 5–6](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/bit-tagged-witness-recovery.md).
The new conclusion combines that construction with one-bit queries and
paired virtual descriptors, charges only a constant number of coefficient
registers, permits the coefficient field to change with the input, and
composes it with the global-block construction. In particular, the
characteristic-two extension and the hash identity themselves are prior
work. The details are restated below to make the tighter space accounting
and varying-field dependencies explicit. This is ordinary mathematics and
exact computation, without new Lean work.

## 1. A computable bound on the number of divisors

Let `N>=2` have binary length `n`. An elementary estimate suffices here;
no prime enumeration or asymptotic formula for the maximal divisor
function is needed.

For `n>=4` set

```math
h=\left\lfloor\frac{\lfloor\log_2 n\rfloor}{2}\right\rfloor,
\quad y=2^h,\quad
q(n)=\min\left(n,\left\lceil\frac n h\right\rceil
                    +y\,\operatorname{bitlength}(n+1)\right).
                                                               \tag{2}
```

For `n<4`, set `q(n)=n`. If `N=prod_p p^{a_p}`, then `a_p<=n`.
There are at most `y` primes at most `y`, so their contribution to
`log2 tau(N)` is at most `y bitlength(n+1)`. For primes greater than
`y`, `sum a_p < n/h`, and `a+1<=2^a` bounds their contribution by
`ceil(n/h)`. Also `tau(N)<=N<2^n`. Therefore

```math
\tau(N)\le D(n):=2^{q(n)},\qquad
q(n)=O(n/\log(n+2)).                                    \tag{3}
```

The second assertion follows from `y<=sqrt(n)` and
`h=Theta(log n)`. Formula (2) uses integer bit lengths and division only.
The bound is deliberately loose; improving its constant does not change
the theorem's exponent.

## 2. A coefficient algebra with short additive characters

Choose the control-field characteristic

```math
\ell=\begin{cases}2,&\operatorname{char}K_n\ne2,\\
                  3,&\operatorname{char}K_n=2.
      \end{cases}                                           \tag{4}
```

When `ell=2`, retain the coefficient field and put `omega=-1`.
When `ell=3`, use the coefficient algebra

```math
R_n=K_n[W]/(W^2+W+1),\qquad \omega=W.                    \tag{5}
```

The latter is free of rank two over `K_n`. It can be a field or a
product algebra; no test of which case holds is needed. Constants embed
faithfully, and

```math
\omega^3=1,\qquad 1+\omega+\omega^2=0.
```

Thus in either case `ell` is a unit, `omega^ell=1`, and
`sum_{j=0}^{ell-1} omega^{aj}` is `ell` for `a=0` and zero for
every nonzero `a in F_ell`. In characteristic two, `ell=3` is the
constant one and is its own inverse. The extension (5) changes scalar
storage and arithmetic costs by a fixed factor under dense arithmetic.
The exact crossing identities and their factors extend from `K_n` to
`R_n` because they use only sums and products after normalization.

Choose the least positive integer `k` such that

```math
Q=\ell^k>nD(n),                                       \tag{6}
```

and construct a control field `H=F_(ell^k)` in a polynomial basis.
Its elements are `k` coordinates over `F_ell`. In particular
`k=O(n/log(n+2))` and `Q=2^{O(n/log(n+2))}`. The control field
is used to generate hashes. It is not embedded in the coefficient algebra.

For control vectors `c,z in F_ell^k`, use their ordinary coordinate dot
product and define

```math
\chi_c(z)=\omega^{c\cdot z}.
```

Summing independently over the `k` coordinates proves

```math
Q^{-1}\sum_{c\in H}\chi_c(z-b)=[z=b]                 \tag{7}
```

in `R_n`. Here `Q^{-1}` means the inverse of the constant `ell^k`,
which exists. No field trace or algebraic-number embedding is required.

### Uniform control-field construction

Enumerate monic degree-`k` polynomials over `F_ell`. Test each by
division by all monic polynomials of degrees from one through
`floor(k/2)`. It is irreducible exactly when none divides it. Stop at
the first irreducible polynomial. Such polynomials exist for every `k`.
For example, the usual factorization of `X^(ell^k)-X` gives
`k I_k >= ell^k-sum_{j<=k/2} ell^j > 0` for `k>=2`, where
`I_k` counts monic irreducibles of degree `k`; degree one is immediate.

This simple search takes `ell^{O(k)} poly(k)` time and `O(k)`
working bits: enumerate candidates and divisors sequentially and perform
schoolbook polynomial division. It fits (1). Faster field construction is
unnecessary for this result. Once constructed, control-field addition,
multiplication, powers, and coordinate dot products use polynomial time
and `O(k)` workspace.

## 3. Hash every divisor without writing it

For an `n`-bit nonnegative integer `d`, including leading zeroes, define

```math
P_d(T)=\sum_{j=0}^{n-1}d_jT^j\in H[T],\qquad
h_a(d)=P_d(a),\quad a\in H.                            \tag{8}
```

Distinct binary integers give distinct polynomials over either control
characteristic. A nonzero difference has degree less than `n`, and hence
at most `n-1` roots. Fix a divisor `d` of `N`. The number of `a` for
which another divisor collides with `d` is at most

```math
(n-1)(\tau(N)-1)<nD(n)<Q.                              \tag{9}
```

Therefore some `a` makes `d` the unique divisor in its bucket
`b=h_a(d)`. The algorithm enumerates all `a,b in H`; it does not have
to recognize an isolating hash beforehand. Every divisor has at least
one isolating descriptor `(a,b)`, stored in `O(k)` bits.

The Fourier phase is bit-separable:

```math
\chi_c(h_a(d))
 =\prod_{j=0}^{n-1}\left(\omega^{c\cdot a^j}\right)^{d_j}.
                                                               \tag{10}
```

This is precisely a product of weights on the first factor's input bits.
The coefficient of bit `j` is generated from `a,c,j` using `O(k)`
control bits and polynomial work. The multiplication network's existing
bit-tag interface can additionally multiply a summand by bit `j` of
either `d` or `N/d`.

## 4. Required weighted queries and exact bucket values

The evaluator interface used here returns

```math
W_t(a,c)=\sum_{d\mid N}\chi_c(h_a(d))\,t(d,N/d),         \tag{11}
```

where `t` is one, one specified first-factor bit, or one specified
second-factor bit. It must compute (11) exactly in `R_n`, with its full
bit cost charged. This is a weighted contraction, not a query that
directly provides a divisor or its hash bucket.

Enumerating `c in H` and applying (7) gives

```math
B_t(a,b)=Q^{-1}\sum_{c\in H}\chi_c(-b)W_t(a,c)
        =\sum_{\substack{d\mid N\\h_a(d)=b}}t(d,N/d).   \tag{12}
```

One bucket query uses `Q` weighted queries, one accumulator, and
`O(k)` control bits. All scalar and control-field operations are charged.
If the bucket contains exactly one divisor, its count is one and its
bit-tag answers are the constant zeroes and ones of the two factors.

An unisolated bucket can have count one in the coefficient characteristic.
It can also have bit-tag answers that resemble bits. Neither event is
accepted as evidence of isolation.

## 5. Verify virtual factors and retain short descriptors

For each descriptor `(a,b)`, query its count. If it is not one, reject
the descriptor. Otherwise regard its two sets of bit-tag answers as
virtual `n`-bit strings. Reject if any answer is neither constant zero
nor constant one. The bit at a requested position is recomputed from
(12), rather than stored in a full candidate integer.

Verify that both virtual factors are positive and that their product is
exactly `N`. Ordinary column multiplication can do this with `O(log n)`
additional bits: for each of the `2n` columns, sum at most `n` products
of queried bits, retain the carry, and compare the resulting bit with
the read-only input or an obligatory zero beyond its length. The carry
and current sum are `O(n)`, so their bit lengths are `O(log n)`.
The check makes polynomially many bit queries and never stores either
candidate integer.

Every accepted candidate is an actual divisor. Conversely an isolating
descriptor from (9) returns the true divisor and cofactor, and passes.
Thus the accepted descriptors represent exactly all divisors, with
possible duplicates. False singleton counts have no effect on correctness.

Scan all `Q^2` descriptors. To obtain the least prime factor, retain a
descriptor for the smallest accepted divisor greater than one. Compare
two candidates from their most significant bits using recomputed queries.
This stores a constant number of `O(k)`-bit descriptors. At the end,
stream the winning divisor's bits. The least divisor greater than one is
prime: a proper factor of a composite divisor would be a smaller such
divisor.

To stream sorted divisors, retain the previous output's descriptor and
repeat the scan, selecting the least accepted value strictly greater than
the previous one. Equal values from different descriptors are skipped by
the comparison. There are at most `D(n)+1` passes, including the final
unsuccessful scan. This avoids a table of
previous outputs or a set of divisor-sized integers. To obtain an exact
prefix or interval count, increment an ordinary binary counter during
this enumeration when the virtual bits meet the read-only condition.
The counter needs only `q(n)+1=O(k)` bits. These conclusions do not rely
on reconstructing an integer from a single modular count.

Complete prime factorization also uses only short descriptors. To test
whether an accepted divisor `p>1` is prime, enumerate pairs of accepted
divisors `a,b>1` and use streamed multiplication to test `ab=p`. A
composite `p` has such a pair, and both members divide `N`. This adds
only a fixed power of `Q` to the query count. For a prime `p`, start with
a descriptor for `q=1` and repeatedly scan for a divisor `r` satisfying
`pq=r`. Retain that descriptor as the new `q`. The number of successful
steps is the exponent of `p` in `N`; the chain stops precisely when its
next power does not divide `N`. Sorted divisor enumeration supplies each
prime once. A constant number of descriptors and an `O(log n)`-bit
exponent counter suffice throughout.

## 6. Full cost

Control-field construction uses `2^{O(k)}` time and `O(k)` space.
There are `Q^2` descriptors, at most `D(n)+1<=Q+1` enumeration passes,
polynomially many bit queries per verification or comparison, and `Q`
weighted queries per bit query. Even charging a fresh verification for
every comparison gives at most `Q^{O(1)} poly(n)` weighted queries.
Finite-field and phase calculations contribute polynomial factors in
`n`, `k`, and the supplied scalar arithmetic cost.

The loops retain only a fixed number of control-field elements,
descriptors, counters, and scalar accumulators in addition to one live
weighted-evaluator invocation. All evaluator calls run sequentially.
This proves (1). In particular, no recursion stack of complete candidate
integers, dense list of buckets, or storage for all `Q` characters is used.

For `0<alpha<=1`, both `n/log(n+2)` and every fixed power of `log n`
are `O(n/(log(n+2))^alpha)`. Hence (1) preserves a supplied time bound
`exp(O(n/(log n)^alpha))` and the matching workspace bound.

## 7. One globally chosen finite field for each input

Suppose a uniform raw crossing family at label length `t` is given over
an explicit finite field `K_t`, with polynomially many pieces, exact
factorizations, nonzero scale, separate widths at most `2t`, and width sum

```math
r_i+c_i\le4\rho t+O(t/\log(t+2)),\qquad
1/2\le\rho<1.                                         \tag{13}
```

Assume complete field data, normalization, piece enumeration, and any
requested factor coordinate take `2^{O(t^h)}` bit time for a fixed
integer `h>=1`, including all representation lengths. The prime
characteristic and extension degree may depend on `t`. Use a polynomial
basis and dense arithmetic, or another explicit representation with the
same charged bound. Normalize before applying the degree-two extension
in (5), so only an inverse in the supplied finite field is needed.

Use the existing streamed grid with its `n^(3/4)` stopping rule. Every
node that performs another grid step has label width at least
`n^(1/4)/8`. The [global-block transfer](global-block-field-transfer.md)
proves the complete accounting for the following choice. For a fixed
sufficiently large integer `B`, set
`m=max(1,iroot_h(floor(ceil(log2(n+2))/B)))` for the whole input.
Every active width `lambda` consists of `floor(lambda/m)` full blocks
and a tail shorter than `m`. Restrict the same `m`-bit family to
zero-padded tail labels. All pieces throughout the computation now use
`K_m`; no maps between different characteristics are required.

Group words of raw piece choices by their multiplicities. The two channel
widths each include the logarithm of the multinomial word count. As in
[the type-grouped transfer](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/rational-simultaneous-transfer.md),
the resulting width sum is

```math
4\rho\lambda+
O\bigl(\lambda/\log(m+2)
       +\lambda\log(m+2)/m+m\bigr).                    \tag{14}
```

The number of type indices has logarithm `polylog(n)`. The tail costs
only `O(m)` width; it uses the same raw field and the same scale-one
crossing. A full-input factorization caps any oversized individual cut
at `2lambda` and can be queried by direct summation of the type words.

Here `m=Theta((log n)^(1/h))`. The relative width error is
`O(1/loglog n)`, and its sum over `O(loglog n)` levels is bounded.
Thus the accumulated multiplicative loss is a constant, which preserves
the exact logarithmic exponent. Type-enumeration costs fit the existing
geometric volume budget. Taking `B` large enough after fixing the raw
family's complexity constant makes its complete cost at most
`O(n^{1/128})`. This pays for field elements, all normalization data,
and multiplication tables; the constant extension (5) changes only
that fixed constant. The shared data and live scalar accumulators fit
below the streamed workspace bound. The referenced global-block proof
uses only effective arithmetic over one nonzero finite ring in its
evaluator steps, so it applies to the algebra (5) as well.

The streamed proof therefore gives the needed weighted evaluator over
`R_n`, with `alpha=log2(1/rho)` and the same simultaneous time/space
bound. Combining it with Sections 1–6 proves the varying-characteristic
transfer through `alpha=1`. This conditional theorem enlarges the class
of admissible constructions. It does not assert a new value of `rho`.

## Verification

The companion checker independently constructs small binary and ternary
control fields, checks character orthogonality in odd characteristic and
in the quadratic characteristic-two extension, compares Fourier bucket
queries with literal divisor buckets, and verifies candidate recovery.
It includes modular-cancellation and false-singleton controls. Direct
divisor enumeration supplies only the finite checker oracle; it is not
the proposed asymptotic weighted evaluator.

Run `python3 tools/check_characteristic_independent_recovery.py` to replay
the checks. The source evaluator is Claude/Anthropic's Appendix B at `16a5ce4`, and
the multiplication reduction is Rudich–Rousseau's at `c91cc3f`, as pinned
in the linked repository notes. The extension and checks here used
substantial OpenAI assistant support.
